#ifndef _ithmm_
#define _ithmm_
#include <boost/format.hpp>
#include <cmath>
#include "tssb.h"

class iTHMM{
public:
	TSSB* _clustering_tssb;
	double _alpha;
	double _gamma;
	double _lambda;
	iTHMM(){
		_alpha = 1;
		_gamma = 1;
		_lambda = 1;
		_clustering_tssb = new TSSB(_alpha, _gamma, _lambda);
		Node* root_on_cluster = _clustering_tssb->_root;
		Node* root_on_htssb = new Node(NULL, root_on_cluster->_identifier);
		root_on_htssb->_htssb_owner_id = root_on_cluster->_identifier;
		root_on_cluster->_transition_tssb = new TSSB(root_on_htssb, _alpha, _gamma, _lambda);
	}
	// クラスタリング用TSSBで子ノードを生成した瞬間全てのHTSSBの同じ位置に子ノードを生成する
	Node* generate_child_node(Node* parent_on_cluster){
		assert(parent_on_cluster != NULL);
		// まずクラスタリング用TSSBで子ノードを作る
		Node* child_on_cluster = parent_on_cluster->generate_child();
		// クラスタリング用TSSBの全ノードを収集
		vector<Node*> nodes;
		_clustering_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(auto node_on_cluster: nodes){
			if(node_on_cluster->_transition_tssb == NULL){
				node_on_cluster->_transition_tssb = _clustering_tssb->copy(node_on_cluster->_identifier);
				node_on_cluster->_transition_tssb_myself = node_on_cluster->_transition_tssb->find_node_with_id(node_on_cluster->_identifier);
				assert(node_on_cluster->_transition_tssb_myself != NULL);
			}else{
				// 遷移確率用TSSBでの同じ位置に子ノードを挿入
				Node* parent_on_htssb = node_on_cluster->_transition_tssb->find_node_with_id(parent_on_cluster->_identifier);
				assert(parent_on_htssb != NULL);
				Node* child_on_htssb = new Node(parent_on_htssb, child_on_cluster->_identifier);
				child_on_htssb->_htssb_owner_id = node_on_cluster->_identifier;
				parent_on_htssb->add_child(child_on_htssb);
			}
		}
		return child_on_cluster;
	}
	Node* sample_node(){
		return _stop_node(_clustering_tssb->_root);
	}
	// 止まるノードを決定する
	Node* _stop_node(Node* node){
		assert(node != NULL);
		double alpha = _alpha * pow(_lambda, node->_depth_v);
		double head = node->compute_expectation_of_clustering_vertical_sbr_ratio(alpha);
		node->_children_stick_length = 1 - node->_stick_length * head;
		double bernoulli = Sampler::uniform(0, 1);
		if(bernoulli <= head){			// 表が出たらこのノードに降りる
			return node;
		}
		// 子ノードがある場合
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			assert(child != NULL);
			double head = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _stop_node(child);
			}
		}
		// ない場合生成しながらコインを投げる
		while(true){
			Node* child = generate_child_node(node);
			double head = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _stop_node(child);
			}
		}
	}
	// [0, 1)の一様分布からノードをサンプリング
	Node* retrospective_sampling(double uniform){
		Node* root = _clustering_tssb->_root;
		double ratio = root->compute_expectation_of_clustering_vertical_sbr_ratio(_alpha);
		double sum_probability = ratio;
		root->_children_stick_length = 1.0 - ratio;
		return _retrospective_sampling(uniform, sum_probability, root);
	}
	Node* _retrospective_sampling(double uniform, double &sum_probability, Node* node){
		if(uniform <= sum_probability){
			return node;
		}
		// 棒の長さとノードの確率の関係に気をつける
		// [<------------------- 棒の長さ -------------------]
		// [<--親ノードの確率--><---子ノードに割り当てる長さ --->]
		//					  [<---子1の棒---><---子2の棒--->]
		assert(node->_children_stick_length > 0);
		double stick_length = node->_children_stick_length;	// 子ノードに割り当てる棒の長さの総和
		double sum_stick_length_over_children = 0;			// 子ノードを走査する時の走査済みの棒の長さ
		Node* last_node = NULL;
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			double ratio_h = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
			double alpha = _alpha * pow(_lambda, child->_depth_v);
			double ratio_v = child->compute_expectation_of_clustering_vertical_sbr_ratio(alpha);
			if(child->_identifier == 118){
				// cout << (boost::format("id=%d, uniform=%f, sum=%f, sum_over=%f, len=%f, len + sum=%f, len * ratio_h=%f, sum + len * ratio_h=%f") % child->_identifier % uniform % sum_probability % sum_stick_length_over_children % child->_stick_length % (sum_probability + child->_stick_length) % (stick_length * ratio_h) % (sum_probability + stick_length * ratio_h)).str() << endl;
			}
			if(child->_identifier == 119){
				// cout << (boost::format("id=%d, uniform=%f, sum=%f, sum_over=%f, len=%f, len + sum=%f, len * ratio_h=%f, sum + len * ratio_h=%f") % child->_identifier % uniform % sum_probability % sum_stick_length_over_children % child->_stick_length % (sum_probability + child->_stick_length) % (stick_length * ratio_h) % (sum_probability + stick_length * ratio_h)).str() << endl;
			}
			if(child->_identifier == 120){
				// exit(0);
			}
			if(uniform <= sum_probability + sum_stick_length_over_children + stick_length * ratio_h){
				// stick_length * ratio_hだけだとこのノードの棒の長さなのでratio_vも掛けてこのノードで止まる確率にする必要がある
				sum_probability += sum_stick_length_over_children + stick_length * ratio_h * ratio_v;
				if(uniform <= sum_probability){
					return child;
				}
				if(child->has_child()){
					return _retrospective_sampling(uniform, sum_probability, child);
				}
				// 生成する
				// cout << "will be " << child->_identifier << "'s child." << endl;
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				//
				return child;
			}
			sum_stick_length_over_children += child->_stick_length;
			stick_length *= 1.0 - ratio_h;
			last_node = child;
		}
		// 見つからなかったら一番右端のノードを返す
		if(last_node != NULL){
			if(last_node->has_child()){
				return _retrospective_sampling(uniform, sum_probability, last_node);
			}
			return last_node;
		}
		return NULL;
	}
	void add_htssb_customer_to_node(Node* node_on_cluster){
		assert(node_on_cluster->_htssb_owner_id == 0);
		double alpha = _alpha * pow(_alpha, node_on_cluster->_depth_v);
		_add_htssb_customer_to_vertical_crp(alpha, node_on_cluster, node_on_cluster->_identifier);
		_add_htssb_customer_to_horizontal_crp(_gamma, node_on_cluster, node_on_cluster->_identifier);
	}
	void add_clustering_customer_to_node(Node* node_on_cluster){
		assert(node_on_cluster->_htssb_owner_id == 0);
		double alpha = _alpha * pow(_alpha, node_on_cluster->_depth_v);
		bool new_table_generated = false;
		node_on_cluster->add_customer_to_vertical_crp(alpha, new_table_generated);
		// node_on_cluster->add_customer_to_horizontal_crp(alpha, new_table_generated);
	}
	void _add_htssb_customer_to_vertical_crp(double alpha, Node* target_on_cluster, int target_id){
		TSSB* transition_tssb = target_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool new_table_generated = false;
		target_on_htssb->add_customer_to_vertical_crp(alpha, new_table_generated);
		if(new_table_generated && target_on_cluster->_parent != NULL){
			_add_htssb_customer_to_vertical_crp(alpha, target_on_cluster->_parent, target_id);
		}
	}
	void _add_htssb_customer_to_horizontal_crp(double gamma, Node* target_on_cluster, int target_id){
		TSSB* transition_tssb = target_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool new_table_generated = false;
		target_on_htssb->add_customer_to_horizontal_crp(gamma, new_table_generated);
		if(new_table_generated && target_on_cluster->_parent != NULL){
			_add_htssb_customer_to_horizontal_crp(gamma, target_on_cluster->_parent, target_id);
		}
	}
	void remove_htssb_customer_from_node(Node* node_on_cluster){
		assert(node_on_cluster->_htssb_owner_id == 0);
		bool empty_table_deleted = false;
		_remove_customer_from_vertical_crp(node_on_cluster, node_on_cluster->_identifier);
		_remove_customer_from_horizontal_crp(node_on_cluster, node_on_cluster->_identifier);
	}
	void _remove_customer_from_vertical_crp(Node* target_on_cluster, int target_id){
		TSSB* transition_tssb = target_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool empty_table_deleted = false;
		target_on_htssb->remove_customer_from_vertical_crp(empty_table_deleted);
		if(empty_table_deleted && target_on_cluster->_parent != NULL){
			_remove_customer_from_vertical_crp(target_on_cluster->_parent, target_id);
		}
	}
	void _remove_customer_from_horizontal_crp(Node* target_on_cluster, int target_id){
		TSSB* transition_tssb = target_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool empty_table_deleted = false;
		target_on_htssb->remove_customer_from_horizontal_crp(empty_table_deleted);
		if(empty_table_deleted && target_on_cluster->_parent != NULL){
			_remove_customer_from_horizontal_crp(target_on_cluster->_parent, target_id);
		}
	}
	// クラスタリング用TSSBから客が消える場合、木構造が変化する可能性があるので専用メソッドを用意
	void remove_clustering_customer_from_node(Node* node_on_cluster){
		assert(node_on_cluster->_htssb_owner_id == 0);
		_remove_clustering_customer_from_vertical_crp_on_node(node_on_cluster);
	}
	// 客を除去
	void _remove_clustering_customer_from_vertical_crp_on_node(Node* target_on_cluster){
		target_on_cluster->dump();
		// cout << "remove_customer_from_vertical_crp: " << tssb_identifier << ", " << node->_identifier << endl;
		Table* table = target_on_cluster->get_vertical_table();
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		target_on_cluster->decrement_vertical_stop_count();
		_decrement_clustering_vertical_pass_counts_on_node(target_on_cluster->_parent);
	}
	void _decrement_clustering_vertical_pass_counts_on_node(Node* parent_on_cluster){
		if(parent_on_cluster == NULL){
			return;
		}
		parent_on_cluster->decrement_vertical_pass_count();
		delete_invalid_children(parent_on_cluster);
		_decrement_clustering_vertical_pass_counts_on_node(parent_on_cluster->_parent);
	}
	// void remove_clustering_customer_from_horizontal_crp_on_node(Node* node){
	// 	// cout << "remove_customer_from_horizontal_crp: " << _identifier << "," << node->_identifier << endl;
	// 	Table* table = get_horizontal_table();
	// 	assert(table != NULL);
	// 	table->remove_customer(empty_table_deleted);
	// 	// 停止回数・通過回数を更新
	// 	Node* stopped_child = this;
	// 	Node* parent = _parent;
	// 	assert(parent);
	// 	while(parent){
	// 		bool found = false;
	// 		for(int i = parent->_children.size() - 1;i >= 0;i--){	// 逆向きに辿らないと通過ノードが先に消えてしまう
	// 			Node* child = parent->_children[i];
	// 			// cout << "foreach: " << child->_identifier << endl;
	// 			// cout << "foreach: " << child->_identifier << endl;

	// 			if(child == stopped_child){
	// 				found = true;
	// 				child->decrement_horizontal_stop_count();
	// 				if(delete_node_if_empty){
	// 					child->delete_node_on_clustering_tssb_if_needed();
	// 				}
	// 				continue;
	// 			}
	// 			if(found){
	// 				child->decrement_horizontal_pass_count();
	// 				child->delete_node_on_clustering_tssb_if_needed();
	// 			}
	// 		}
	// 		stopped_child = parent;
	// 		parent = parent->_parent;
	// 	}
	// 	stopped_child->decrement_horizontal_stop_count();
	// 	if(delete_node_if_empty){
	// 		stopped_child->delete_node_on_clustering_tssb_if_needed();
	// 	}
	// }
	void delete_invalid_children(Node* parent){
		vector<Node*> &children = parent->_children;
		for(int i = children.size() - 1;i >= 0;i--){
			Node* child = children[i];
			bool success = delete_node_on_clustering_tssb_if_needed(child);
			if(success == false){	// 失敗したらそれ以上は消さない
				break;
			}
		}
	}
	bool delete_node_on_clustering_tssb_if_needed(Node* node_on_cluster){
		if(node_on_cluster->_depth_v == 0){
			return false;
		}
		cout << "deleting ... ";
		node_on_cluster->dump();
		assert(node_on_cluster->_parent != NULL);
		int delete_id = node_on_cluster->_identifier;
		int parent_id = node_on_cluster->_parent->_identifier;
		int pass_count_v = node_on_cluster->get_vertical_pass_count();
		int stop_count_v = node_on_cluster->get_vertical_stop_count();
		int pass_count_h = node_on_cluster->get_horizontal_pass_count();
		int stop_count_h = node_on_cluster->get_horizontal_stop_count();
		if(pass_count_v != 0){
			return false;
		}
		if(stop_count_v != 0){
			return false;
		}
		if(pass_count_h != 0){
			return false;
		}
		if(stop_count_h != 0){
			return false;
		}
		Node* delete_node = node_on_cluster->_parent->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		// クラスタリング用TSSBの全ノードを収集
		vector<Node*> nodes;
		_clustering_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(auto node_on_cluster: nodes){
			assert(node_on_cluster->_transition_tssb != NULL);
			// 遷移確率用TSSBでの同じ位置の子ノードを削除
			Node* parent_on_htssb = node_on_cluster->_transition_tssb->find_node_with_id(parent_id);
			assert(parent_on_htssb != NULL);
			assert(parent_on_htssb->_htssb_owner_id != 0);
			Node* delete_node = parent_on_htssb->delete_child_node(delete_id);
			if(delete_node != NULL){
				TSSB* delete_tssb = delete_node->_transition_tssb;
				delete delete_node;
				delete delete_tssb;
			}
		}
		return true;
	}
};

#endif