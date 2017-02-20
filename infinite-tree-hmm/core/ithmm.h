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
		_add_customer_to_vertical_crp(alpha, node_on_cluster, node_on_cluster->_identifier);
		_add_customer_to_horizontal_crp(_gamma, node_on_cluster, node_on_cluster->_identifier);
	}
	void add_clustering_customer_to_node(Node* node_on_cluster){
		assert(node_on_cluster->_htssb_owner_id == 0);
		double alpha = _alpha * pow(_alpha, node_on_cluster->_depth_v);
		bool new_table_generated = false;
		node_on_cluster->add_customer_to_vertical_crp(alpha, new_table_generated);
	}
	void _add_customer_to_vertical_crp(double alpha, Node* iterator_on_cluster, int target_id){
		TSSB* transition_tssb = iterator_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool new_table_generated = false;
		target_on_htssb->add_customer_to_vertical_crp(alpha, new_table_generated);
		if(new_table_generated && iterator_on_cluster->_parent != NULL){
			_add_customer_to_vertical_crp(alpha, iterator_on_cluster->_parent, target_id);
		}
	}
	void _add_customer_to_horizontal_crp(double gamma, Node* iterator_on_cluster, int target_id){
		TSSB* transition_tssb = iterator_on_cluster->_transition_tssb;
		assert(transition_tssb != NULL);
		Node* target_on_htssb = transition_tssb->find_node_with_id(target_id);
		assert(target_on_htssb != NULL);
		bool new_table_generated = false;
		target_on_htssb->add_customer_to_horizontal_crp(gamma, new_table_generated);
		if(new_table_generated && iterator_on_cluster->_parent != NULL){
			_add_customer_to_horizontal_crp(gamma, iterator_on_cluster->_parent, target_id);
		}
	}
	void remove_customer_from_node(Node* node_on_cluster, bool remove_proxy_customer = false){
		node_on_cluster->remove_customer_from_vertical_crp(remove_proxy_customer);
		node_on_cluster->remove_customer_from_horizontal_crp(remove_proxy_customer);
	}
};

#endif