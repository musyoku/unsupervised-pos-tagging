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
		node_on_cluster->add_customer_to_horizontal_crp(alpha, new_table_generated);
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
		_remove_clustering_customer_from_horizontal_crp_on_node(node_on_cluster);
	}
	// 客を除去
	void _remove_clustering_customer_from_vertical_crp_on_node(Node* target_on_cluster){
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
	void _remove_clustering_customer_from_horizontal_crp_on_node(Node* target_on_cluster){
		// cout << "remove_customer_from_horizontal_crp: " << tssb_identifier << ", " << node->_identifier << endl;
		Table* table = target_on_cluster->get_horizontal_table();
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		// 通過回数・停止回数を減らす
		Node* stopped_child = target_on_cluster;
		Node* parent = target_on_cluster->_parent;
		while(parent){
			bool found = false;
			for(int i = parent->_children.size() - 1;i >= 0;i--){	// 逆向きに辿らないと通過ノードが先に消えてしまう
				Node* child = parent->_children[i];
				if(child == stopped_child){
					found = true;
					child->decrement_horizontal_stop_count();
					continue;
				}
				if(found){
					child->decrement_horizontal_pass_count();
				}
			}
			delete_invalid_children(parent);
			stopped_child = parent;
			parent = parent->_parent;
		}
		// ルートノードのカウントを減らす
		stopped_child->decrement_horizontal_stop_count();
	}
	// 縦の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double compute_expectation_of_htssb_vertical_sbr_ratio_on_node(Node* target_on_cluster){
		// cout << "target" << endl << "	";
		// target_on_cluster->dump();
		double sbr_ratio = 0;
		int depth_v = target_on_cluster->_depth_v;
		int num_parents = depth_v;
		// 遷移確率用TSSBのトップレベルノードから対象ノードまでのノードの水平方向のインデックスを格納
		// ルートノードから子ノードをたどって到達できるようにするには水平方向の位置が分かればよい
		int* node_horizontal_indices = target_on_cluster->_node_indices;
		Node* parent_on_cluster = target_on_cluster;
		node_horizontal_indices[depth_v - 1] = target_on_cluster->_depth_h;
		for(int n = 0;n < num_parents - 1;n++){
			parent_on_cluster = parent_on_cluster->_parent;
			node_horizontal_indices[depth_v - n - 2] = parent_on_cluster->_depth_h;
		}

		// cout << "horizontal indices:" << endl;
		// for(int n = 0;n < num_parents;n++){
		// 	cout << node_horizontal_indices[n] << " -> ";
		// }
		// cout << endl;

		// クラスタリング用TSSBでの基準となるノードを選ぶ
		Node** pointer_nodes_on_cluster_top_to_bottom = target_on_cluster->_pointer_nodes_v;
		Node* pointer_on_cluster = target_on_cluster;
		pointer_nodes_on_cluster_top_to_bottom[depth_v] = pointer_on_cluster;
		for(int n = 0;n < num_parents;n++){
			pointer_on_cluster = pointer_on_cluster->_parent;
			assert(pointer_on_cluster != NULL);
			pointer_nodes_on_cluster_top_to_bottom[depth_v - n - 1] = pointer_on_cluster;
		}

		// ノードの深さをdとすると(d+1)^2回計算する必要がある
		double* stop_ratio_over_parent = target_on_cluster->_stop_ratio_v_over_parent;
		double* stop_probability_over_parent = target_on_cluster->_stop_probability_v_over_parent;
		for(int n = 0;n < num_parents + 1;n++){
			// cout << "n = " << n << endl;
			// クラスタリング用TSSBでの基準となるノードを選ぶ
			// nが増えるごとに下に降りていく
			Node* pointer_on_cluster = pointer_nodes_on_cluster_top_to_bottom[n];
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			Node* pointer_on_htssb = pointer_on_cluster->_transition_tssb->_root;	// 遷移確率用TSSBのルートから始める
			assert(pointer_on_htssb);
			for(int m = 0;m < num_parents + 1;m++){
				// cout << "m = " << m << endl;
				// pointer_on_htssb->dump();
				// cout << "pointer_on_htssb = " << pointer_on_htssb->_identifier << endl;
				if(pointer_on_cluster->_depth_v == 0){	// クラスタリング用TSSBの親ノードの場合
					int pass_count = pointer_on_htssb->get_vertical_pass_count();
					int stop_count = pointer_on_htssb->get_vertical_stop_count();
					double ratio_v = (1.0 + stop_count) / (1.0 + _alpha + stop_count + pass_count);
					// cout << "ratio_v = " << ratio_v << endl;
					stop_ratio_over_parent[m] = ratio_v;
				}else{	// 親の遷移確率用HTSSBから生成
					int pass_count = pointer_on_htssb->get_vertical_pass_count();
					int stop_count = pointer_on_htssb->get_vertical_stop_count();
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					double alpha = _alpha * pow(_lambda, pointer_on_htssb->_depth_v);
					// cout << "alpha = " << alpha << endl;
					double ratio_v = (alpha * parent_stop_probability + stop_count) / (alpha * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
					// cout << "ratio_v = " << ratio_v << endl;
					stop_ratio_over_parent[m] = ratio_v;
					sum_parent_stop_probability += parent_stop_probability;
					sbr_ratio = ratio_v;
				}
				// 親から子へ降りていく
				// 水平方向の位置が分ればアクセス可能
				if(m < num_parents){
					// cout << "next index: " << node_horizontal_indices[m] << endl;
					assert(node_horizontal_indices[m] < pointer_on_htssb->_children.size());
					pointer_on_htssb = pointer_on_htssb->_children[node_horizontal_indices[m]];
					assert(pointer_on_htssb != NULL);
				}
			}
			// 計算した棒を折る比率から確率を計算
			if(n < num_parents){
				double rest_stick_length = 1;
				for(int m = 0;m < num_parents + 1;m++){
					double ratio_v = stop_ratio_over_parent[m];
					double stop_probability = rest_stick_length * ratio_v;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_v;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					stop_probability_over_parent[m] = stop_probability;
				}
			}
		}
		return sbr_ratio;
	}
	// 横の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double compute_expectation_of_htssb_horizontal_sbr_ratio_on_node(Node* target_on_cluster){
		int depth_v = target_on_cluster->_depth_v;
		int depth_h = target_on_cluster->_depth_h;
		if(depth_v == 0){	// ルートノードなら必ず止まる
			return 1;
		}
		double sbr_ratio = 0;
		int num_parents = depth_v;
		// 遷移確率用TSSBのトップレベルノードから対象ノードまでのノードの水平方向のインデックスを格納
		// ルートノードから子ノードをたどって到達できるようにするには水平方向の位置が分かればよい
		int* node_horizontal_indices = target_on_cluster->_node_indices;
		Node* parent_on_cluster = target_on_cluster;
		node_horizontal_indices[depth_v - 1] = target_on_cluster->_depth_h;
		for(int n = 0;n < num_parents - 1;n++){
			parent_on_cluster = parent_on_cluster->_parent;
			node_horizontal_indices[depth_v - n - 2] = parent_on_cluster->_depth_h;
		}
		// クラスタリング用TSSBのルートノード
		Node* root_on_cluster = parent_on_cluster->_parent;

		// クラスタリング用TSSBで辿れる全ての親ノードが持つ遷移確率用TSSB上での対象ノードを持っている親ノードへのポインタ
		Node** parents_on_htssb_contain_target = target_on_cluster->_pointer_nodes_v;
		Node* pointer_on_cluster = target_on_cluster;
		Node* parent_contains_target_on_htssb = pointer_on_cluster->_transition_tssb_myself->_parent;
		assert(parent_contains_target_on_htssb);
		parents_on_htssb_contain_target[depth_v] = parent_contains_target_on_htssb;
		// cout << depth_v << " <- " << endl << "	";
		// parent_contains_target_on_htssb->dump();
		for(int n = 0;n < num_parents;n++){
			pointer_on_cluster = pointer_on_cluster->_parent;
			Node* pointer_on_htssb = pointer_on_cluster->_transition_tssb->_root;
			// cout << "searching ..." << endl;
			// クラスタリング用TSSBのそれぞれのノードが持つ遷移確率用TSSBのルートノードからたどっていき目的のノードの親ノードを見つける
			for(int m = 0;m < num_parents - 1;m++){
				pointer_on_htssb = pointer_on_htssb->_children[node_horizontal_indices[m]];
				// cout << "	";
				// pointer_on_htssb->dump();
			}
			assert(pointer_on_htssb != NULL);
			parents_on_htssb_contain_target[depth_v - n - 1] = pointer_on_htssb;
			// cout << depth_v - n - 1 << " <- " << endl << "	";
			// pointer_on_htssb->dump();
		}

		double* stop_ratio_over_parent = target_on_cluster->_stop_ratio_h_over_parent;
		double* stop_probability_over_parent = target_on_cluster->_stop_probability_h_over_parent;
		pointer_on_cluster = root_on_cluster;
		for(int n = 0;n < num_parents + 1;n++){
			// クラスタリング用TSSBでの基準となるノードを選ぶ
			// nが増えるごとに下に降りていく
			// cout << "pointer" << endl << "	";
			// pointer_on_cluster->dump();
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			parent_contains_target_on_htssb = parents_on_htssb_contain_target[n];
			// cout << "parent contains" << endl << "	";
			// parent_contains_target_on_htssb->dump();
			assert(parent_contains_target_on_htssb);
			for(int m = 0;m < depth_h + 1;m++){
				Node* child_on_htssb = parent_contains_target_on_htssb->_children[m];
				// cout << "m = " << m << endl << "	";
				// child_on_htssb->dump();
				if(pointer_on_cluster->_depth_v == 0){	// 親ノードの場合
					int pass_count = child_on_htssb->get_horizontal_pass_count();
					int stop_count = child_on_htssb->get_horizontal_stop_count();
					double ratio_h = (1.0 + stop_count) / (1.0 + _gamma + stop_count + pass_count);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
				}else{
					int pass_count = child_on_htssb->get_horizontal_pass_count();
					int stop_count = child_on_htssb->get_horizontal_stop_count();
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					double ratio_h = (_gamma * parent_stop_probability + stop_count) / (_gamma * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
					sum_parent_stop_probability += parent_stop_probability;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					sbr_ratio = ratio_h;
				}

			}
			// 親から子へ降りていく
			// 水平方向の位置が分ればアクセス可能
			if(n < num_parents){
				// cout << "next index: " << node_horizontal_indices[n] << endl;
				assert(node_horizontal_indices[n] < pointer_on_cluster->_children.size());
				pointer_on_cluster = pointer_on_cluster->_children[node_horizontal_indices[n]];
				assert(pointer_on_cluster != NULL);
			}
			if(n < num_parents){
				double rest_stick_length = 1;
				for(int m = 0;m < depth_h + 1;m++){
					double ratio_h = stop_ratio_over_parent[m];
					double stop_probability = rest_stick_length * ratio_h;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_h;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					stop_probability_over_parent[m] = stop_probability;
				}
			}
		}
		return sbr_ratio;
	}
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