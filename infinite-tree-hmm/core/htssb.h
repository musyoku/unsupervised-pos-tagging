#ifndef _htssb_
#define _htssb_
#include <boost/format.hpp>
#include <cmath>
#include "node.h"

class HTSSB{
public:
	Node* _root;
	double _alpha;
	double _gamma;
	double _lambda;
	HTSSB(){
		_root = new Node(NULL);
		_root->_stick_length = 1;
		_alpha = 1;
		_gamma = 1;
		_lambda = 1;
	}
	void add_customer_to_clustering_node(Node* node){
		double alpha = _alpha * pow(_alpha, node->_depth_v);
		node->add_customer_to_clustering_vertical_crp(alpha);
		node->add_customer_to_clustering_horizontal_crp(_gamma);
	}
	void remove_customer_from_clustering_node(Node* node){
		node->remove_customer_from_clustering_vertical_crp();
		node->remove_customer_from_clustering_horizontal_crp();
	}
	void add_customer_to_htssb_node(Node* node){
		double alpha = _alpha * pow(_alpha, node->_depth_v);
		node->add_customer_to_htssb_vertical_crp(alpha, node->_identifier, node);
		node->add_customer_to_htssb_horizontal_crp(_gamma, node->_identifier, node);
	}
	void remove_customer_from_htssb_node(Node* node){
		// node->remove_customer_from_htssb_vertical_crp();
		// node->remove_customer_from_htssb_horizontal_crp();
	}
	// [0, 1)の一様分布からノードをサンプリング
	Node* retrospective_sampling(double uniform){
		double ratio = _root->compute_expectation_of_clustering_vertical_sbr_ratio(_alpha);
		double sum_probability = ratio;
		_root->_children_stick_length = 1.0 - ratio;
		return _retrospective_sampling(uniform, sum_probability, _root);
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
	void update_stick_length(){
		double ratio_v = _root->compute_expectation_of_clustering_vertical_sbr_ratio(_gamma);
		double sum_probability = ratio_v;
		_root->_stick_length = 1;
		_root->_children_stick_length = 1.0 - ratio_v;
		_root->_probability = ratio_v;
		_update_stick_length(sum_probability, _root);
		cout << "sum_probability: " << sum_probability << endl;
	}
	void _update_stick_length(double &sum_probability, Node* node){
		assert(node->_children_stick_length > 0);
		double rest_stick_length = node->_children_stick_length;
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			double ratio_h = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
			child->_stick_length = rest_stick_length * ratio_h;
			double ratio_v = child->compute_expectation_of_clustering_vertical_sbr_ratio(_gamma);
			child->_probability = child->_stick_length * ratio_v;
			sum_probability += child->_probability;
			rest_stick_length *= 1.0 - ratio_h;
			double alpha = _alpha * pow(_lambda, child->_depth_v);
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
			if(child->has_child()){
				_update_stick_length(sum_probability, child);
			}
		}
	}
	Node* sample_node(){
		return _stop_node(_root);
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
			Node* child = node->generate_child();
			double head = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _stop_node(child);
			}
		}
	}
	void enumerate_nodes_from_left_to_right(vector<Node*> &nodes){
		_enumerate_nodes_from_left_to_right(_root, nodes);
	}
	void _enumerate_nodes_from_left_to_right(Node* node, vector<Node*> &nodes){
		nodes.push_back(node);
		for(const auto child: node->_children){
			_enumerate_nodes_from_left_to_right(child, nodes);
		}
	}
	int get_num_nodes(){
		return _get_num_children(_root) + 1;
	}
	int _get_num_children(Node* node){
		int sum = node->_children.size();
		for(const auto child: node->_children){
			sum += _get_num_children(child);
		}
		return sum;
	}
	int get_max_depth(){
		return _get_max_depth(_root);
	}
	int _get_max_depth(Node* node){
		int max_depth = node->_depth_v;
		for(const auto child: node->_children){
			int depth = _get_max_depth(child);
			if(depth > max_depth){
				max_depth = depth;
			}
		}
		return max_depth;
	}
	void dump_tssb(int identifier){
		_dump_tssb(_root, identifier);
	}
	void _dump_tssb(Node* node, int identifier){
		string tab = "";
		for(int i = 0;i < node->_depth_v;i++){
			tab += "	";
		}
		cout << tab;
		int pass_count_v = node->get_vertical_pass_count_with_id(identifier);
		int stop_count_v = node->get_vertical_stop_count_with_id(identifier);
		int pass_count_h = node->get_horizontal_pass_count_with_id(identifier);
		int stop_count_h = node->get_horizontal_stop_count_with_id(identifier);
		cout << (boost::format("%d [vp:%d,vs:%d,hp:%d,hs:%d][len:%f,self:%f,ch:%f,p:%f]") % node->_identifier % pass_count_v % stop_count_v % pass_count_h % stop_count_h % node->_stick_length % (node->_stick_length - node->_children_stick_length) % node->_children_stick_length % node->_probability).str() << endl;
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			_dump_tssb(child, identifier);
		}
	}
};

#endif