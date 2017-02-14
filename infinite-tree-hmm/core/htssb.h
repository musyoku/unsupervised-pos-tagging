#ifndef _htssb_
#define _htssb_
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
	void add_customer_to_node(Node* node){
		double alpha = _alpha * pow(_alpha, node->_depth_v);
		node->add_customer_to_clustering_vertical_crp(alpha);
		// node->add_customer_to_horizontal_crp(_gamma);
	}
	Node* sample_node(){
		return _stop_node(_root);
	}
	// 止まるノードを決定する
	Node* _stop_node(Node* node){
		assert(node != NULL);
		double alpha = _alpha * pow(_lambda, node->_depth_v);
		double head = node->compute_expectation_of_vertical_sbr_ratio(alpha);
		node->_children_stick_length = 1 - node->_stick_length * head;
		double bernoulli = Sampler::uniform(0, 1);
		if(bernoulli <= head){			// 表が出たらこのノードに降りる
			return node;
		}
		// 子ノードがある場合
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			assert(child != NULL);
			double head = child->compute_expectation_of_horizontal_sbr_ratio(_gamma);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _stop_node(child);
			}
		}
		// ない場合生成しながらコインを投げる
		while(true){
			Node* child = node->generate_child();
			double head = child->compute_expectation_of_horizontal_sbr_ratio(_gamma);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _stop_node(child);
			}
		}
	}
	int get_max_depth(){
		return _get_max_depth(_root);
	}
	int _get_max_depth(Node* node){
		int max_depth = node->_depth_v;
		for(const auto &child: node->_children){
			int depth = _get_max_depth(child);
			if(depth > max_depth){
				max_depth = depth;
			}
		}
		return max_depth;
	}
	void dump_nodes(){
		int max_depth = get_max_depth();
		cout << _root->compute_stop_probability() << endl;
		for(int depth = 1;depth <= max_depth;depth++){
			_dump_node(_root, depth);
			cout << endl;
		}
	}
	void _dump_node(Node* node, int target_depth){
		for(const auto &child: node->_children){
			if(child->_depth_v == target_depth){
				cout << child->compute_stop_probability() << " ";
			}
			_dump_node(child, target_depth);
		}
	}
};

#endif