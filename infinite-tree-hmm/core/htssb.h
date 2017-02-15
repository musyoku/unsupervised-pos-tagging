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
	void add_customer_to_node(Node* node){
		double alpha = _alpha * pow(_alpha, node->_depth_v);
		node->add_customer_to_clustering_vertical_crp(alpha);
		node->add_customer_to_clustering_horizontal_crp(_gamma);
	}
	void remove_customer_from_node(Node* node){
		node->remove_customer_from_clustering_vertical_crp();
		node->remove_customer_from_clustering_horizontal_crp();
	}
	Node* sample_node(){
		return _stop_node(_root);
	}
	// 止まるノードを決定する
	Node* _stop_node(Node* node){
		assert(node != NULL);
		double alpha = _alpha * pow(_lambda, node->_depth_v);
		double head = node->compute_expectation_of_clustering_vertical_sbr_ratio(alpha);
		if(head == 0){
			head = 16 * 0 - head;
		}
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
		cout << (boost::format("%d [vp:%d,vs:%d,hp:%d,hs:%d][len:%f]") % node->_identifier % pass_count_v % stop_count_v % pass_count_h % stop_count_h % node->_stick_length).str() << endl;
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			_dump_tssb(child, identifier);
		}
	}
};

#endif