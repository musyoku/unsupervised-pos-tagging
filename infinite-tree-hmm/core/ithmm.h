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
	}
	void udpate_htssb_structure(Node* new_node){
		// クラスタリング用TSSBの全ノードを収集
		vector<Node*> nodes;
		_clustering_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(auto node: nodes){
			if(node->_transition_tssb == NULL){

			}
		}
	}
};

#endif