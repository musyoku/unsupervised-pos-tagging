#include  <iostream>
#include  <chrono>
#include "core/ithmm.h"
#include "core/cprintf.h"
using namespace std;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->sample_node();
		model->_clustering_tssb->add_customer_to_node(node);
	}
}

void test1(iTHMM* model){
	add_customer(model, 100);
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();

	vector<Node*> nodes;
	model->_clustering_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
	}
}

int main(){
	iTHMM* model = new iTHMM();
	test1(model);
	return 0;
}