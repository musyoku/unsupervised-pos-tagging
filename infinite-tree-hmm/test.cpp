#include  <iostream>
#include  <chrono>
#include "core/ithmm.h"
#include "core/cprintf.h"
using namespace std;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->sample_node();
		model->add_clustering_customer_to_node(node);
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

void test2(iTHMM* model){
	add_customer(model, 20);
	Node* target_on_cluster = model->_clustering_tssb->find_node_with_id(11);
	assert(target_on_cluster != NULL);
	for(int n = 0;n < 100;n++){
		model->add_htssb_customer_to_node(target_on_cluster);
	}
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
	Node* parent = target_on_cluster;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}

	for(int n = 0;n < 100;n++){
		model->remove_htssb_customer_from_node(target_on_cluster);
	}

	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
	parent = target_on_cluster;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}

	double ratio = target_on_cluster->_transition_tssb->compute_expectation_of_htssb_vertical_sbr_ratio(target_on_cluster->_transition_tssb_myself);
	cout << ratio << endl;
}

int main(){
	iTHMM* model = new iTHMM();
	test2(model);
	return 0;
}