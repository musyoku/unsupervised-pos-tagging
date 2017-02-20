#include  <iostream>
#include  <chrono>
#include <algorithm>
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
}

void test4(iTHMM* model){
	add_customer(model, 20);
	Node* target_on_cluster = model->_clustering_tssb->find_node_with_id(14);
	assert(target_on_cluster != NULL);
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
	Node* parent = target_on_cluster;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
	model->remove_clustering_customer_from_node(target_on_cluster);
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
	for(const auto child: model->_clustering_tssb->_root->_children){
		c_printf("[*]%d\n", child->_identifier);
		child->_transition_tssb->dump();
	}
}

void test5(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 10000;n++){
		Node* node = model->sample_node();
		model->add_clustering_customer_to_node(node);
		nodes.push_back(node);
	}
	for(auto node: nodes){
		model->remove_clustering_customer_from_node(node);
	}
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
}

void test6(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node();
		model->add_clustering_customer_to_node(node);
		nodes.push_back(node);
	}
	int rand_index = Sampler::uniform_int(0, nodes.size());
	Node* back = nodes[rand_index];
	nodes.erase(nodes.begin() + rand_index);
	nodes.pop_back();
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_htssb_customer_to_node(node);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_htssb_customer_to_node(node);
		}
	}
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();

	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_htssb_customer_from_node(node);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_htssb_customer_from_node(node);
		}
	}

	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		model->remove_clustering_customer_from_node(node);
	}
	c_printf("[*]%s\n", "back");
	back->dump();
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();

	Node* parent = back;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test7(iTHMM* model){
	add_customer(model, 30);
	c_printf("[*]%s\n", "cluster");
	model->_clustering_tssb->dump();
	Node* target_on_cluster = model->_clustering_tssb->find_node_with_id(19);
	assert(target_on_cluster != NULL);
	c_printf("[*]%s\n", "target");
	target_on_cluster->dump();
	for(int i = 0;i < 100;i++){
		model->add_htssb_customer_to_node(target_on_cluster);
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 10000;i++){
		ratio = model->compute_expectation_of_htssb_vertical_sbr_ratio_on_node(target_on_cluster);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;
}

int main(){
	iTHMM* model = new iTHMM();
	test7(model);
	return 0;
}