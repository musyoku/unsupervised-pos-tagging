#include  <iostream>
#include  <chrono>
#include <algorithm>
#include "core/ithmm.h"
#include "core/cprintf.h"
using namespace std;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->sample_node_on_structure_tssb();
		model->add_customer_to(node->_transition_tssb_myself);
	}
}

void test1(iTHMM* model){
	add_customer(model, 100);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
	}
}

void test2(iTHMM* model){
	add_customer(model, 20);
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(11);
	assert(target_on_structure != NULL);
	for(int n = 0;n < 100;n++){
		model->add_customer_to(target_on_structure);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = target_on_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}

	for(int n = 0;n < 100;n++){
		model->remove_htssb_customer_from_node(target_on_structure);
	}

	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	parent = target_on_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test4(iTHMM* model){
	add_customer(model, 20);
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(14);
	assert(target_on_structure != NULL);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = target_on_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
	model->remove_clustering_customer_from_node(target_on_structure);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	for(const auto child: model->_structure_tssb->_root->_children){
		c_printf("[*]%d\n", child->_identifier);
		child->_transition_tssb->dump();
	}
}

void test5(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 10000;n++){
		Node* node = model->sample_node_on_structure_tssb();
		model->add_customer_to(node);
		nodes.push_back(node);
	}
	for(auto node: nodes){
		model->remove_clustering_customer_from_node(node);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
}

void test6(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node_on_structure_tssb();
		model->add_customer_to(node);
		nodes.push_back(node);
	}
	int rand_index = Sampler::uniform_int(0, nodes.size());
	Node* back = nodes[rand_index];
	nodes.erase(nodes.begin() + rand_index);
	nodes.pop_back();
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to(node);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to(node);
		}
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

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
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

	Node* parent = back;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test7(iTHMM* model){
	add_customer(model, 20);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(14);
	assert(target_on_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_on_structure->dump();
	target_on_structure->_transition_tssb_myself->dump();
	for(int i = 0;i < 100;i++){
		model->add_customer_to(target_on_structure->_transition_tssb_myself);
	}
	c_printf("[*]%s\n", "transition");
	target_on_structure->_transition_tssb->dump();
	Node* parent_on_structure = target_on_structure->_parent;
	while(parent_on_structure){
		c_printf("[*]%d\n", parent_on_structure->_identifier);
		parent_on_structure->_transition_tssb->dump();
		parent_on_structure = parent_on_structure->_parent;
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_htssb_vertical_sbr_ratio_on_node(target_on_structure->_transition_tssb_myself);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;
}

void test8(iTHMM* model){
	add_customer(model, 20);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(14);
	assert(target_on_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_on_structure->dump();
	for(int i = 0;i < 100;i++){
		model->add_customer_to(target_on_structure);
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_htssb_horizontal_sbr_ratio_on_node(target_on_structure);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;
}

void test9(iTHMM* model){
	add_customer(model, 3000);
	double uniform = 0;
	model->_structure_tssb->update_stick_length();
	model->_structure_tssb->dump();
	vector<Node*> nodes_true;
	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes_true);
	int num_nodes_true = model->_structure_tssb->get_num_nodes();
	int num_nodes = 0;
	int prev_id = -1;
	for(int i = 0;i < 100000000;i++){
		uniform = i / 100000000.0;
		Node* node = model->retrospective_sampling_on_tssb(uniform, model->_structure_tssb);
		if(node->_identifier != prev_id){
			cout << uniform << ": " << node->_identifier << endl;
			prev_id = node->_identifier;
			num_nodes += 1;
			nodes.push_back(node);
		}
	}
	int limit = (nodes_true.size() > nodes.size()) ? nodes_true.size() : nodes.size();
	for(int i = 0;i < limit;i++){
		int identifier = (i < nodes.size()) ? nodes[i]->_identifier : -1;
		int identifier_true = (i < nodes_true.size()) ? nodes_true[i]->_identifier : -1;
		cout << identifier << " : " << identifier_true << endl;
	}
	cout << num_nodes << " == " << num_nodes_true << endl;
	assert(num_nodes == num_nodes_true);
}

void test10(iTHMM* model){
	add_customer(model, 10);
	double uniform = 0.91;
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(12);
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
	target->_transition_tssb->update_stick_length();
	Node* node = model->retrospective_sampling_on_tssb(uniform, target->_transition_tssb);
	c_printf("[*]%s\n", "sampled");
	node->dump();
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test11(iTHMM* model){
	for(int n = 0;n < 10;n++){
		Node* node = model->sample_node_on_structure_tssb();
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = model->_structure_tssb->find_node_with_id(12);
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		Node* myself = parent->_transition_tssb_myself;
		myself->dump();
		if(myself->_identifier != 1){
			assert(myself->_parent_transition_tssb_myself != NULL);
			myself->_parent_transition_tssb_myself->dump();
		}
		assert(myself->_structure_tssb_myself != NULL);
		myself->_structure_tssb_myself->dump();
		parent = parent->_parent;
	}
}

void test12(iTHMM* model){
	for(int n = 0;n < 10;n++){
		Node* node = model->sample_node_on_structure_tssb();
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(12);
	model->delete_node_if_needed(target);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = model->_structure_tssb->find_node_with_id(14);
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		Node* myself = parent->_transition_tssb_myself;
		myself->dump();
		if(myself->_identifier != 1){
			assert(myself->_parent_transition_tssb_myself != NULL);
			myself->_parent_transition_tssb_myself->dump();
		}
		assert(myself->_structure_tssb_myself != NULL);
		myself->_structure_tssb_myself->dump();
		parent = parent->_parent;
	}
}

int main(){
	iTHMM* model = new iTHMM();
	test11(model);
	return 0;
}