#include  <iostream>
#include  <chrono>
#include "core/htssb.h"
using namespace std;

void add_customer(HTSSB* model){
	for(int n = 0;n < 5000;n++){
		Node* node = model->sample_node();
		model->add_customer_to_clustering_node(node);
	}
}

void test1(HTSSB* model){
	Node* target = NULL;
	vector<Node*> nodes;
	for(int n = 0;n < 20;n++){
		Node* node = model->sample_node();
		model->add_customer_to_clustering_node(node);
		nodes.push_back(node);
		if(node->_identifier == 14){
			target = node;
		}
	}
	assert(target != NULL);
	model->add_customer_to_htssb_node(target);
	model->add_customer_to_htssb_node(target);
	model->add_customer_to_htssb_node(target);
	model->add_customer_to_htssb_node(target);
	model->add_customer_to_htssb_node(target);
	for(const auto &elem: target->_stop_count_v){
		int identifier = elem.first;
		c_printf("[*]%d\n", identifier);
		model->dump_tssb(identifier);
	}

	cout << "depth: " << model->get_max_depth() << endl;
	model->dump_tssb(CLUSTERING_TSSB_ID);
	// for(int n = 0;n < nodes.size();n++){
	// 	Node* node = nodes[n];
	// 	model->remove_customer_from_clustering_node(node);
	// 	// cout << "removing " << node->_identifier << " ..." << endl;
	// 	// model->dump_tssb(CLUSTERING_TSSB_ID);
	// }
	// cout << "depth: " << model->get_max_depth() << endl;
	// model->dump_tssb(CLUSTERING_TSSB_ID);
}
void test2(HTSSB* model){
	add_customer(model);
	double uniform = 0;
	model->update_stick_length();
	model->dump_tssb(CLUSTERING_TSSB_ID);
	vector<Node*> nodes_true;
	vector<Node*> nodes;
	model->enumerate_nodes_from_left_to_right(nodes_true);
	int num_nodes_true = model->get_num_nodes();
	int num_nodes = 0;
	int prev_id = -1;
	for(int i = 0;i < 10000000;i++){
		uniform = i / 10000000.0;
		Node* node = model->retrospective_sampling(uniform);
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
void test3(HTSSB* model){
	for(int n = 0;n < 10;n++){
		Node* node = model->sample_node();
		model->add_customer_to_clustering_node(node);
	}
	Node* target = model->find_node_with_id(11);
	assert(target != NULL);
	for(int n = 0;n < 100;n++){
		model->add_customer_to_htssb_node(target);
	}
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		model->dump_tssb(parent->_identifier);
		parent = parent->_parent;
	}

	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 10000;i++){
		ratio = model->compute_expectation_of_htssb_vertical_sbr_ratio(target);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;

}
int main(){
	HTSSB* model = new HTSSB();
	test3(model);
	return 0;
}