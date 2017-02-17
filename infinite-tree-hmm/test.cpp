#include  <iostream>
#include "core/htssb.h"
using namespace std;

void add_customer(HTSSB* model){
	for(int n = 0;n < 10;n++){
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
	double uniform = Sampler::uniform(0, 1);
	uniform = 0.75 + 0.041667 + 0.009259 + 0.001;
	Node* node = model->retrospective_sampling(uniform);
	cout << uniform << endl;
	cout << node->_identifier << endl;
	model->update_stick_length();
	model->dump_tssb(CLUSTERING_TSSB_ID);
}
int main(){
	HTSSB* model = new HTSSB();
	add_customer(model);
	test2(model);
	return 0;
}