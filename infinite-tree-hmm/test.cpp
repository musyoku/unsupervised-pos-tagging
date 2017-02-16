#include  <iostream>
#include "core/htssb.h"
using namespace std;

int main(){
	HTSSB* model = new HTSSB();
	Node* target = NULL;
	vector<Node*> nodes;
	for(int n = 0;n < 5;n++){
		Node* node = model->sample_node();
		model->add_customer_to_clustering_node(node);
		nodes.push_back(node);
		if(node->_identifier == 6){
			target = node;
		}
	}
	assert(target != NULL);
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
	return 0;
}