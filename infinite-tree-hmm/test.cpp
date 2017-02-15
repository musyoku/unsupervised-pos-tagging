#include  <iostream>
#include "core/htssb.h"
using namespace std;

int main(){
	HTSSB* model = new HTSSB();
	vector<Node*> nodes;
	for(int n = 0;n < 1;n++){
		Node* node = model->sample_node();
		model->add_customer_to_node(node);
		nodes.push_back(node);
	}
	cout << "depth: " << model->get_max_depth() << endl;
	model->dump_tssb(CLUSTERING_TSSB_ID);
	for(int n = 0;n < nodes.size();n++){
		Node* node = nodes[n];
		model->remove_customer_from_node(node);
	}
	cout << "depth: " << model->get_max_depth() << endl;
	model->dump_tssb(CLUSTERING_TSSB_ID);
	return 0;
}