#include  <iostream>
#include "core/htssb.h"
using namespace std;

int main(){
	HTSSB* model = new HTSSB();
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node();
		cout << node << endl;
		model->add_customer_to_node(node);
		cout << node->_identifier << endl;
	}
	cout << "depth: " << model->get_max_depth() << endl;
	model->dump_nodes();
	return 0;
}