#include  <iostream>
#include  <chrono>
#include "core/ithmm.h"
using namespace std;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->_clustering_tssb->sample_node();
		model->_clustering_tssb->add_customer_to_node(node);
	}
}

void test1(iTHMM* model){
	add_customer(model, 100);
	model->_clustering_tssb->dump();
}

int main(){
	iTHMM* model = new iTHMM();
	test1(model);
	return 0;
}