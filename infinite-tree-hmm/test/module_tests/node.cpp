#include  <iostream>
#include  <cassert>
#include  <string>
#include "../../src/ithmm/node.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;
using std::vector;

void test_depth_v(){
	Node* node_1 = new Node();
	Node* node_2 = new Node(node_1);
	Node* node_3 = new Node(node_2);
	assert(node_1->_depth_v == 0);
	assert(node_2->_depth_v == 1);
	assert(node_3->_depth_v == 2);
}
void test_depth_h(){
	Node* node = new Node();
	for(int i = 0 ;i < 100;i++){
		node->generate_child();
	}
	assert(node->_children.size() == 100);	
	for(int i = 0 ;i < 100;i++){
		Node* child = node->_children[i];
		assert(child->_depth_h == i);
	}
}

void test_add_customer_to_vertical_crp(){
	Node* node = new Node();
	double concentration = 0.2;
	double g0 = 0.001;
	bool new_table_generated;
	for(int i = 0 ;i < 100;i++){
		node->add_customer_to_vertical_crp(concentration, g0, new_table_generated);
	}
	assert(node->_pass_count_v == 0);
	assert(node->_pass_count_h == 0);
	assert(node->_stop_count_v == 100);
	assert(node->_stop_count_h == 0);
	assert(node->_table_v->_num_customers == 100);
}

void test_add_customer_to_vertical_crp_recursively(){
	Node* parent = new Node();
	Node* child = new Node(parent);
	double concentration = 0.2;
	double g0 = 0.001;
	bool new_table_generated;
	for(int i = 0 ;i < 100;i++){
		child->add_customer_to_vertical_crp(concentration, g0, new_table_generated);
	}
	assert(child->_pass_count_v == 0);
	assert(child->_pass_count_h == 0);
	assert(child->_stop_count_v == 100);
	assert(child->_stop_count_h == 0);
	assert(child->_table_v->_num_customers == 100);
	assert(parent->_pass_count_v == 100);
	assert(parent->_pass_count_h == 0);
	assert(parent->_stop_count_v == 0);
	assert(parent->_stop_count_h == 0);
	assert(parent->_table_v->_num_customers == 0);
}

void test_add_customer_to_horizontal_crp(){
	Node* node = new Node();
	double concentration = 0.2;
	double g0 = 0.001;
	bool new_table_generated;
	for(int i = 0 ;i < 100;i++){
		node->add_customer_to_horizontal_crp(concentration, g0, new_table_generated);
	}
	assert(node->_pass_count_v == 0);
	assert(node->_pass_count_h == 0);
	assert(node->_stop_count_v == 0);
	assert(node->_stop_count_h == 100);
	assert(node->_table_h->_num_customers == 100);
}

void test_add_customer_to_horizontal_crp_recursively(){
	Node* parent = new Node();
	parent->generate_child();
	parent->generate_child();
	Node* child = new Node(parent);
	parent->add_child(child);
	parent->generate_child();
	double concentration = 0.2;
	double g0 = 0.001;
	bool new_table_generated;
	for(int i = 0 ;i < 100;i++){
		child->add_customer_to_horizontal_crp(concentration, g0, new_table_generated);
	}
	assert(child->_pass_count_v == 0);
	assert(child->_pass_count_h == 0);
	assert(child->_stop_count_v == 0);
	assert(child->_stop_count_h == 100);
	assert(parent->_pass_count_v == 0);
	assert(parent->_pass_count_h == 0);
	assert(parent->_stop_count_v == 0);
	assert(parent->_stop_count_h == 100);
	assert(parent->_children[0]->_pass_count_v == 0);
	assert(parent->_children[0]->_pass_count_h == 100);
	assert(parent->_children[0]->_stop_count_v == 0);
	assert(parent->_children[0]->_stop_count_h == 0);
	assert(parent->_children[1]->_pass_count_v == 0);
	assert(parent->_children[1]->_pass_count_h == 100);
	assert(parent->_children[1]->_stop_count_v == 0);
	assert(parent->_children[1]->_stop_count_h == 0);
	assert(parent->_children[3]->_pass_count_v == 0);
	assert(parent->_children[3]->_pass_count_h == 0);
	assert(parent->_children[3]->_stop_count_v == 0);
	assert(parent->_children[3]->_stop_count_h == 0);
	child->generate_child();
	child->generate_child();
	Node* grandson = new Node(child);
	child->add_child(grandson);
	child->generate_child();
	for(int i = 0 ;i < 100;i++){
		grandson->add_customer_to_horizontal_crp(concentration, g0, new_table_generated);
	}
	cout << child->_pass_count_v << endl;
	cout << child->_pass_count_h << endl;
	cout << child->_stop_count_v << endl;
	cout << child->_stop_count_h << endl;
	cout << child->_table_h->_num_customers << endl;
	cout << parent->_pass_count_v << endl;
	cout << parent->_pass_count_h << endl;
	cout << parent->_stop_count_v << endl;
	cout << parent->_stop_count_h << endl;
	cout << parent->_table_h->_num_customers << endl;
	assert(grandson->_pass_count_v == 0);
	assert(grandson->_pass_count_h == 0);
	assert(grandson->_stop_count_v == 0);
	assert(grandson->_stop_count_h == 100);
	assert(child->_pass_count_v == 0);
	assert(child->_pass_count_h == 0);
	assert(child->_stop_count_v == 0);
	assert(child->_stop_count_h == 200);
	assert(parent->_pass_count_v == 0);
	assert(parent->_pass_count_h == 0);
	assert(parent->_stop_count_v == 0);
	assert(parent->_stop_count_h == 200);
	assert(parent->_children[0]->_pass_count_v == 0);
	assert(parent->_children[0]->_pass_count_h == 200);
	assert(parent->_children[0]->_stop_count_v == 0);
	assert(parent->_children[0]->_stop_count_h == 0);
	assert(parent->_children[1]->_pass_count_v == 0);
	assert(parent->_children[1]->_pass_count_h == 200);
	assert(parent->_children[1]->_stop_count_v == 0);
	assert(parent->_children[1]->_stop_count_h == 0);
	assert(parent->_children[3]->_pass_count_v == 0);
	assert(parent->_children[3]->_pass_count_h == 0);
	assert(parent->_children[3]->_stop_count_v == 0);
	assert(parent->_children[3]->_stop_count_h == 0);
	assert(child->_children[0]->_pass_count_v == 0);
	assert(child->_children[0]->_pass_count_h == 100);
	assert(child->_children[0]->_stop_count_v == 0);
	assert(child->_children[0]->_stop_count_h == 0);
	assert(child->_children[1]->_pass_count_v == 0);
	assert(child->_children[1]->_pass_count_h == 100);
	assert(child->_children[1]->_stop_count_v == 0);
	assert(child->_children[1]->_stop_count_h == 0);
	assert(child->_children[3]->_pass_count_v == 0);
	assert(child->_children[3]->_pass_count_h == 0);
	assert(child->_children[3]->_stop_count_v == 0);
	assert(child->_children[3]->_stop_count_h == 0);
}

int main(){
	test_depth_v();
	cout << "OK" << endl;
	test_depth_h();
	cout << "OK" << endl;
	test_add_customer_to_vertical_crp();
	cout << "OK" << endl;
	test_add_customer_to_vertical_crp_recursively();
	cout << "OK" << endl;
	test_add_customer_to_horizontal_crp();
	cout << "OK" << endl;
	test_add_customer_to_horizontal_crp_recursively();
	cout << "OK" << endl;
	return 0;
}