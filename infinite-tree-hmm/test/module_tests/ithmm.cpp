#include  <iostream>
#include  <cassert>
#include  <string>
#include "../../src/ithmm/ithmm.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;
using std::vector;

void test_copy_children_in_structure_to_transition_tssb(){
	Node* root_in_structure = new Node(NULL);

	Node* root_in_htssb = new Node(NULL, root_in_structure->_identifier);
	root_in_htssb->_htssb_owner_node_in_structure = root_in_structure;
	root_in_htssb->_myself_in_structure_tssb = root_in_structure;

	root_in_structure->_transition_tssb = new TSSB(root_in_htssb);
	root_in_structure->_transition_tssb->_owner = root_in_structure;
	root_in_structure->_transition_tssb->_is_htssb = true;
	root_in_structure->_myself_in_transition_tssb = root_in_htssb;

	for(int i = 0;i < 100;i++){
		root_in_structure->generate_child();
	}

	Node* child_in_structure = root_in_structure->_children[50];
	for(int i = 0;i < 100;i++){
		child_in_structure->generate_child();
	}

	copy_children_in_structure_to_transition_tssb(root_in_structure, root_in_htssb);
	assert(root_in_htssb->_children.size() == 100);

	for(int i = 0;i < 100;i++){
		assert(root_in_structure->_children[i]->_identifier == root_in_htssb->_children[i]->_identifier);
	}
	assert(root_in_htssb->_children[50]->_children.size() == 100);
	for(int i = 0;i < 100;i++){
		assert(child_in_structure->_children[i]->_identifier == root_in_htssb->_children[50]->_children[i]->_identifier);
	}
}

void test_generate_transition_tssb_belonging_to(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	for(int i = 0;i < 100;i++){
		root_in_structure->generate_child();
	}
	TSSB* tssb = ithmm->generate_transition_tssb_belonging_to(root_in_structure);
	Node* root_in_htssb = tssb->_root;
	assert(root_in_structure->_children.size() == 100);
	assert(root_in_htssb->_children.size() == 100);
	assert(root_in_htssb != NULL);
	for(int i = 0;i < 100;i++){
		assert(root_in_structure->_children[i]->_identifier == root_in_htssb->_children[i]->_identifier);
	}

	ithmm->generate_and_add_new_child_to(root_in_structure);
}

void test_generate_and_add_new_child_to(){
	iTHMM* ithmm = new iTHMM();
	Node* node = ithmm->generate_and_add_new_child_to(ithmm->_root_in_htssb);
	assert(node->_depth_v == 1);
}

int main(){
	test_copy_children_in_structure_to_transition_tssb();
	cout << "OK" << endl;
	test_generate_transition_tssb_belonging_to();
	cout << "OK" << endl;
	return 0;
}