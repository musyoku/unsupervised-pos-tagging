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
	root_in_structure->_transition_tssb->_owner_id = root_in_structure->_identifier;
	root_in_structure->_myself_in_transition_tssb = root_in_htssb;

}

void test_generate_transition_tssb_belonging_to(){
	iTHMM* ithmm = new iTHMM();
	TSSB* tssb = ithmm->generate_transition_tssb_belonging_to(ithmm->_root_in_structure);
	assert(tssb->_owner == ithmm->_root_in_structure);
	assert(tssb->_root->_identifier == ithmm->_root_in_structure->_identifier);
	assert(tssb->_root->get_htssb_owner_node_id() == ithmm->_root_in_structure->_identifier);
}

void test_generate_and_add_new_child_to(){
	iTHMM* ithmm = new iTHMM();
	Node* node = ithmm->generate_and_add_new_child_to(ithmm->_root_in_htssb);
	assert(node->_depth_v == 1);
}

int main(){
	test_copy_children_in_structure_to_transition_tssb();
	cout << "OK" << endl;
	return 0;
}