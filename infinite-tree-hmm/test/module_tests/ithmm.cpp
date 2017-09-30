#include <iostream>
#include <cassert>
#include <string>
#include <memory>
#include "../../src/ithmm/ithmm.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;
using std::vector;
using std::unique_ptr;

void test_copy_children_in_structure_to_transition_tssb(){
	Node* root_in_structure = new Node(NULL);
	root_in_structure->set_as_structure_node();

	Node* root_in_htssb = new Node(NULL, root_in_structure->_identifier);
	root_in_htssb->set_as_htssb_node();
	root_in_htssb->set_htssb_owner_node_in_structure(root_in_structure);
	root_in_htssb->set_myself_in_structure_tssb(root_in_structure);

	TSSB* transition_tssb = new TSSB(root_in_htssb);;
	root_in_structure->set_transition_tssb(transition_tssb);
	transition_tssb->_owner = root_in_structure;
	transition_tssb->_is_htssb = true;
	root_in_structure->set_myself_in_transition_tssb(root_in_htssb);

	for(int i = 0;i < 100;i++){
		Node* child = root_in_structure->generate_child();
		child->set_as_structure_node();
		root_in_structure->add_child(child);
	}

	Node* child_in_structure = root_in_structure->_children[50];
	for(int i = 0;i < 100;i++){
		Node* grandson = child_in_structure->generate_child();
		grandson->set_as_structure_node();
		child_in_structure->add_child(grandson);
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
	Node::_delete_all_children(root_in_htssb);
	delete root_in_htssb;
	Node::_delete_all_children(root_in_structure);
	delete root_in_structure;
}
void test_generate_and_add_child_to_parent_in_structure(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	{
		Node* child_in_structure = ithmm->_generate_and_add_child_to_parent_in_structure(root_in_structure);
		assert(child_in_structure->is_structure_node() == true);
	}
	{
		Node* child_in_structure = ithmm->_generate_and_add_child_to_parent_in_structure(root_in_htssb);
		assert(child_in_structure->is_structure_node() == true);
	}
	{
		Node* child_in_structure = ithmm->_generate_and_add_child_to_parent_in_structure(root_in_bos);
		assert(child_in_structure->is_structure_node() == true);
	}
	delete ithmm;
}

void test_generate_transition_tssb_belonging_to(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;

	Node* child_1_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	assert(root_in_structure->_children.size() == 1);
	assert(root_in_htssb->_children.size() == 1);
	assert(root_in_bos->_children.size() == 1);
	assert(child_1_in_structure->is_structure_node());
	assert(child_1_in_structure->get_transition_tssb() != NULL);
	assert(child_1_in_structure->_parent == root_in_structure);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_in_htssb->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == root_in_bos->_children[i]->_identifier);
	}
	TSSB* root_tssb = root_in_structure->get_transition_tssb();
	assert(root_tssb != NULL);
	TSSB* child_1_tssb = child_1_in_structure->get_transition_tssb();
	assert(child_1_tssb != NULL);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_1_tssb->_root->_children[i]->_identifier);
	}

	Node* child_2_in_htssb = ithmm->generate_and_add_new_child_to(root_in_htssb);
	assert(root_in_structure->_children.size() == 2);
	assert(root_in_htssb->_children.size() == 2);
	assert(root_in_bos->_children.size() == 2);
	assert(child_2_in_htssb->is_htssb_node());
	assert(child_2_in_htssb->_parent == root_in_htssb);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_in_htssb->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == root_in_bos->_children[i]->_identifier);
	}
	Node* child_2_in_structure = child_2_in_htssb->get_myself_in_structure_tssb();
	assert(child_2_in_structure != NULL);
	TSSB* child_2_tssb = child_2_in_structure->get_transition_tssb();
	assert(child_2_tssb != NULL);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_1_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_2_tssb->_root->_children[i]->_identifier);
	}

	Node* child_3_in_bos = ithmm->generate_and_add_new_child_to(root_in_bos);
	assert(root_in_structure->_children.size() == 3);
	assert(root_in_htssb->_children.size() == 3);
	assert(root_in_bos->_children.size() == 3);
	assert(child_3_in_bos->is_bos_tssb_node());
	assert(child_3_in_bos->_parent == root_in_bos);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_in_htssb->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == root_in_bos->_children[i]->_identifier);
	}
	Node* child_3_in_structure = child_3_in_bos->get_myself_in_structure_tssb();
	assert(child_3_in_structure != NULL);
	TSSB* child_3_tssb = child_3_in_structure->get_transition_tssb();
	assert(child_3_tssb != NULL);
	for(int i = 0;i < root_in_structure->_children.size();i++){
		assert(root_in_structure->_children[i]->_identifier == root_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_1_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_2_tssb->_root->_children[i]->_identifier);
		assert(root_in_structure->_children[i]->_identifier == child_3_tssb->_root->_children[i]->_identifier);
	}

	Node* grandson_1_in_structure = ithmm->generate_and_add_new_child_to(child_1_in_structure);
	assert(root_in_structure->_children.size() == 3);
	assert(root_in_htssb->_children.size() == 3);
	assert(root_in_bos->_children.size() == 3);
	assert(child_1_in_structure->_children.size() == 1);
	assert(grandson_1_in_structure->is_structure_node());
	assert(grandson_1_in_structure->get_transition_tssb() != NULL);

	TSSB* grandson_1_tssb = grandson_1_in_structure->get_transition_tssb();
	assert(grandson_1_tssb != NULL);
	assert(grandson_1_tssb->_root->_children[0]->_children.size() == 1);
	assert(child_1_tssb->_root->_children[0]->_children.size() == 1);
	assert(child_2_tssb->_root->_children[0]->_children.size() == 1);
	assert(child_3_tssb->_root->_children[0]->_children.size() == 1);
	for(int i = 0;i < child_1_in_structure->_children.size();i++){
		assert(child_1_in_structure->_children[i]->_identifier == grandson_1_tssb->_root->_children[0]->_children[i]->_identifier);
		assert(child_1_in_structure->_children[i]->_identifier == child_1_tssb->_root->_children[0]->_children[i]->_identifier);
		assert(child_1_in_structure->_children[i]->_identifier == child_2_tssb->_root->_children[0]->_children[i]->_identifier);
		assert(child_1_in_structure->_children[i]->_identifier == child_3_tssb->_root->_children[0]->_children[i]->_identifier);
	}

	assert(root_in_structure->_hpylm != NULL);
	assert(child_1_in_structure->_hpylm != NULL);
	assert(child_2_in_structure->_hpylm != NULL);
	assert(child_3_in_structure->_hpylm != NULL);
	assert(grandson_1_in_structure->_hpylm != NULL);

	delete ithmm;
}

void test_pointers(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;

	Node* child_1_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_2_in_htssb = ithmm->generate_and_add_new_child_to(root_in_htssb);
	Node* child_2_in_structure = child_2_in_htssb->get_myself_in_structure_tssb();

	Node* child_3_in_bos = ithmm->generate_and_add_new_child_to(root_in_bos);
	Node* child_3_in_structure = ithmm->_structure_tssb->find_node_by_tracing_horizontal_indices(child_3_in_bos);
	assert(child_3_in_structure != NULL);
	Node* child_3_in_htssb = child_3_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(child_3_in_structure);
	assert(child_3_in_htssb != NULL);

	assert(child_3_in_structure->_parent == root_in_structure);
	assert(child_3_in_structure->get_myself_in_transition_tssb() == child_3_in_htssb);
	assert(child_3_in_structure->get_myself_in_bos_tssb() == child_3_in_bos);
	assert(child_3_in_htssb->get_myself_in_structure_tssb() == child_3_in_structure);
	assert(child_3_in_htssb->get_myself_in_parent_transition_tssb() == child_3_in_structure->_parent->get_transition_tssb()->find_node_by_tracing_horizontal_indices(child_3_in_structure));
	assert(child_3_in_structure->get_transition_tssb() != NULL);

	Node* grandson_1_in_structure = ithmm->generate_and_add_new_child_to(child_1_in_structure);
	Node* grandson_1_in_htssb = grandson_1_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_1_in_structure);
	Node* grandson_1_in_bos = ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(grandson_1_in_structure);
	assert(grandson_1_in_htssb != NULL);
	assert(grandson_1_in_bos != NULL);

	assert(grandson_1_in_structure->_parent == child_1_in_structure);
	assert(grandson_1_in_structure->get_myself_in_transition_tssb() == grandson_1_in_htssb);
	assert(grandson_1_in_structure->get_myself_in_bos_tssb() == grandson_1_in_bos);
	assert(grandson_1_in_htssb->get_myself_in_structure_tssb() == grandson_1_in_structure);
	assert(grandson_1_in_htssb->get_myself_in_parent_transition_tssb() == grandson_1_in_structure->_parent->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_1_in_structure));
	assert(grandson_1_in_structure->get_transition_tssb() != NULL);
	delete ithmm;
}

void test_add_customer_to_hpylm(){
	iTHMM* ithmm = new iTHMM();
	ithmm->set_word_g0(0.001);
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->add_customer_to_hpylm(grandson_in_structure, 1);
	assert(grandson_in_structure->_hpylm->get_num_customers() == 1);
	assert(child_in_structure->_hpylm->get_num_customers() == 1);
	assert(root_in_structure->_hpylm->get_num_customers() == 1);
	delete ithmm;
}

void test_remove_customer_from_hpylm(){
	iTHMM* ithmm = new iTHMM();
	ithmm->set_word_g0(0.001);
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->add_customer_to_hpylm(grandson_in_structure, 1);
	ithmm->remove_customer_from_hpylm(grandson_in_structure, 1);
	assert(grandson_in_structure->_hpylm->get_num_customers() == 0);
	assert(child_in_structure->_hpylm->get_num_customers() == 0);
	assert(root_in_structure->_hpylm->get_num_customers() == 0);
	assert(grandson_in_structure->_hpylm->get_num_tables() == 0);
	assert(child_in_structure->_hpylm->get_num_tables() == 0);
	assert(root_in_structure->_hpylm->get_num_tables() == 0);
	delete ithmm;
}

void test_hpylm_concentration(){
	iTHMM* ithmm = new iTHMM();
	ithmm->set_word_g0(0.001);
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->_hpylm_d_m[0] = 0.1;
	ithmm->_hpylm_d_m[1] = 0.1;
	ithmm->_hpylm_d_m[2] = 0.1;

	ithmm->_hpylm_theta_m[0] = 10000;
	ithmm->_hpylm_theta_m[1] = 10000;
	ithmm->_hpylm_theta_m[2] = 10000;
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_hpylm(grandson_in_structure, 1);
	}
	assert(grandson_in_structure->_hpylm->get_num_customers() == 100);
	assert(child_in_structure->_hpylm->get_num_customers() == grandson_in_structure->_hpylm->get_num_tables());
	assert(root_in_structure->_hpylm->get_num_customers() == child_in_structure->_hpylm->get_num_tables());
	int num_tables_in_child_huge = child_in_structure->_hpylm->get_num_tables();
	int num_tables_in_grandson_huge = grandson_in_structure->_hpylm->get_num_tables();
	for(int i = 0;i < 100;i++){
		ithmm->remove_customer_from_hpylm(grandson_in_structure, 1);
	}
	assert(grandson_in_structure->_hpylm->get_num_customers() == 0);
	assert(child_in_structure->_hpylm->get_num_customers() == 0);
	assert(root_in_structure->_hpylm->get_num_customers() == 0);
	assert(grandson_in_structure->_hpylm->get_num_tables() == 0);
	assert(child_in_structure->_hpylm->get_num_tables() == 0);
	assert(root_in_structure->_hpylm->get_num_tables() == 0);

	ithmm->_hpylm_theta_m[0] = 1;
	ithmm->_hpylm_theta_m[1] = 1;
	ithmm->_hpylm_theta_m[2] = 1;
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_hpylm(grandson_in_structure, 1);
	}
	assert(grandson_in_structure->_hpylm->get_num_customers() == 100);
	assert(child_in_structure->_hpylm->get_num_customers() == grandson_in_structure->_hpylm->get_num_tables());
	assert(root_in_structure->_hpylm->get_num_customers() == child_in_structure->_hpylm->get_num_tables());
	int num_tables_in_child_small = child_in_structure->_hpylm->get_num_tables();
	int num_tables_in_grandson_small = grandson_in_structure->_hpylm->get_num_tables();
	for(int i = 0;i < 100;i++){
		ithmm->remove_customer_from_hpylm(grandson_in_structure, 1);
	}
	assert(grandson_in_structure->_hpylm->get_num_customers() == 0);
	assert(child_in_structure->_hpylm->get_num_customers() == 0);
	assert(root_in_structure->_hpylm->get_num_customers() == 0);
	assert(grandson_in_structure->_hpylm->get_num_tables() == 0);
	assert(child_in_structure->_hpylm->get_num_tables() == 0);
	assert(root_in_structure->_hpylm->get_num_tables() == 0);
	assert(num_tables_in_child_huge > num_tables_in_child_small);
	assert(num_tables_in_grandson_huge > num_tables_in_grandson_small);
	delete ithmm;
}

int main(){
	test_copy_children_in_structure_to_transition_tssb();
	cout << "OK" << endl;
	test_generate_and_add_child_to_parent_in_structure();
	cout << "OK" << endl;
	test_generate_transition_tssb_belonging_to();
	cout << "OK" << endl;
	test_pointers();
	cout << "OK" << endl;
	test_add_customer_to_hpylm();
	cout << "OK" << endl;
	test_remove_customer_from_hpylm();
	cout << "OK" << endl;
	test_hpylm_concentration();
	cout << "OK" << endl;
	return 0;
}