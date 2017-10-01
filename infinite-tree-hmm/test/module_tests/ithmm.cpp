#include <iostream>
#include <cassert>
#include <string>
#include <memory>
#include "../../src/ithmm/ithmm.h"
#include "../../src/ithmm/sampler.h"
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

void test_add_customer_to_tssb_node(){
	iTHMM* ithmm = new iTHMM();
	ithmm->set_word_g0(0.001);
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_tssb_node(grandson_in_structure);
	}
	assert(ithmm->_structure_tssb->get_num_customers() == 200);
	assert(ithmm->_structure_tssb->_num_customers == 200);
	assert(ithmm->_bos_tssb->get_num_customers() == 0);
	Node* grandson_in_bos = ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(grandson_in_structure);
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_tssb_node(grandson_in_bos);
	}
	assert(ithmm->_structure_tssb->get_num_customers() == 200);
	assert(ithmm->_structure_tssb->_num_customers == 200);
	assert(ithmm->_bos_tssb->get_num_customers() == 200);
	assert(ithmm->_bos_tssb->_num_customers == 200);
}

void test_remove_customer_from_tssb_node(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_tssb_node(grandson_in_structure);
	}
	Node* grandson_in_bos = ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(grandson_in_structure);
	for(int i = 0;i < 100;i++){
		ithmm->add_customer_to_tssb_node(grandson_in_bos);
	}
	assert(ithmm->_structure_tssb->get_num_customers() == 200);
	assert(ithmm->_structure_tssb->_num_customers == 200);
	assert(ithmm->_bos_tssb->get_num_customers() == 200);
	assert(ithmm->_bos_tssb->_num_customers == 200);
	for(int i = 0;i < 100;i++){
		ithmm->remove_customer_from_tssb_node(grandson_in_structure);
		ithmm->remove_customer_from_tssb_node(grandson_in_bos);
	}
	assert(ithmm->_structure_tssb->get_num_customers() == 0);
	assert(ithmm->_structure_tssb->_num_customers == 0);
	assert(ithmm->_bos_tssb->get_num_customers() == 0);
	assert(ithmm->_bos_tssb->_num_customers == 0);
}

void test_htssb_concentration_vertical(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	ithmm->_strength_v = 1;
	for(int i = 0;i < 100;i++){
		ithmm->_add_customer_to_htssb_vertical_crp(grandson_in_htssb);
	}
	int num_customers_small_1_small = grandson_in_htssb->_table_v->get_num_customers();
	int num_customers_small_2_small = grandson_in_child_htssb->_table_v->get_num_customers();
	int num_customers_small_3_small = grandson_in_root_htssb->_table_v->get_num_customers();
	assert(grandson_in_htssb->_table_v->get_num_customers() == 100);
	assert(grandson_in_child_htssb->_table_v->get_num_customers() == grandson_in_htssb->_table_v->get_num_tables());
	assert(grandson_in_root_htssb->_table_v->get_num_customers() == grandson_in_child_htssb->_table_v->get_num_tables());
	assert(grandson_in_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_tables() == 0);

	for(int i = 0;i < 100;i++){
		ithmm->_remove_customer_from_htssb_vertical_crp(grandson_in_htssb);
	}
	assert(grandson_in_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_tables() == 0);

	ithmm->_strength_v = 10000;
	for(int i = 0;i < 100;i++){
		ithmm->_add_customer_to_htssb_vertical_crp(grandson_in_htssb);
	}
	int num_customers_small_1_huge = grandson_in_htssb->_table_v->get_num_customers();
	int num_customers_small_2_huge = grandson_in_child_htssb->_table_v->get_num_customers();
	int num_customers_small_3_huge = grandson_in_root_htssb->_table_v->get_num_customers();
	assert(grandson_in_htssb->_table_v->get_num_customers() == 100);
	assert(grandson_in_child_htssb->_table_v->get_num_customers() == grandson_in_htssb->_table_v->get_num_tables());
	assert(grandson_in_root_htssb->_table_v->get_num_customers() == grandson_in_child_htssb->_table_v->get_num_tables());
	assert(num_customers_small_1_huge == num_customers_small_1_small);
	assert(num_customers_small_2_huge > num_customers_small_2_small);
	assert(num_customers_small_3_huge > num_customers_small_3_small);
	assert(grandson_in_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_tables() == 0);
	delete ithmm;
}

void test_htssb_concentration_horizontal(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	ithmm->_strength_h = 1;
	for(int i = 0;i < 100;i++){
		ithmm->_add_customer_to_htssb_horizontal_crp(grandson_in_htssb);
	}
	int num_customers_small_1_small = grandson_in_htssb->_table_h->get_num_customers();
	int num_customers_small_2_small = grandson_in_child_htssb->_table_h->get_num_customers();
	int num_customers_small_3_small = grandson_in_root_htssb->_table_h->get_num_customers();
	assert(grandson_in_htssb->_table_h->get_num_customers() == 100);
	assert(grandson_in_child_htssb->_table_h->get_num_customers() == grandson_in_htssb->_table_h->get_num_tables());
	assert(grandson_in_root_htssb->_table_h->get_num_customers() == grandson_in_child_htssb->_table_h->get_num_tables());
	assert(grandson_in_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_tables() == 0);

	for(int i = 0;i < 100;i++){
		ithmm->_remove_customer_from_htssb_horizontal_crp(grandson_in_htssb);
	}
	assert(grandson_in_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_h->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_h->get_num_tables() == 0);

	ithmm->_strength_h = 10000;
	for(int i = 0;i < 100;i++){
		ithmm->_add_customer_to_htssb_horizontal_crp(grandson_in_htssb);
	}
	int num_customers_small_1_huge = grandson_in_htssb->_table_h->get_num_customers();
	int num_customers_small_2_huge = grandson_in_child_htssb->_table_h->get_num_customers();
	int num_customers_small_3_huge = grandson_in_root_htssb->_table_h->get_num_customers();
	assert(grandson_in_htssb->_table_h->get_num_customers() == 100);
	assert(grandson_in_child_htssb->_table_h->get_num_customers() == grandson_in_htssb->_table_h->get_num_tables());
	assert(grandson_in_root_htssb->_table_h->get_num_customers() == grandson_in_child_htssb->_table_h->get_num_tables());
	assert(num_customers_small_1_huge == num_customers_small_1_small);
	assert(num_customers_small_2_huge > num_customers_small_2_small);
	assert(num_customers_small_3_huge > num_customers_small_3_small);
	assert(grandson_in_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_customers() == 0);
	assert(grandson_in_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_child_htssb->_table_v->get_num_tables() == 0);
	assert(grandson_in_root_htssb->_table_v->get_num_tables() == 0);
	delete ithmm;
}

void test_nodes_from_root_to_myself(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	assert(grandson_in_structure->_nodes_from_root_to_myself[0] == root_in_structure);
	assert(grandson_in_structure->_nodes_from_root_to_myself[1] == child_in_structure);
}

void test_compute_expectation_of_vertical_htssb_sbr_ratio(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	double* stop_ratio_over_parent = new double[3];
	double* stop_probability_over_parent = new double[3];
	double ratio_v, rest_stick_length, stop_probability, sum_parent_stop_probability, parent_stop_probability;
	int stop_count, pass_count;

	for(int i = 0;i < 10;i++){
		ithmm->_lambda_alpha = 1;
		ithmm->_strength_v = sampler::uniform_int(1, 1000) / 10.0;

		double alpha = ithmm->_alpha;
		double strength_v = ithmm->_strength_v;

		grandson_in_htssb->_pass_count_v = 20;
		grandson_in_htssb->_stop_count_v = 40;
		grandson_in_htssb->_parent->_pass_count_v = 100;
		grandson_in_htssb->_parent->_stop_count_v = 50;
		grandson_in_htssb->_parent->_parent->_pass_count_v = 80;
		grandson_in_htssb->_parent->_parent->_stop_count_v = 160;

		grandson_in_child_htssb->_pass_count_v = 80;
		grandson_in_child_htssb->_stop_count_v = 40;
		grandson_in_child_htssb->_parent->_pass_count_v = 60;
		grandson_in_child_htssb->_parent->_stop_count_v = 120;
		grandson_in_child_htssb->_parent->_parent->_pass_count_v = 60;
		grandson_in_child_htssb->_parent->_parent->_stop_count_v = 30;

		grandson_in_root_htssb->_pass_count_v = 10;
		grandson_in_root_htssb->_stop_count_v = 30;
		grandson_in_root_htssb->_parent->_pass_count_v = 120;
		grandson_in_root_htssb->_parent->_stop_count_v = 40;
		grandson_in_root_htssb->_parent->_parent->_pass_count_v = 50;
		grandson_in_root_htssb->_parent->_parent->_stop_count_v = 150;

		// root
		stop_count = grandson_in_root_htssb->_parent->_parent->_stop_count_v;
		pass_count = grandson_in_root_htssb->_parent->_parent->_pass_count_v;
		ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
		stop_ratio_over_parent[0] = ratio_v;

		stop_count = grandson_in_root_htssb->_parent->_stop_count_v;
		pass_count = grandson_in_root_htssb->_parent->_pass_count_v;
		ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
		stop_ratio_over_parent[1] = ratio_v;

		stop_count = grandson_in_root_htssb->_stop_count_v;
		pass_count = grandson_in_root_htssb->_pass_count_v;
		ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
		stop_ratio_over_parent[2] = ratio_v;
		
		rest_stick_length = 1;
		for(int m = 0;m < 3;m++){
			ratio_v = stop_ratio_over_parent[m];
			stop_probability = rest_stick_length * ratio_v;
			rest_stick_length *= 1.0 - ratio_v;
			stop_probability_over_parent[m] = stop_probability;
		}

		// child
		sum_parent_stop_probability = 0;

		stop_count = grandson_in_child_htssb->_parent->_parent->_stop_count_v;
		pass_count = grandson_in_child_htssb->_parent->_parent->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[0];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[0] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		stop_count = grandson_in_child_htssb->_parent->_stop_count_v;
		pass_count = grandson_in_child_htssb->_parent->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[1];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[1] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		stop_count = grandson_in_child_htssb->_stop_count_v;
		pass_count = grandson_in_child_htssb->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[2];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[2] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		rest_stick_length = 1;
		for(int m = 0;m < 3;m++){
			ratio_v = stop_ratio_over_parent[m];
			stop_probability = rest_stick_length * ratio_v;
			rest_stick_length *= 1.0 - ratio_v;
			stop_probability_over_parent[m] = stop_probability;
		}

		// grandson
		sum_parent_stop_probability = 0;

		stop_count = grandson_in_htssb->_parent->_parent->_stop_count_v;
		pass_count = grandson_in_htssb->_parent->_parent->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[0];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[0] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		stop_count = grandson_in_htssb->_parent->_stop_count_v;
		pass_count = grandson_in_htssb->_parent->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[1];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[1] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		stop_count = grandson_in_htssb->_stop_count_v;
		pass_count = grandson_in_htssb->_pass_count_v;
		parent_stop_probability = stop_probability_over_parent[2];
		ratio_v = (strength_v * parent_stop_probability + stop_count) / (strength_v * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_v > 0);
		stop_ratio_over_parent[2] = ratio_v;
		sum_parent_stop_probability += parent_stop_probability;

		rest_stick_length = 1;
		for(int m = 0;m < 3;m++){
			ratio_v = stop_ratio_over_parent[m];
			stop_probability = rest_stick_length * ratio_v;
			rest_stick_length *= 1.0 - ratio_v;
			stop_probability_over_parent[m] = stop_probability;
		}

		double ratio = ithmm->compute_expectation_of_vertical_htssb_sbr_ratio(grandson_in_htssb);
		assert(ratio == ratio_v);
	}
	delete[] stop_ratio_over_parent;
	delete[] stop_probability_over_parent;
	delete ithmm;
}

void test_compute_expectation_of_horizontal_htssb_sbr_ratio(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	double* stop_ratio_over_parent = new double[3];
	double* stop_probability_over_parent = new double[3];
	double ratio_h, rest_stick_length, stop_probability, sum_parent_stop_probability, parent_stop_probability, parent_ratio_h;
	int stop_count, pass_count;
	Node* parent_in_htssb;
	Node* child_in_htssb;

	for(int i = 0;i < 10;i++){
		ithmm->_lambda_gamma = 1;
		ithmm->_strength_h = sampler::uniform_int(1, 1000) / 10.0;

		for(Node* child: grandson_in_htssb->_parent->_children){
			child->_pass_count_h = sampler::uniform_int(10, 1000);
			child->_stop_count_h = sampler::uniform_int(10, 1000);
		}
		for(Node* child: grandson_in_child_htssb->_parent->_children){
			child->_pass_count_h = sampler::uniform_int(10, 1000);
			child->_stop_count_h = sampler::uniform_int(10, 1000);
		}
		for(Node* child: grandson_in_root_htssb->_parent->_children){
			child->_pass_count_h = sampler::uniform_int(10, 1000);
			child->_stop_count_h = sampler::uniform_int(10, 1000);
		}

		double gamma = ithmm->_gamma;
		double strength_h = ithmm->_strength_h;

		// root
		parent_in_htssb = grandson_in_root_htssb->_parent;

		child_in_htssb = parent_in_htssb->_children[0];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		ratio_h = (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
		stop_ratio_over_parent[0] = ratio_h;

		child_in_htssb = parent_in_htssb->_children[1];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		ratio_h = (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
		stop_ratio_over_parent[1] = ratio_h;

		child_in_htssb = parent_in_htssb->_children[2];
		assert(child_in_htssb == grandson_in_root_htssb);
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		ratio_h = (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
		stop_ratio_over_parent[2] = ratio_h;

		rest_stick_length = 1;
		for(int m = 0;m < 3;m++){
			double ratio_h = stop_ratio_over_parent[m];
			double stop_probability = rest_stick_length * ratio_h;
			rest_stick_length *= 1.0 - ratio_h;
			stop_probability_over_parent[m] = stop_probability;
		}

		// child
		sum_parent_stop_probability = 0;
		parent_in_htssb = grandson_in_child_htssb->_parent;

		child_in_htssb = parent_in_htssb->_children[0];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[0];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[0] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;

		child_in_htssb = parent_in_htssb->_children[1];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[1];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[1] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;

		child_in_htssb = parent_in_htssb->_children[2];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[2];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[2] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;

		rest_stick_length = 1;
		for(int m = 0;m < 3;m++){
			double ratio_h = stop_ratio_over_parent[m];
			double stop_probability = rest_stick_length * ratio_h;
			rest_stick_length *= 1.0 - ratio_h;
			stop_probability_over_parent[m] = stop_probability;
		}

		// grandson
		sum_parent_stop_probability = 0;
		parent_in_htssb = grandson_in_htssb->_parent;

		child_in_htssb = parent_in_htssb->_children[0];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[0];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[0] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;

		child_in_htssb = parent_in_htssb->_children[1];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[1];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[1] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;

		child_in_htssb = parent_in_htssb->_children[2];
		pass_count = child_in_htssb->_pass_count_h;
		stop_count = child_in_htssb->_stop_count_h;
		parent_stop_probability = stop_probability_over_parent[2];
		ratio_h = (strength_h * parent_stop_probability + stop_count) / (strength_h * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
		assert(ratio_h > 0);
		stop_ratio_over_parent[2] = ratio_h;
		sum_parent_stop_probability += parent_stop_probability;
		
		double ratio = ithmm->compute_expectation_of_horizontal_htssb_sbr_ratio(grandson_in_htssb);
		assert(ratio == ratio_h);
	}

	delete[] stop_ratio_over_parent;
	delete[] stop_probability_over_parent;
	delete ithmm;
}

void test_compute_concentration_vertical_htssb_sbr_ratio(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	ithmm->_lambda_alpha = 1;

	grandson_in_htssb->_pass_count_v = 100;
	grandson_in_htssb->_stop_count_v = 100;
	grandson_in_htssb->_parent->_pass_count_v = 100;
	grandson_in_htssb->_parent->_stop_count_v = 100;
	grandson_in_htssb->_parent->_parent->_pass_count_v = 100;
	grandson_in_htssb->_parent->_parent->_stop_count_v = 100;

	grandson_in_child_htssb->_pass_count_v = 10;
	grandson_in_child_htssb->_stop_count_v = 100;
	grandson_in_child_htssb->_parent->_pass_count_v = 100;
	grandson_in_child_htssb->_parent->_stop_count_v = 10;
	grandson_in_child_htssb->_parent->_parent->_pass_count_v = 100;
	grandson_in_child_htssb->_parent->_parent->_stop_count_v = 10;

	grandson_in_root_htssb->_pass_count_v = 1;
	grandson_in_root_htssb->_stop_count_v = 100;
	grandson_in_root_htssb->_parent->_pass_count_v = 100;
	grandson_in_root_htssb->_parent->_stop_count_v = 1;
	grandson_in_root_htssb->_parent->_parent->_pass_count_v = 100;
	grandson_in_root_htssb->_parent->_parent->_stop_count_v = 1;

	ithmm->_strength_v = 1;
	double ratio_1 = ithmm->compute_expectation_of_vertical_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_v = 10;
	double ratio_2 = ithmm->compute_expectation_of_vertical_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_v = 100;
	double ratio_3 = ithmm->compute_expectation_of_vertical_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_v = 1000;
	double ratio_4 = ithmm->compute_expectation_of_vertical_htssb_sbr_ratio(grandson_in_htssb);
	assert(ratio_4 > ratio_3 && ratio_3 > ratio_2 && ratio_2 > ratio_1);

	delete ithmm;
}

void test_compute_concentration_horizontal_htssb_sbr_ratio(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);

	ithmm->_lambda_gamma = 1;

	for(Node* child: grandson_in_htssb->_parent->_children){
		child->_pass_count_h = 100;
		child->_stop_count_h = 100;
	}
	for(Node* child: grandson_in_child_htssb->_parent->_children){
		child->_pass_count_h = 10;
		child->_stop_count_h = 100;
	}
	for(Node* child: grandson_in_root_htssb->_parent->_children){
		child->_pass_count_h = 1;
		child->_stop_count_h = 100;
	}

	ithmm->_strength_h = 1;
	double ratio_1 = ithmm->compute_expectation_of_horizontal_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_h = 10;
	double ratio_2 = ithmm->compute_expectation_of_horizontal_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_h = 100;
	double ratio_3 = ithmm->compute_expectation_of_horizontal_htssb_sbr_ratio(grandson_in_htssb);
	ithmm->_strength_h = 1000;
	double ratio_4 = ithmm->compute_expectation_of_horizontal_htssb_sbr_ratio(grandson_in_htssb);
	assert(ratio_4 > ratio_3 && ratio_3 > ratio_2 && ratio_2 > ratio_1);

	delete ithmm;
}

void test_sample_node_in_tssb_by_iterating_node(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);
	ithmm->_lambda_alpha = 0.1;
	ithmm->_lambda_gamma = 1;
	ithmm->_strength_h = 10;
	ithmm->_strength_v = 10;
	bool new_table_generated;

	// vertical
	for(int i = 0;i < 100;i++){
		grandson_in_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 1;i++){
		grandson_in_child_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}

	// horizontal
	for(int i = 0;i < 1;i++){
		grandson_in_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(Node* child: grandson_in_root_htssb->_parent->_children){
		for(int i = 0;i < 100;i++){
			child->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
		}
	}

	double avg_depth_1 = 0;
	ithmm->_strength_h = 1;
	ithmm->_strength_v = 1;
	for(int i = 0;i < 100;i++){
		Node* node = ithmm->sample_node_in_htssb(grandson_in_structure->get_transition_tssb(), false);
		avg_depth_1 += node->_depth_v / 100.0;
	}
	double avg_depth_2 = 0;
	ithmm->_strength_h = 10;
	ithmm->_strength_v = 10;
	for(int i = 0;i < 100;i++){
		Node* node = ithmm->sample_node_in_htssb(grandson_in_structure->get_transition_tssb(), false);
		avg_depth_2 += node->_depth_v / 100.0;
	}
	double avg_depth_3 = 0;
	ithmm->_strength_h = 100;
	ithmm->_strength_v = 100;
	for(int i = 0;i < 100;i++){
		Node* node = ithmm->sample_node_in_htssb(grandson_in_structure->get_transition_tssb(), false);
		avg_depth_3 += node->_depth_v / 100.0;
	}
	double avg_depth_4 = 0;
	ithmm->_strength_h = 1000;
	ithmm->_strength_v = 1000;
	for(int i = 0;i < 100;i++){
		Node* node = ithmm->sample_node_in_htssb(grandson_in_structure->get_transition_tssb(), false);
		avg_depth_4 += node->_depth_v / 100.0;
	}
	assert(avg_depth_1 > avg_depth_2 && avg_depth_2 > avg_depth_3 && avg_depth_3 > avg_depth_4);
}

void test_update_stick_length_of_tssb(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);
	ithmm->_lambda_alpha = 0.1;
	ithmm->_lambda_gamma = 1;
	ithmm->_strength_h = 1;
	ithmm->_strength_v = 1;
	bool new_table_generated;

	// vertical
	for(int i = 0;i < 100;i++){
		grandson_in_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 1;i++){
		grandson_in_child_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}

	// horizontal
	for(int i = 0;i < 1;i++){
		grandson_in_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(Node* child: grandson_in_root_htssb->_parent->_children){
		for(int i = 0;i < 100;i++){
			child->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
		}
	}

	double total_stick_length = 0.5;
	TSSB* htssb = grandson_in_structure->get_transition_tssb();
	double* sum_child_length = new double[4];
	for(int i = 0;i < 4;i++){
		ithmm->_strength_h *= 10;
		ithmm->_strength_v *= 10;
		ithmm->update_stick_length_of_tssb(htssb, total_stick_length);
		assert(std::abs(htssb->_root->_probability + htssb->_root->_children_stick_length - total_stick_length) < 1e-8);
		double root_length = htssb->_root->_probability;
		sum_child_length[i] = 0;
		for(Node* child: htssb->_root->_children){
			assert(child->_stick_length > 0);
			sum_child_length[i] += child->_stick_length;
			double sum_grandson_length = 0;
			for(Node* grandson: child->_children){
				assert(grandson->_stick_length > 0);
				sum_grandson_length += grandson->_stick_length;
			}
			assert(child->_children_stick_length > sum_grandson_length);
		}
		assert(total_stick_length > root_length + sum_child_length[i]);
	}
	assert(sum_child_length[0] > sum_child_length[1] && sum_child_length[1] > sum_child_length[2] && sum_child_length[2] > sum_child_length[3]);
	delete[] sum_child_length;
}

void test_retrospective_sampling(){
	iTHMM* ithmm = new iTHMM();
	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* root_in_bos = ithmm->_root_in_bos;
	ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(root_in_structure);
	Node* child_in_structure = ithmm->generate_and_add_new_child_to(root_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_structure = ithmm->generate_and_add_new_child_to(child_in_structure);
	Node* grandson_in_htssb = grandson_in_structure->get_myself_in_transition_tssb();
	assert(grandson_in_htssb != NULL);
	Node* grandson_in_child_htssb = child_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_child_htssb != NULL);
	Node* grandson_in_root_htssb = root_in_structure->get_transition_tssb()->find_node_by_tracing_horizontal_indices(grandson_in_htssb);
	assert(grandson_in_root_htssb != NULL);
	ithmm->_lambda_alpha = 0.1;
	ithmm->_lambda_gamma = 1;
	ithmm->_strength_h = 1;
	ithmm->_strength_v = 1;
	bool new_table_generated;

	// vertical
	for(int i = 0;i < 100;i++){
		grandson_in_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->_parent->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 1;i++){
		grandson_in_child_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
		grandson_in_root_htssb->add_customer_to_vertical_crp(1, 0.5, new_table_generated);
	}

	// horizontal
	for(int i = 0;i < 1;i++){
		grandson_in_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(int i = 0;i < 10;i++){
		grandson_in_child_htssb->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
	}
	for(Node* child: grandson_in_root_htssb->_parent->_children){
		for(int i = 0;i < 100;i++){
			child->add_customer_to_horizontal_crp(1, 0.5, new_table_generated);
		}
	}

	double total_stick_length = 1.0;
	for(int depth = 1;depth < 10;depth++){
		ithmm->_depth_limit = depth;
		TSSB* htssb = grandson_in_structure->get_transition_tssb();
		for(int i = 1;i < 1000;i++){
			double uniform = i / 1000.0;
			Node* node = ithmm->retrospective_sampling(uniform, htssb, total_stick_length);
			assert(node->_depth_v <= depth);
		}
		double root_length = htssb->_root->_probability;
		double sum_child_length = 0;
		for(Node* child: htssb->_root->_children){
			assert(child->_stick_length > 0);
			sum_child_length += child->_stick_length;
		}
		assert(total_stick_length > root_length + sum_child_length);
	}
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
	test_add_customer_to_tssb_node();
	cout << "OK" << endl;
	test_remove_customer_from_tssb_node();
	cout << "OK" << endl;
	test_htssb_concentration_vertical();
	cout << "OK" << endl;
	test_htssb_concentration_horizontal();
	cout << "OK" << endl;
	test_nodes_from_root_to_myself();
	cout << "OK" << endl;
	test_compute_expectation_of_vertical_htssb_sbr_ratio();
	cout << "OK" << endl;
	test_compute_expectation_of_horizontal_htssb_sbr_ratio();
	cout << "OK" << endl;
	test_compute_concentration_vertical_htssb_sbr_ratio();
	cout << "OK" << endl;
	test_compute_concentration_horizontal_htssb_sbr_ratio();
	cout << "OK" << endl;
	test_sample_node_in_tssb_by_iterating_node();
	cout << "OK" << endl;
	test_update_stick_length_of_tssb();
	cout << "OK" << endl;
	test_retrospective_sampling();
	cout << "OK" << endl;
	return 0;
}