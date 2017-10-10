#include  <iostream>
#include  <string>
#include "../../src/ithmm/sampler.h"
#include "../../src/python/model.h"
#include "../../src/python/dataset.h"
#include "../../src/python/dictionary.h"
#include "../../src/python/trainer.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;

void compare_table(Table* a, Table* b){
	assert(a->_num_customers == b->_num_customers);
	assert(a->_token_id == b->_token_id);
	assert(a->_last_added_index == b->_last_added_index);
	assert(a->_arrangement.size() == b->_arrangement.size());
	for(int i = 0;i < a->_arrangement.size();i++){
		assert(a->_arrangement[i] == b->_arrangement[i]);
	}
}

void compare_hpylm(HPYLM* a, HPYLM* b){
	assert(a->_num_tables == b->_num_tables);
	assert(a->_num_customers == b->_num_customers);
	assert(a->_depth == b->_depth);
	assert(a->_arrangement.size() == b->_arrangement.size());
	for(int i = 0;i < a->_arrangement.size();i++){
		assert(a->_arrangement[i] == b->_arrangement[i]);
	}
}

void compare_node(Node* a, Node* b){
	assert(a->_depth_v == b->_depth_v);
	assert(a->_depth_h == b->_depth_h);
	assert(a->_pass_count_v == b->_pass_count_v);
	assert(a->_stop_count_v == b->_stop_count_v);
	assert(a->_pass_count_h == b->_pass_count_h);
	assert(a->_stop_count_h == b->_stop_count_h);
	assert(a->_num_transitions_to_eos == b->_num_transitions_to_eos);
	assert(a->_num_transitions_to_other == b->_num_transitions_to_other);
	assert(a->_ref_count == b->_ref_count);
	compare_table(a->_table_v, b->_table_v);
	compare_table(a->_table_h, b->_table_h);
	compare_hpylm(a->_hpylm, b->_hpylm);
	for(int v = 0;v < a->_depth_v;v++){
		if(a->_horizontal_indices_from_root[v] != b->_horizontal_indices_from_root[v]){
			a->dump();
			b->dump();
		}
		assert(a->_horizontal_indices_from_root[v] == b->_horizontal_indices_from_root[v]);
	}
	assert(a->has_child() == b->has_child());
	assert(a->_children.size() == b->_children.size());
	for(int i = 0;i < a->_children.size();i++){
		compare_node(a->_children[i], b->_children[i]);
	}
}

void compare(iTHMM* a, iTHMM* b){
	compare_node(a->_root_in_structure, b->_root_in_structure);
}

void check_distripution(){
	iTHMM* ithmm = new iTHMM(sampler::uniform(iTHMM_ALPHA_MIN, iTHMM_ALPHA_MAX),
							sampler::uniform(iTHMM_GAMMA_MIN, iTHMM_GAMMA_MAX),
							sampler::uniform(iTHMM_LAMBDA_ALPHA_MIN, iTHMM_LAMBDA_ALPHA_MAX),
							sampler::uniform(iTHMM_LAMBDA_GAMMA_MAX, iTHMM_LAMBDA_GAMMA_MAX),
							sampler::uniform(ITHMM_SBP_CONCENTRATION_HORIZONTAL_STRENGTH_MIN, ITHMM_SBP_CONCENTRATION_HORIZONTAL_STRENGTH_MAX),
							sampler::uniform(ITHMM_SBP_CONCENTRATION_VERTICAL_STRENGTH_MIN, ITHMM_SBP_CONCENTRATION_VERTICAL_STRENGTH_MAX),
							iTHMM_TAU_0,
							iTHMM_TAU_1,
							1.0 / 2,
							1);

	// ithmm->_conc_v = 1;
	// ithmm->_conc_h = 1;

	Node* root_in_structure = ithmm->_root_in_structure;
	Node* root_in_htssb = ithmm->_root_in_htssb;
	Node* child;
	Node* child_in_htssb;
	TSSB* root_transition_tssb = root_in_structure->get_transition_tssb();

	for(int i = 0;i < 50;i++){
		ithmm->add_customer_to_htssb_node(root_in_htssb);
	}
	
	Node* target = ithmm->generate_and_add_new_child_to(root_in_structure);
	TSSB* target_transition_tssb = target->get_transition_tssb();
	for(int i = 0;i < 10;i++){
		ithmm->add_customer_to_htssb_node(target_transition_tssb->_root);
		ithmm->add_customer_to_tssb_node(root_in_structure);
		ithmm->add_customer_to_hpylm(root_in_structure, 2);
	}
	
	child = ithmm->generate_and_add_new_child_to(root_in_structure);
	child_in_htssb = target_transition_tssb->find_node_by_tracing_horizontal_indices(child);
	for(int i = 0;i < 6;i++){
		ithmm->add_customer_to_htssb_node(child_in_htssb);
		ithmm->add_customer_to_tssb_node(child);
		ithmm->add_customer_to_hpylm(child, 2);
	}

	for(int i = 0;i < 200;i++){
		ithmm->add_customer_to_htssb_node(root_transition_tssb->find_node_by_tracing_horizontal_indices(child));
	}
	
	child = ithmm->generate_and_add_new_child_to(root_in_structure);
	child_in_htssb = target_transition_tssb->find_node_by_tracing_horizontal_indices(child);
	for(int i = 0;i < 40;i++){
		ithmm->add_customer_to_htssb_node(child_in_htssb);
		ithmm->add_customer_to_tssb_node(child);
		ithmm->add_customer_to_hpylm(child, 2);
	}

	for(int i = 0;i < 5;i++){
		ithmm->add_customer_to_htssb_node(root_transition_tssb->find_node_by_tracing_horizontal_indices(child));
	}
	
	child = ithmm->generate_and_add_new_child_to(root_in_structure);
	child_in_htssb = target_transition_tssb->find_node_by_tracing_horizontal_indices(child);
	for(int i = 0;i < 3;i++){
		ithmm->add_customer_to_htssb_node(child_in_htssb);
		ithmm->add_customer_to_tssb_node(child);
		ithmm->add_customer_to_hpylm(child, 2);
	}

	for(int i = 0;i < 40;i++){
		ithmm->add_customer_to_htssb_node(root_transition_tssb->find_node_by_tracing_horizontal_indices(child));
	}

	ithmm->update_stick_length_of_tssb(root_transition_tssb, 1.0);
	root_transition_tssb->dump();
	ithmm->update_stick_length_of_tssb(target_transition_tssb, 1.0);
	target_transition_tssb->dump();
	exit(0);
}

void run_train_loop(){
	std::string filename = "../../../text/test.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	int seed = 0;
	Dataset* dataset = new Dataset(corpus, 0.9, 1, seed);
	Model* model = new Model(dataset, 1, 1, 1, 1, 1, 1, 1, 100, 3);

	// _alpha = sampler::uniform(iTHMM_ALPHA_MIN, iTHMM_ALPHA_MAX);
	// _gamma = sampler::uniform(iTHMM_GAMMA_MIN, iTHMM_GAMMA_MAX);
	// _lambda_alpha = sampler::uniform(iTHMM_LAMBDA_ALPHA_MIN, iTHMM_LAMBDA_ALPHA_MAX);
	// _lambda_gamma = sampler::uniform(iTHMM_LAMBDA_GAMMA_MAX, iTHMM_LAMBDA_GAMMA_MAX);
	// _conc_h = sampler::uniform(ITHMM_SBP_CONCENTRATION_HORIZONTAL_STRENGTH_MIN, ITHMM_SBP_CONCENTRATION_HORIZONTAL_STRENGTH_MAX);
	// _conc_v = sampler::uniform(ITHMM_SBP_CONCENTRATION_VERTICAL_STRENGTH_MIN, ITHMM_SBP_CONCENTRATION_VERTICAL_STRENGTH_MAX);
	// _tau0 = iTHMM_TAU_0;
	// _tau1 = iTHMM_TAU_1;
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("ithmm.dict");
	Trainer* trainer = new Trainer(dataset, model);

	model->show_assigned_words_for_each_tag(dictionary, 10, false);

	for(int i = 0;i < 1000;i++){
		trainer->gibbs();
		trainer->update_hyperparameters();
		cout << "\r" << i << flush;
		if(i % 10 == 0){
			model->show_assigned_words_for_each_tag(dictionary, 10, false);
			cout << trainer->compute_log_p_dataset_train() << endl;
			model->save("ithmm.model");
			Model* _model = new Model("ithmm.model");
			compare(model->_ithmm, _model->_ithmm);
			delete _model;
		}
	}
	model->save("ithmm.model");
	delete corpus;
	delete dataset;
	delete trainer;
	delete model;
}

int main(){
	for(int i = 0;i < 10;i++){
		run_train_loop();
	}
	return 0;
}