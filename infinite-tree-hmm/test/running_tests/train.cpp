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

// $1 [vp:139,vs:53,hp:0,hs:190,ref:0][len:1.000000,p:0.261018,ch:0.738982,p:0.261018,sp:0.261018][ow:$1,dv:0,dh:0][][eos:0,other:0]
// 	$2 [vp:0,vs:0,hp:137,hs:0,ref:0][len:0.005351,p:0.005351,ch:0.000000,p:0.005351,sp:0.266369][ow:$1,dv:1,dh:0][0,][eos:0,other:0]
// 	$3 [vp:0,vs:24,hp:113,hs:24,ref:0][len:0.132810,p:0.132810,ch:0.000000,p:0.132810,sp:0.399179][ow:$1,dv:1,dh:1][1,][eos:0,other:0]
// 	$4 [vp:0,vs:68,hp:49,hs:64,ref:0][len:0.342279,p:0.342279,ch:0.000000,p:0.342279,sp:0.741458][ow:$1,dv:1,dh:2][2,][eos:0,other:0]
// 	$8 [vp:0,vs:47,hp:0,hs:49,ref:0][len:0.258036,p:0.258036,ch:0.000000,p:0.258036,sp:0.999493][ow:$1,dv:1,dh:3][3,][eos:0,other:0]
// TSSB[ow:$8,#c:137,#ct:98]
// $1 [vp:39,vs:10,hp:0,hs:49,ref:0][len:1.000000,p:0.241996,ch:0.758004,p:0.241996,sp:0.241996][ow:$8,dv:0,dh:0][][eos:0,other:0]
// 	$2 [vp:0,vs:0,hp:39,hs:0,ref:0][len:0.003830,p:0.003830,ch:0.000000,p:0.003830,sp:0.245826][ow:$8,dv:1,dh:0][0,][eos:0,other:0]
// 	$3 [vp:0,vs:6,hp:33,hs:6,ref:0][len:0.130302,p:0.130302,ch:0.000000,p:0.130302,sp:0.376128][ow:$8,dv:1,dh:1][1,][eos:0,other:0]
// 	$4 [vp:0,vs:24,hp:9,hs:24,ref:0][len:0.385953,p:0.385953,ch:0.000000,p:0.385953,sp:0.762081][ow:$8,dv:1,dh:2][2,][eos:0,other:0]
// 	$8 [vp:0,vs:9,hp:0,hs:9,ref:0][len:0.237556,p:0.237556,ch:0.000000,p:0.237556,sp:0.999637][ow:$8,dv:1,dh:3][3,][eos:0,other:0]

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

int main(){
	// check_distripution();
	std::string filename = "../../../text/test.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	int seed = 0;
	Dataset* dataset = new Dataset(corpus, 0.9, 1, seed);
	Model* model = new Model(dataset, 2);

	model->set_alpha(1);
	model->set_gamma(1);
	model->set_lambda_alpha(0.01);
	model->set_lambda_gamma(0.01);
	model->set_concentration_v(1);
	model->set_concentration_h(1);

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

	model->show_assigned_words_for_each_tag(dictionary, 10);

	for(int i = 0;i < 500000;i++){
		trainer->gibbs();
		trainer->update_hyperparameters();
		cout << "\r" << i << flush;
		if(i % 100 == 0){
			model->show_assigned_words_for_each_tag(dictionary, 10, false);
			cout << trainer->compute_log_p_dataset_train() << endl;
		}
	}
	model->save("ithmm.model");
	return 0;
}