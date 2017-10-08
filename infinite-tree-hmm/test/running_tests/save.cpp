#include  <iostream>
#include  <string>
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
		assert(a->_horizontal_indices_from_root[v] == b->_horizontal_indices_from_root[v]);
	}
	assert(a->has_child() == b->has_child());
	if(a->has_child()){
		assert(a->_children.size() == b->_children.size());
		for(int i = 0;i < a->_children.size();i++){
			compare_node(a->_children[i], b->_children[i]);
		}
	}
}

void compare(iTHMM* a, iTHMM* b){
	compare_node(a->_root_in_structure, b->_root_in_structure);
}

int main(){
	std::string filename = "../../../text/test.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	Dataset* dataset = new Dataset(corpus, 0.9, 1);
	Model* model = new Model(dataset, 1);
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("ithmm.dict");
	Trainer* trainer = new Trainer(dataset, model);

	model->show_assigned_words_for_each_tag(dictionary, 10);

	for(int i = 0;i < 1000;i++){
		trainer->gibbs();
		// trainer->update_hyperparameters();
		model->show_assigned_words_for_each_tag(dictionary, 10, false);
	}
	model->save("ithmm.model");
	double p_dataset_1 = trainer->compute_log_p_dataset_train();
	double p_dataset_2 = trainer->compute_log_p_dataset_train();

	cout << "loading ..." << endl;
	Model* _model = new Model("ithmm.model");

	cout << model->_ithmm->_structure_tssb->_num_customers << endl;
	cout << _model->_ithmm->_structure_tssb->_num_customers << endl;
	model->_ithmm->_root_in_structure->dump();
	_model->_ithmm->_root_in_structure->dump();
	compare(model->_ithmm, _model->_ithmm);

	cout << p_dataset_1 << endl;
	cout << p_dataset_2 << endl;

	trainer->set_model(model);
	double _p_dataset = trainer->compute_log_p_dataset_train();

	cout << _p_dataset << endl;
	delete model;
	delete trainer;
	delete corpus;
	delete dataset;
	return 0;
}