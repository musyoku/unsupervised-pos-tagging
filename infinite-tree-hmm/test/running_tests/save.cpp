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
	double p_dataset = trainer->compute_log_p_dataset_train();
	cout << model->_ithmm->_structure_tssb->_num_customers << endl;
	delete model;

	model = new Model("ithmm.model");
	trainer->_model = model;
	double _p_dataset = trainer->compute_log_p_dataset_train();
	cout << model->_ithmm->_structure_tssb->_num_customers << endl;
	cout << p_dataset << endl;
	cout << _p_dataset << endl;

	delete model;
	delete trainer;
	delete corpus;
	delete dataset;
	return 0;
}