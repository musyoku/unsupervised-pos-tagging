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

void run_train_loop(){
	std::string filename = "../../../text/test.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	int seed = 0;
	Dataset* dataset = new Dataset(corpus, 0.9, 0, seed);
	Model* model = new Model(dataset, 1, 1, 1, 1, 1, 1, 1, 100, 1);

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
		}
	}
	model->save("ithmm.model");
	delete corpus;
	delete dataset;
	delete trainer;
	delete model;
}

int main(){
	run_train_loop();
	return 0;
}