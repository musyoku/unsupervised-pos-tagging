#include  <iostream>
#include  <string>
#include "../src/python/model.h"
#include "../src/python/dataset.h"
#include "../src/python/dictionary.h"
#include "../src/python/trainer.h"
using namespace bhmm;
using std::cout;
using std::flush;
using std::endl;

int main(){
	int num_tags = 10;
	std::string filename = "../alice.txt";
	Dictionary* dictionary = new Dictionary();
	Dataset* dataset = new Dataset(dictionary);
	dataset->add_textfile(filename, 0.8);
	dictionary->save("bhmm.dict");
	Model* model = new Model();
	boost::python::list py_Wt;
	Trainer* trainer = new Trainer(dataset, model, dictionary);

	for(int i = 0;i < 100;i++){
		trainer->perform_gibbs_sampling();
		trainer->update_hyperparameters();
		cout << "\repoch:" << i << flush;
		if(i % 100 == 0){
			cout << "epoch:" << i << endl;
			// model->_ithmm->_structure_tssb->dump();
			model->show_assigned_words_for_each_tag(dictionary, 20, false);
			cout << "alpha: " << model->_ithmm->_alpha << endl;
			cout << "gamma: " << model->_ithmm->_gamma << endl;
			cout << "lambda_alpha: " << model->_ithmm->_lambda_alpha << endl;
			cout << "lambda_gamma: " << model->_ithmm->_lambda_gamma << endl;
			cout << "strength: " << model->_ithmm->_strength << endl;
			cout << "log_Pdata: " << trainer->compute_log_Pdataset_test() << ", " << trainer->compute_log_Pdataset_train() << endl;
			cout << "PPL: " << trainer->compute_perplexity_test() << ", " << trainer->compute_perplexity_train() << endl;
			cout << "MH: " << model->_ithmm->_num_mh_acceptance / (double)(model->_ithmm->_num_mh_acceptance + model->_ithmm->_num_mh_rejection) << endl;;
			for(int i = 0;i <= model->_ithmm->_current_max_depth;i++){
				cout << "d[" << i << "] = " << model->_ithmm->_hpylm_d_m[i] << endl;
				cout << "theta[" << i << "] = " << model->_ithmm->_hpylm_theta_m[i] << endl;
			}
			model->save("ithmm.model");
		}
	}
	return 0;
}