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

	for(int i = 0;i < 1000000;i++){
		trainer->gibbs();
		// trainer->update_hyperparameters();
		cout << "\r" << i << flush;
		if(i % 100 == 0){
			model->show_assigned_words_for_each_tag(dictionary, 10, false);
			cout << trainer->compute_log_p_dataset_train() << endl;
		}
	}
	return 0;
}