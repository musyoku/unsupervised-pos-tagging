#include  <iostream>
#include  <string>
#include "../src/bhmm/utils.h"
#include "../src/python/model.h"
#include "../src/python/dataset.h"
#include "../src/python/dictionary.h"
#include "../src/python/trainer.h"
using namespace bhmm;
using std::cout;
using std::flush;
using std::endl;

int main(){
	int num_tags = 20;
	std::string filename = "../../text/ptb.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("bhmm.dict");
	Model* model = new Model(num_tags);
	model->set_temperature(1.5);
	model->set_minimum_temperature(0.08);
	std::vector<int> Wt;
	for(int p = 0;p < num_tags;p++){
		Wt.push_back(1000);
	}
	Trainer* trainer = new Trainer(dataset, model, Wt);

	for(int i = 1;i <= 10000;i++){
		trainer->perform_gibbs_sampling();
		model->anneal_temperature(0.99989);
		trainer->update_hyperparameters();
		if(i % 100 == 0){
			cout << trainer->compute_log_p_dataset_train() << endl;
			model->save("bhmm.model");
		}
		if(i % 1000 == 0){
			trainer->show_typical_words_of_each_tag(10);
		}

	}
	return 0;
}