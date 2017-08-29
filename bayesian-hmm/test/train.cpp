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
	int num_tags = 10;
	std::string filename = "../../text/ptb.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("bhmm.dict");
	Model* model = new Model(num_tags);
	std::vector<int> Wt;
	for(int p = 0;p < num_tags;p++){
		Wt.push_back(1000);
	}
	Trainer* trainer = new Trainer(dataset, model, Wt);

	for(int i = 1;i <= 10000;i++){
		trainer->perform_gibbs_sampling();
		// trainer->update_hyperparameters();
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