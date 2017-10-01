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
	Model* model = new Model(dataset);
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("ithmm.dict");
	Trainer* trainer = new Trainer(dataset, model);

	for(int i = 0;i < 100;i++){
		trainer->perform_gibbs_sampling();
	}
	return 0;
}