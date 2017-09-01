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
	std::string filename = "../../text/alice.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.9);
	Dictionary* dictionary = new Dictionary();
	dictionary->load("bhmm.dict");
	Model* model = new Model("bhmm.model");
	std::vector<int> sampled_state_sequence;
	for(auto sentence: dataset->_word_sequences_dev){
		model->viterbi_decode(sentence, sampled_state_sequence);
	}
	return 0;
}