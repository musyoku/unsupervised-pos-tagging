#include  <iostream>
#include  <string>
#include "../src/bhmm/utils.h"
#include "../src/python/model.h"
#include "../src/python/dataset.h"
#include "../src/python/dictionary.h"
#include "../src/python/trainer.h"
using namespace bhmm;
using std::wcout;
using std::endl;

int main(){
	std::string filename = "../../text/alice.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	Dataset* dataset = new Dataset(corpus, 0.9, 1);
	Dictionary* dictionary = new Dictionary();
	dictionary->load("bhmm.dict");
	Model* model = new Model("bhmm.model");
	std::vector<int> sampled_state_sequence;
	for(auto sentence: dataset->_word_sequences_train){
		model->viterbi_decode(sentence, sampled_state_sequence);
		for(int i = 0;i < sentence.size();i++){
			wcout << dictionary->word_id_to_string(sentence[i]->_id) << ", " << sampled_state_sequence[i] << endl;
		}
	}
	return 0;
}