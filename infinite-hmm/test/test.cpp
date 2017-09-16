#include  <iostream>
#include  <string>
#include "../src/ihmm/ihmm.h"
#include "../src/ihmm/utils.h"
#include "../src/python/model.h"
#include "../src/python/corpus.h"
#include "../src/python/dataset.h"
#include "../src/python/dictionary.h"
#include "../src/python/trainer.h"
using namespace ihmm;
using namespace std;

void test1(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	int num_tags = ihmm->get_num_tags();
	for(int tag = num_tags;tag >= 1;tag--){
		ihmm->_delete_tag(tag);
	}
	for(int i = 0;i < 100;i++){
		ihmm->_add_new_tag();
		ihmm->_delete_tag(1);
	}
	cout << ihmm->get_num_tags() << endl;
	for(int i = 0;i < 10;i++){
		ihmm->_add_new_tag();
	}
	num_tags = ihmm->get_num_tags();
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_increment_tag_bigram_count(context_tag, tag);
			}
		}
	}
	for(int context_tag = 1;context_tag <= num_tags + 1;context_tag++){
		for(int tag = 1;tag <= num_tags + 1;tag++){
			double prob = ihmm->compute_p_tag_given_context(tag, context_tag);
			cout << "p(" << tag << "|" << context_tag << ") = " << prob << endl;
		}
	}
	for(int tag = 1;tag <= num_tags + 1;tag++){
		for(id word_id = 0;word_id < ihmm->get_num_words();word_id++){
			double prob = ihmm->compute_p_word_given_tag(word_id, tag);
			cout << "p(" << word_id << "|" << tag << ") = " << prob << endl;
		}
	}
	num_tags = ihmm->get_num_tags();
	for(int tag = 0;tag <= num_tags;tag++){
		cout << "tag: " << tag << endl;
		cout << ihmm->_sum_n_i_over_j[tag] << endl;
		cout << ihmm->_oracle_n_j_counts[tag] << endl;
	}
	cout << "oracle: " << ihmm->_oracle_sum_n_over_j << endl;
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_decrement_tag_bigram_count(context_tag, tag);
			}
		}
	}
	num_tags = ihmm->get_num_tags();
	for(int tag = 0;tag <= num_tags;tag++){
		cout << "tag: " << tag << endl;
		cout << ihmm->_sum_n_i_over_j[tag] << endl;
		cout << ihmm->_oracle_n_j_counts[tag] << endl;
	}
	for(int tag = num_tags;tag >= 1;tag--){
		ihmm->_delete_tag(tag);
	}
	cout << ihmm->_sum_n_i_over_j.size() << endl;
	cout << ihmm->_oracle_n_j_counts.size() << endl;
	cout << ihmm->_sum_m_i_over_q.size() << endl;
	cout << "oracle: " << ihmm->_oracle_sum_n_over_j << endl;
	delete ihmm;
}

void test2(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	ihmm->_perform_gibbs_sampling_on_markov_blanket(1, 2, 0);
}

void test3(int num_iterations){
	int num_tags = 10;
	std::string filename = "../../text/alice.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	Dataset* dataset = new Dataset(corpus, 0.9, 1, 0);
	Dictionary* dictionary = dataset->_dict;
	dictionary->save("ihmm.dict");
	Model* model = new Model(num_tags, dataset);
	Trainer* trainer = new Trainer(dataset, model);

	for(int i = 1;i <= num_iterations;i++){
		trainer->perform_gibbs_sampling();
		cout << "\r" << i << flush;
		if(i % 100 == 0){
			cout << "\r" << flush;
			cout << trainer->compute_log_p_dataset_train() << ", " << trainer->compute_log_p_dataset_dev() << endl;
			model->save("ihmm.model");
			cout << model->_hmm->get_num_tags() << endl;
		}
	}
	delete corpus;
	delete dataset;
	delete dictionary;
	delete model;
	delete trainer;
}

int main(){
	for(int i = 0;i < 10;i++){
		test1();
	}
	test2();
	test3(1000);
	return 0;
}