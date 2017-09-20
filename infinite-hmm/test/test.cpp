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

void test_tag_word_count(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	for(int tag = 1;tag <= ihmm->get_num_tags();tag++){
		for(int word_id = 0;word_id < 100;word_id++){
			for(int i = 0;i < 100;i++){
				ihmm->_increment_tag_word_count(tag, word_id);
			}
			for(int i = 0;i < 100;i++){
				ihmm->_decrement_tag_word_count(tag, word_id);
			}
			assert(ihmm->_oracle_sum_m_over_q == 0);
			assert(ihmm->_sum_m_i_over_q[tag] == 0);
			assert(ihmm->_m_iq_tables[tag][word_id]->get_num_customers() == 0);
		}
	}
}

void test_tag_bigram_count(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	for(int tag = 1;tag <= ihmm->get_num_tags();tag++){
		for(int context_tag = 1;context_tag <= ihmm->get_num_tags();context_tag++){
			for(int i = 0;i < 100;i++){
				ihmm->_increment_tag_bigram_count(context_tag, tag);
			}
			for(int i = 0;i < 100;i++){
				ihmm->_decrement_tag_bigram_count(context_tag, tag);
			}
			assert(ihmm->_oracle_sum_n_over_j == 0);
			assert(ihmm->_sum_n_i_over_j[tag] == 0);
			assert(ihmm->_n_ij_tables[context_tag][tag]->get_num_customers() == 0);
		}
	}
}

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
		for(int word_id = 0;word_id < ihmm->get_num_words();word_id++){
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
	int num_tags = 1;
	std::string filename = "../../text/test.txt";
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
			cout << model->_hmm->get_num_valid_tags() << endl;
		}
	}
	delete corpus;
	delete dataset;
	delete model;
	delete trainer;
}

void test4(){
	int num_tags = 1;
	std::string filename = "../../text/test.txt";
	Corpus* corpus = new Corpus();
	corpus->add_textfile(filename);
	Dataset* dataset = new Dataset(corpus, 0.9, 1, 0);
	Model* model = new Model(num_tags, dataset);
	model->_hmm->_remove_all_training_dataset(dataset->_word_sequences_train);

	delete corpus;
	delete dataset;
	delete model;
}

void test5(){
	Table* table = new Table();
	bool a = true;
	for(int i = 0;i < 1000;i++){
		table->add_customer(1, a);
	}
	for(int i = 0;i < 1000;i++){
		table->remove_customer(a);
	}
}

void test6(){
	std::vector<Word*> word_vec_1;
	{
		Word* word = new Word();
		word->_id = 0;
		word->_tag = 0;
		word_vec_1.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 1;
		word->_tag = 1;
		word_vec_1.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 2;
		word->_tag = 2;
		word_vec_1.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 0;
		word->_tag = 0;
		word_vec_1.push_back(word);
	}
	std::vector<Word*> word_vec_2;
	{
		Word* word = new Word();
		word->_id = 0;
		word->_tag = 0;
		word_vec_2.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 3;
		word->_tag = 3;
		word_vec_2.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 4;
		word->_tag = 4;
		word_vec_2.push_back(word);
	}
	{
		Word* word = new Word();
		word->_id = 0;
		word->_tag = 0;
		word_vec_2.push_back(word);
	}
	std::vector<std::vector<Word*>> word_sequences;
	for(int i = 0;i < 10;i++){
		word_sequences.push_back(word_vec_1);
	}
	for(int i = 0;i < 5;i++){
		word_sequences.push_back(word_vec_2);
	}

	InfiniteHMM* hmm = new InfiniteHMM(10, 5);
	hmm->initialize_with_training_dataset(word_sequences);

	cout << "oracle:" << endl;
	for(int tag = 0;tag <= hmm->get_num_tags();tag++){
		cout << tag << ":" << hmm->_oracle_n_j_counts[tag] << endl;
	}

	for(int i = 0;i < 2;i++){
		std::vector<Word*> &word_vec = word_sequences[i];
		hmm->perform_gibbs_sampling_with_sequence(word_vec);
		for(int i = 1;i < word_vec.size() - 1;i++){
			cout << word_vec[i]->_tag << endl;
		}
	}
}

int main(){
	// for(int i = 0;i < 10;i++){
	// 	test1();
	// }
	// test2();
	// test3(1000000);
	test_tag_word_count();
	test_tag_bigram_count();
	test6();
	return 0;
}