#include "ihmm.h"

namespace ihmm {
	InfiniteHMM::InfiniteHMM(){
		_initial_num_tags = 0;
		_num_words = 0;
		_alpha = 1;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_sum_oracle_words_count = 0;
		_sum_oracle_tags_count = 0;
	}
	InfiniteHMM::InfiniteHMM(int initial_num_tags, int num_words){
		InfiniteHMM();
		_initial_num_tags = initial_num_tags;
		_num_words = num_words;
	}
	InfiniteHMM::~InfiniteHMM(){
		for(int tag: _bigram_tag_table){
			unordered_map<int, Table*> &tables = _bigram_tag_table[tag];
			for(auto &table: tables){
				delete table.second;
			}
		}
		for(int tag: _tag_word_table){
			unordered_map<int, Table*> &tables = _bigram_tag_table[tag];
			for(auto &table: tables){
				delete table.second;
			}
		}
	}
	void initialize_with_training_corpus(std::vector<std::vector<Word*>> &dataset){
		
	}
}