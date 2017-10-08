#include <iostream>
#include <fstream>
#include <unordered_set>
#include "dataset.h"
#include "../ithmm/utils.h"
#include "../ithmm/sampler.h"

namespace ithmm {
	Dataset::Dataset(Corpus* corpus, double train_split, int unknown_count, int seed){
		_dict = new Dictionary();
		_max_num_words_in_line = corpus->_max_num_words_in_line;
		_min_num_words_in_line = corpus->_min_num_words_in_line;

		std::vector<int> rand_indices;
		for(int i = 0;i < corpus->_word_sequences.size();i++){
			rand_indices.push_back(i);
		}
		sampler::set_seed(seed);
		shuffle(rand_indices.begin(), rand_indices.end(), sampler::mt);	// データをシャッフル
		train_split = std::min(1.0, std::max(0.0, train_split));
		int num_train_data = corpus->_word_sequences.size() * train_split;
		for(int i = 0;i < rand_indices.size();i++){
			std::vector<std::wstring> &word_str_vec = corpus->_word_sequences[rand_indices[i]];
			if(i < num_train_data){
				_add_words_to_dataset(word_str_vec, _word_sequences_train, corpus, unknown_count);
			}else{
				_add_words_to_dataset(word_str_vec, _word_sequences_dev, corpus, unknown_count);
			}
		}
	}
	Dataset::~Dataset(){
		for(int n = 0;n < _word_sequences_train.size();n++){
			std::vector<Word*> &words = _word_sequences_train[n];
			for(int m = 0;m < words.size();m++){
				Word* word = words[m];
				delete word;
			}
		}
		for(int n = 0;n < _word_sequences_dev.size();n++){
			std::vector<Word*> &words = _word_sequences_dev[n];
			for(int m = 0;m < words.size();m++){
				Word* word = words[m];
				delete word;
			}
		}
		delete _dict;
	}
	void Dataset::_add_words_to_dataset(std::vector<std::wstring> &word_str_vec, std::vector<std::vector<Word*>> &dataset, Corpus* corpus, int unknown_count){
		assert(word_str_vec.size() > 0);
		std::vector<Word*> words;
		// 単語列
		for(auto word_str: word_str_vec){
			if(word_str.size() == 0){
				continue;
			}
			Word* word = new Word();
			int count = corpus->get_count_of_word(word_str);
			if(count <= unknown_count){
				word->_id = ID_UNK;
			}else{
				word->_id = _dict->add_word_string(word_str);
				_word_count[word->_id] += 1;
			}
			word->_state = NULL;
			words.push_back(word);
		}
		// 追加
		dataset.push_back(words);

		if((int)words.size() > _max_num_words_in_line){
			_max_num_words_in_line = words.size();
		}
		if((int)words.size() < _min_num_words_in_line || _min_num_words_in_line == -1){
			_min_num_words_in_line = words.size();
		}
	}
	int Dataset::get_num_words(){
		return _dict->get_vocabrary_size();
	}
	Dictionary &Dataset::get_dict_obj(){
		return *_dict;
	}
}