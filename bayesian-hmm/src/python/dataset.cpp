#include <fstream>
#include "dataset.h"
#include "../bhmm/utils.h"
#include "../bhmm/sampler.h"

namespace bhmm {
	Dataset::Dataset(Dictionary* dict){
		_dict = dict;
		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
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
	}
	void Dataset::add_textfile(std::string filename, double train_split_ratio){
		std::wifstream ifs(filename.c_str());
		std::wstring sentence_str;
		assert(ifs.fail() == false);
		std::vector<std::wstring> sentence_vec;
		while (getline(ifs, sentence_str) && !sentence_str.empty()){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			sentence_vec.push_back(sentence_str);
		}
		train_split_ratio = std::min(1.0, std::max(0.0, train_split_ratio));
		int train_split = sentence_vec.size() * train_split_ratio;
		std::vector<int> rand_indices;
		for(int i = 0;i < sentence_vec.size();i++){
			rand_indices.push_back(i);
		}
		shuffle(rand_indices.begin(), rand_indices.end(), sampler::mt);	// データをシャッフル
		for(int i = 0;i < rand_indices.size();i++){
			std::wstring &sentence_str = sentence_vec[rand_indices[i]];
			if(i < train_split){
				add_sentence_str_train(sentence_str);
			}else{
				add_sentence_str_dev(sentence_str);
			}
		}
	}
	void Dataset::add_sentence_str_train(std::wstring sentence_str){
		std::vector<std::wstring> word_str_vec;
		utils::split_word_by(sentence_str, L' ', word_str_vec);	// スペースで分割
		_add_words_to_dataset(word_str_vec, _word_sequences_train);
	}
	void Dataset::add_sentence_str_dev(std::wstring sentence_str){
		std::vector<std::wstring> word_str_vec;
		utils::split_word_by(sentence_str, L' ', word_str_vec);	// スペースで分割
		_add_words_to_dataset(word_str_vec, _word_sequences_dev);
	}
	void Dataset::_before_python_add_sentence_str(boost::python::list &py_word_str_list, std::vector<std::wstring> &word_str_vec){
		int num_words = boost::python::len(py_word_str_list);
		for(int i = 0;i < num_words;i++){
			std::wstring word = boost::python::extract<std::wstring>(py_word_str_list[i]);
			word_str_vec.push_back(word);
		}
	}
	void Dataset::python_add_words_train(boost::python::list py_word_str_list){
		std::vector<std::wstring> word_str_vec;
		_before_python_add_sentence_str(py_word_str_list, word_str_vec);
		_add_words_to_dataset(word_str_vec, _word_sequences_train);
	}
	void Dataset::python_add_words_dev(boost::python::list py_word_str_list){
		std::vector<std::wstring> word_str_vec;
		_before_python_add_sentence_str(py_word_str_list, word_str_vec);
		_add_words_to_dataset(word_str_vec, _word_sequences_dev);
	}
	void Dataset::_add_words_to_dataset(std::vector<std::wstring> &word_str_vec, std::vector<std::vector<Word*>> &dataset){
		assert(word_str_vec.size() > 0);
		std::vector<Word*> words;

		for(auto word_str: word_str_vec){
			if(word_str.size() == 0){
				continue;
			}
			Word* word = new Word();
			word->_id = _dict->add_word_string(word_str);
			word->_state = 0;
			words.push_back(word);
			_word_count[word->_id] += 1;
		}

		Word* eos = new Word();
		eos->_id = ID_EOS;
		eos->_state = 0;
		words.push_back(eos);
		_word_count[ID_EOS] += 1;

		dataset.push_back(words);

		if((int)words.size() > _max_num_words_in_line){
			_max_num_words_in_line = words.size();
		}
		if((int)words.size() < _min_num_words_in_line || _min_num_words_in_line == -1){
			_min_num_words_in_line = words.size();
		}
	}
	void Dataset::mark_low_frequency_words_as_unknown(int threshold){
		_mark_low_frequency_words_as_unknown(threshold, _word_sequences_train);
		_mark_low_frequency_words_as_unknown(threshold, _word_sequences_dev);
	}
	void Dataset::_mark_low_frequency_words_as_unknown(int threshold, std::vector<std::vector<Word*>> &word_sequence_vec){
		for(int data_index = 0;data_index < word_sequence_vec.size();data_index++){
			std::vector<Word*> &data = word_sequence_vec[data_index];
			for(auto word = data.begin(), end = data.end();word != end;word++){
				id word_id = (*word)->_id;
				int count = get_count_of_word(word_id);
				if(count <= threshold){
					(*word)->_id = ID_UNK;
				}
			}
		}
	}
	int Dataset::get_num_words(){
		return _word_count.size();
	}
	int Dataset::get_count_of_word(id word_id){
		auto itr = _word_count.find(word_id);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
}