#include <iostream>
#include <fstream>
#include <unordered_set>
#include "corpus.h"
#include "../ithmm/utils.h"
#include "../ithmm/sampler.h"

namespace ithmm {
	Corpus::Corpus(){
		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
	}
	void Corpus::add_textfile(std::string filename){
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
		for(std::wstring &sentence_str: sentence_vec){
			std::vector<std::wstring> word_str_vec;
			utils::split_word_by(sentence_str, L' ', word_str_vec);	// スペースで分割
			_add_words_to_corpus(word_str_vec);
		}
	}
	void Corpus::_before_python_add_sentence_str(boost::python::list &py_word_str_list, std::vector<std::wstring> &word_str_vec){
		int num_words = boost::python::len(py_word_str_list);
		for(int i = 0;i < num_words;i++){
			std::wstring word = boost::python::extract<std::wstring>(py_word_str_list[i]);
			word_str_vec.push_back(word);
		}
	}
	void Corpus::python_add_words(boost::python::list py_word_str_list){
		std::vector<std::wstring> word_str_vec;
		_before_python_add_sentence_str(py_word_str_list, word_str_vec);
		assert(word_str_vec.size() > 0);
		_add_words_to_corpus(word_str_vec);
	}
	void Corpus::_add_words_to_corpus(std::vector<std::wstring> &word_str_vec){
		assert(word_str_vec.size() > 0);
		// 頻度をカウント
		for(std::wstring word_str: word_str_vec){
			auto itr = _word_count.find(word_str);
			if(itr == _word_count.end()){
				_word_count[word_str] = 1;
				continue;
			}
			_word_count[word_str] += 1;
		}
		// コーパスに追加
		_word_sequences.push_back(word_str_vec);
		// 行あたりの最大単語数を更新
		if((int)word_str_vec.size() > _max_num_words_in_line){
			_max_num_words_in_line = word_str_vec.size();
		}
		if((int)word_str_vec.size() < _min_num_words_in_line || _min_num_words_in_line == -1){
			_min_num_words_in_line = word_str_vec.size();
		}
	}
	int Corpus::get_num_words(){
		return _word_count.size();
	}
	int Corpus::get_count_of_word(std::wstring word_str){
		auto itr = _word_count.find(word_str);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
}