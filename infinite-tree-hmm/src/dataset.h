#pragma once
#include "dictionary.h"

class Dataset{
public:
	Dictionary* _dict;
	std::unordered_map<id, int> _word_count;
	std::vector<std::vector<Word*>> _word_sequences_train;
	std::vector<std::vector<Word*>> _word_sequences_test;
	std::vector<int> _rand_indices;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
	Dataset(Dictionary* dict){
		_dict = dict;
		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
	}
	~Dataset(){
		for(int n = 0;n < _word_sequences_train.size();n++){
			std::vector<Word*> &words = _word_sequences_train[n];
			for(int m = 0;m < words.size();m++){
				Word* word = words[m];
				delete word;
			}
		}
		for(int n = 0;n < _word_sequences_test.size();n++){
			std::vector<Word*> &words = _word_sequences_test[n];
			for(int m = 0;m < words.size();m++){
				Word* word = words[m];
				delete word;
			}
		}
	}
	void add_textfile(std::string filename, double train_split_ratio){
		// c_printf("[*]%s\n", (boost::format("%sを読み込んでいます ...") % filename.c_str()).str().c_str());
		std::wifstream ifs(filename.c_str());
		std::wstring line_str;
		if (ifs.fail()){
			c_printf("[R]%s [*]%s", "エラー", (boost::format("%sを開けません.") % filename.c_str()).str().c_str());
			exit(1);
		}
		std::vector<std::wstring> lines;
		while (getline(ifs, line_str) && !line_str.empty()){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			lines.push_back(line_str);
		}
		train_split_ratio = std::min(1.0, std::max(0.0, train_split_ratio));
		int train_split = lines.size() * train_split_ratio;
		std::vector<int> rand_indices;
		for(int i = 0;i < lines.size();i++){
			rand_indices.push_back(i);
		}
		shuffle(rand_indices.begin(), rand_indices.end(), Sampler::mt);	// データをシャッフル
		for(int i = 0;i < rand_indices.size();i++){
			std::wstring &line_str = lines[rand_indices[i]];
			if(i < train_split){
				add_train_data(line_str);
			}else{
				add_test_data(line_str);
			}
		}
		// std::cout << "train: " << _word_sequences_train.size() << std::endl;
		// std::cout << "test:  " << _word_sequences_test.size() << std::endl;
		// c_printf("[*]%s\n", (boost::format("%sを読み込みました.") % filename.c_str()).str().c_str());
	}
	void add_train_data(std::wstring line_str){
		_add_data_to(line_str, _word_sequences_train);
	}
	void add_test_data(std::wstring line_str){
		_add_data_to(line_str, _word_sequences_test);
	}
	void _add_data_to(std::wstring &line_str, std::vector<std::vector<Word*>> &dataset){
		std::vector<std::wstring> word_strs;
		split_word_by(line_str, L' ', word_strs);	// スペースで分割
		if(word_strs.size() > 0){
			std::vector<Word*> words;

			for(auto word_str: word_strs){
				if(word_str.size() == 0){
					continue;
				}
				Word* word = new Word();
				word->_id = _dict->add_string(word_str);
				word->_state = NULL;
				words.push_back(word);
				_word_count[word->_id] += 1;
			}

			Word* eos = new Word();
			eos->_id = ID_EOS;
			eos->_state = NULL;
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
	}
	void mark_low_frequency_words_as_unknown(int threshold = 1){
		for(int data_index = 0;data_index < _word_sequences_train.size();data_index++){
			std::vector<Word*> &data = _word_sequences_train[data_index];
			for(auto word = data.begin(), end = data.end();word != end;word++){
				id word_id = (*word)->_id;
				int count = get_count_for_word(word_id);
				if(count <= threshold){
					(*word)->_id = ID_UNK;
				}
			}
		}
	}
	int get_num_words(){
		return _word_count.size();
	}
	int get_count_for_word(id word_id){
		auto itr = _word_count.find(word_id);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
};