#include <cassert>
#include <iostream>
#include "../ithmm/sampler.h"
#include "../ithmm/utils.h"
#include "trainer.h"

namespace ithmm {
	Trainer::Trainer(Dataset* dataset, Model* model){
		_model = model;
		_dict = dataset->_dict;
		_dataset = dataset;
	}
	void Trainer::perform_gibbs_sampling(){
		std::vector<std::vector<Word*>> &dataset = _dataset->_word_sequences_train;
		if(_rand_indices.size() != dataset.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < dataset.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		shuffle(_rand_indices.begin(), _rand_indices.end(), sampler::mt);	// データをシャッフル
		for(int n = 0;n < dataset.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			std::vector<Word*> &word_vec = dataset[data_index];
			_model->_hmm->perform_gibbs_sampling_with_sequence(word_vec);
			
		}
	}
	struct value_comparator {
		bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
			return a.second > b.second;
		}   
	};
	void Trainer::_before_viterbi_decode(){
		assert(_dataset->_max_num_words_in_line > 0);
		_decode_table = new double*[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_decode_table[i] = new double[_model->_hmm->get_num_tags() + 1];
		}
	}
	void Trainer::_after_viterbi_decode(){
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			delete[] _decode_table[i];
		}
		delete[] _decode_table;
	}
	void Trainer::_before_compute_log_p_dataset(){
		// 計算用のテーブルを確保
		assert(_dataset->_max_num_words_in_line > 0);
		_forward_table = new double*[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_forward_table[i] = new double[_model->_hmm->get_num_tags() + 1];
		}
	}
	void Trainer::_after_compute_log_p_dataset(){
		// 計算用のテーブルを解放
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			delete[] _forward_table[i];
		}
		delete[] _forward_table;
	}
	// データセット全体の対数尤度を計算
	double Trainer::compute_log_p_dataset_train(){
		return _compute_log_p_dataset(_dataset->_word_sequences_train);
	}
	double Trainer::compute_log_p_dataset_dev(){
		return _compute_log_p_dataset(_dataset->_word_sequences_dev);
	}
	double Trainer::_compute_log_p_dataset(std::vector<std::vector<Word*>> &dataset){
		_before_compute_log_p_dataset();
		// データごとの対数尤度を足していく
		double log_p_dataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			std::vector<Word*> &sentence = dataset[data_index];
			double p_x = _model->compute_p_sentence(sentence, _forward_table);
			if(p_x > 0){
				log_p_dataset += log(p_x);
			}
		}
		_after_compute_log_p_dataset();
		return log_p_dataset;
	}
	void Trainer::set_model(Model* model){
		_model = model;
	}
}