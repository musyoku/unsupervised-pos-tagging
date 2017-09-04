#include <cassert>
#include <iostream>
#include "../bhmm/sampler.h"
#include "../bhmm/utils.h"
#include "trainer.h"

namespace bhmm {
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
	void Trainer::update_hyperparameters(){
		sample_new_alpha();
		sample_new_beta();
	}
	void Trainer::sample_new_alpha(){
		HMM* hmm = _model->_hmm;
		double old_alpha = hmm->_alpha;
		double old_log_p_x = compute_log_p_dataset_train();

		double new_alpha = sampler::normal(old_alpha, 0.1 * old_alpha);
		hmm->_alpha = new_alpha;
		double new_log_p_x = compute_log_p_dataset_train();
		
		double sigma_old_alpha = 0.1 * old_alpha;
		double sigma_new_alpha = 0.1 * new_alpha;
		double var_old_alpha = sigma_old_alpha * sigma_old_alpha;
		double var_new_alpha = sigma_new_alpha * sigma_new_alpha;

		double correcting_term = (old_alpha / new_alpha) * exp(
			  0.5 * (new_alpha - old_alpha) * (new_alpha - old_alpha) / var_old_alpha
			+ 0.5 * (old_alpha - new_alpha) * (old_alpha - new_alpha) / var_new_alpha
		);

		// 採択率
		double adoption_rate = std::min(1.0, exp(new_log_p_x - old_log_p_x) * correcting_term);
		double bernoulli = sampler::uniform(0, 1);
		if(bernoulli <= adoption_rate){
			hmm->_alpha = new_alpha;
		}else{
			hmm->_alpha = old_alpha;
		}
	}
	void Trainer::sample_new_beta(){
		HMM* hmm = _model->_hmm;
		for(int tag = 1;tag <= hmm->_num_tags;tag++){
			double old_beta = _model->_hmm->_beta[tag];
			double old_log_p_x = compute_log_p_dataset_train();

			double new_beta = sampler::normal(old_beta, 0.1 * old_beta);
			_model->_hmm->_beta[tag] = new_beta;
			double new_log_p_x = compute_log_p_dataset_train();
			
			double sigma_old_beta = 0.1 * old_beta;
			double sigma_new_beta = 0.1 * new_beta;
			double var_old_beta = sigma_old_beta * sigma_old_beta;
			double var_new_beta = sigma_new_beta * sigma_new_beta;

			double correcting_term = (old_beta / new_beta) * exp(
				  0.5 * (new_beta - old_beta) * (new_beta - old_beta) / var_old_beta
				+ 0.5 * (old_beta - new_beta) * (old_beta - new_beta) / var_new_beta
			);

			// 採択率
			double adoption_rate = std::min(1.0, exp(new_log_p_x - old_log_p_x) * correcting_term);
			double bernoulli = sampler::uniform(0, 1);
			if(bernoulli <= adoption_rate){
				_model->_hmm->_beta[tag] = new_beta;
			}else{
				_model->_hmm->_beta[tag] = old_beta;
			}
		}
	}
	struct value_comparator {
		bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
			return a.second > b.second;
		}   
	};
	boost::python::list Trainer::python_get_all_words_of_each_tag(int threshold){
		std::vector<boost::python::list> result;
		HMM* hmm = _model->_hmm;
		for(int tag = 1;tag <= hmm->_num_tags;tag++){
			std::vector<boost::python::tuple> words;
			std::multiset<std::pair<int, int>, value_comparator> ranking;
			for(id word_id = 0;word_id < hmm->_num_words;word_id++){
				int count = hmm->_tag_word_counts[tag][word_id];
				if(count > 0){
					ranking.insert(std::make_pair(word_id, count));
				}
			}
			for(auto elem: ranking){
				if(elem.second <= threshold){
					continue;
				}
				std::wstring word = _dict->word_id_to_string(elem.first);
				words.push_back(boost::python::make_tuple(word, elem.second));
			}
			result.push_back(utils::list_from_vector(words));
		}
		return utils::list_from_vector(result);
	}
	void Trainer::show_typical_words_of_each_tag(int number_to_show){
		_model->show_typical_words_of_each_tag(number_to_show, _dict);
	}
	void Trainer::_before_viterbi_decode(){
		assert(_dataset->_max_num_words_in_line > 0);
		_decode_table = new double**[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_decode_table[i] = new double*[_model->_hmm->_num_tags + 1];
			for(int k = 0;k <=_model->_hmm->_num_tags;k++){
				_decode_table[i][k] = new double[_model->_hmm->_num_tags + 1];
			}
		}
	}
	void Trainer::_after_viterbi_decode(){
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			for(int k = 0;k <= _model->_hmm->_num_tags;k++){
				delete[] _decode_table[i][k];
			}
			delete[] _decode_table[i];
		}
		delete[] _decode_table;
	}
	void Trainer::_before_compute_log_p_dataset(){
		// 計算用のテーブルを確保
		assert(_dataset->_max_num_words_in_line > 0);
		_forward_table = new double**[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_forward_table[i] = new double*[_model->_hmm->_num_tags + 1];
			for(int k = 0;k <= _model->_hmm->_num_tags;k++){
				_forward_table[i][k] = new double[_model->_hmm->_num_tags + 1];
			}
		}
	}
	void Trainer::_after_compute_log_p_dataset(){
		// 計算用のテーブルを解放
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			for(int k = 0;k <= _model->_hmm->_num_tags;k++){
				delete[] _forward_table[i][k];
			}
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
	void Trainer::anneal_temperature(double temperature){
		_model->_hmm->anneal_temperature(temperature);
	}
}