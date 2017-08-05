#include <cassert>
#include "trainer.h"

namespace ithmm {
	Trainer::Trainer(Dataset* dataset, Model* model, Dictionary* dict){
		_dataset = dataset;
		_model = model;
		_dict = dict;
		_forward_table = NULL;
		_decode_table = NULL;
		_model->_ithmm->set_word_g0(1.0 / _dataset->_word_count.size());
		_model->_ithmm->initialize_with_training_dataset(_dataset->_word_sequences_train);
	}
	void Trainer::remove_all_data(){
		_model->_ithmm->remove_all_data(_dataset->_word_sequences_train);
		_model->_ithmm->delete_invalid_children();
	}
	void Trainer::perform_gibbs_sampling(){
		if(_rand_indices.size() != _dataset->_word_sequences_train.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < _dataset->_word_sequences_train.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		_model->_ithmm->_num_mh_acceptance = 0;
		_model->_ithmm->_num_mh_rejection = 0;
		shuffle(_rand_indices.begin(), _rand_indices.end(), sampler::mt);	// データをシャッフル
		for(int n = 0;n < _dataset->_word_sequences_train.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			std::vector<Word*> &data = _dataset->_word_sequences_train[data_index];
			_model->_ithmm->perform_gibbs_sampling_data(data);
		}
		_model->_ithmm->delete_invalid_children();
	}
	void Trainer::_before_viterbi_decode(std::vector<Node*> &nodes){
		_before_compute_log_Pdataset(nodes);
		assert(_dataset->_max_num_words_in_line > 0);
		_decode_table = new double*[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_decode_table[i] = new double[nodes.size()];
		}
	}
	void Trainer::_after_viterbi_decode(){
		_after_compute_log_Pdataset();
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			delete[] _decode_table[i];
		}
		delete[] _decode_table;
	}
	void Trainer::_before_compute_log_Pdataset(std::vector<Node*> &nodes){
		// あらかじめ全HTSSBの棒の長さを計算しておく
		_model->enumerate_all_states(nodes);
		_model->precompute_all_stick_lengths(nodes);
		// 計算用のテーブルを確保
		assert(_dataset->_max_num_words_in_line > 0);
		_forward_table = new double*[_dataset->_max_num_words_in_line];
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			_forward_table[i] = new double[nodes.size()];
		}
	}
	void Trainer::_after_compute_log_Pdataset(){
		// 計算用のテーブルを解放
		assert(_dataset->_max_num_words_in_line > 0);
		for(int i = 0;i < _dataset->_max_num_words_in_line;i++){
			delete[] _forward_table[i];
		}
		delete[] _forward_table;
	}
	// データセット全体の対数尤度を計算
	double Trainer::compute_log_Pdataset_train(){
		return _compute_log_Pdataset(_dataset->_word_sequences_train);
	}
	double Trainer::compute_log_Pdataset_dev(){
		return _compute_log_Pdataset(_dataset->_word_sequences_dev);
	}
	double Trainer::_compute_log_Pdataset(std::vector<std::vector<Word*>> &dataset){
		std::vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			std::vector<Word*> &data = dataset[data_index];
			double Px = _model->compute_Pdata(data, nodes, _forward_table);
			if(Px > 0){
				log_Pdataset += log(Px);
			}
		}
		_after_compute_log_Pdataset();
		return log_Pdataset;
	}
	double Trainer::compute_log2_Pdataset_train(){
		return _compute_log2_Pdataset(_dataset->_word_sequences_train);
	}
	double Trainer::compute_log2_Pdataset_dev(){
		return _compute_log2_Pdataset(_dataset->_word_sequences_dev);
	}
	double Trainer::_compute_log2_Pdataset(std::vector<std::vector<Word*>> &dataset){
		std::vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			std::vector<Word*> &data = dataset[data_index];
			double Px = _model->compute_Pdata(data, nodes, _forward_table);
			if(Px > 0){
				log_Pdataset += log2(Px);
			}
		}
		_after_compute_log_Pdataset();
		return log_Pdataset;
	}
	double Trainer::compute_perplexity_train(){
		return _compute_perplexity(_dataset->_word_sequences_train);
	}
	double Trainer::compute_perplexity_dev(){
		return _compute_perplexity(_dataset->_word_sequences_dev);
	}
	double Trainer::_compute_perplexity(std::vector<std::vector<Word*>> &dataset){
		std::vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			std::vector<Word*> &data = dataset[data_index];
			double Px = _model->compute_Pdata(data, nodes, _forward_table);
			if(Px > 0){
				log_Pdataset += log2(Px) / data.size();
			}
		}
		_after_compute_log_Pdataset();
		return pow(2.0, -log_Pdataset / (double)dataset.size());
	}
	void Trainer::update_hyperparameters(){
		_model->update_hyperparameters();
	}
	void Trainer::show_assigned_words_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag, bool show_probability){
		_model->show_assigned_words_for_each_tag(dict, number_to_show_for_each_tag, show_probability);
	}
}