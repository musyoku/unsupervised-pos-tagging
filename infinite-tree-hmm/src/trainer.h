#pragma once

class Trainer{
public:
	Dataset* _dataset;
	Model* _model;
	Dictionary* _dict;
	std::vector<int> _rand_indices;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
	double** _forward_table;		// 前向き確率計算用
	double** _decode_table;			// viterbiデコーディング用
	Trainer(Dataset* dataset, Model* model, Dictionary* dict){
		_dataset = dataset;
		_model = model;
		_dict = dict;
		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
		_forward_table = NULL;
		_decode_table = NULL;
		_model->_ithmm->set_word_g0(1.0 / _dataset->_word_count.size());
		_model->_ithmm->initialize_data(_dataset->_word_sequences_train);
	}
	void remove_all_data(){
		_model->_ithmm->remove_all_data(_dataset->_word_sequences_train);
		_model->_ithmm->delete_invalid_children();
	}
	void perform_gibbs_sampling(){
		if(_rand_indices.size() != _dataset->_word_sequences_train.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < _dataset->_word_sequences_train.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		_model->_ithmm->_num_mh_acceptance = 0;
		_model->_ithmm->_num_mh_rejection = 0;
		shuffle(_rand_indices.begin(), _rand_indices.end(), Sampler::mt);	// データをシャッフル
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
	void _before_viterbi_decode(std::vector<Node*> &nodes){
		_before_compute_log_Pdataset(nodes);
		_decode_table = new double*[_max_num_words_in_line];
		for(int i = 0;i < _max_num_words_in_line;i++){
			_decode_table[i] = new double[nodes.size()];
		}
	}
	void _after_viterbi_decode(){
		_after_compute_log_Pdataset();
		for(int i = 0;i < _max_num_words_in_line;i++){
			delete[] _decode_table[i];
		}
		delete[] _decode_table;
	}
	boost::python::list viterbi_decode_test(){
		boost::python::list result_list;
		std::vector<Node*> nodes;
		_before_viterbi_decode(nodes);
		std::vector<Node*> series;
		for(int data_index = 0;data_index < _dataset->_word_sequences_test.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return result_list;
			}
			std::vector<Word*> &data = _dataset->_word_sequences_test[data_index];
			boost::python::list tuple_list;
			_model->viterbi_decode_data(data, nodes, series, _forward_table, _decode_table);
			for(int i = 0;i < data.size();i++){
				boost::python::list tuple;
				std::wstring word = _dict->_id_to_str[data[i]->_id];
				std::wstring tag = L"[" + series[i]->_wdump_indices() + L"]";
				tuple.append(word);
				tuple.append(tag);
				tuple_list.append(tuple);
			}
			result_list.append(tuple_list);
		}
		_after_viterbi_decode();
		return result_list;
	}
	boost::python::list viterbi_decode_train(){
		boost::python::list result_list;
		std::vector<Node*> nodes;
		_before_viterbi_decode(nodes);
		std::vector<Node*> series;
		for(int data_index = 0;data_index < _dataset->_word_sequences_train.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return result_list;
			}
			std::vector<Word*> &data = _dataset->_word_sequences_train[data_index];
			boost::python::list tuple_list;
			_model->viterbi_decode_data(data, nodes, series, _forward_table, _decode_table);
			for(int i = 0;i < data.size();i++){
				boost::python::list tuple;
				std::wstring word = _dict->_id_to_str[data[i]->_id];
				std::wstring tag = L"[" + series[i]->_wdump_indices() + L"]";
				tuple.append(word);
				tuple.append(tag);
				tuple_list.append(tuple);
			}
			result_list.append(tuple_list);
		}
		_after_viterbi_decode();
		return result_list;
	}
	void _before_compute_log_Pdataset(std::vector<Node*> &nodes){
		// あらかじめ全HTSSBの棒の長さを計算しておく
		_model->_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(auto node: nodes){
			double Peos_given_s = node->compute_transition_probability_to_eos(_model->_ithmm->_tau0, _model->_ithmm->_tau1);
			double total_stick_length = 1.0 - Peos_given_s;	// <eos>以外に遷移する確率をTSSBで分配する
			_model->_ithmm->update_stick_length_of_tssb(node->_transition_tssb, total_stick_length, true);
		}
		_model->_ithmm->update_stick_length_of_tssb(_model->_ithmm->_bos_tssb, 1.0, false);
		// 計算用のテーブルを確保
		_forward_table = new double*[_max_num_words_in_line];
		for(int i = 0;i < _max_num_words_in_line;i++){
			_forward_table[i] = new double[nodes.size()];
		}
	}
	void _after_compute_log_Pdataset(){
		// 計算用のテーブルを解放
		for(int i = 0;i < _max_num_words_in_line;i++){
			delete[] _forward_table[i];
		}
		delete[] _forward_table;
	}
	// データセット全体の対数尤度を計算
	double compute_log_Pdataset_train(){
		return _compute_log_Pdataset(_dataset->_word_sequences_train);
	}
	double compute_log_Pdataset_test(){
		return _compute_log_Pdataset(_dataset->_word_sequences_test);
	}
	double _compute_log_Pdataset(std::vector<std::vector<Word*>> &dataset){
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
	double compute_log2_Pdataset_train(){
		return _compute_log2_Pdataset(_dataset->_word_sequences_train);
	}
	double compute_log2_Pdataset_test(){
		return _compute_log2_Pdataset(_dataset->_word_sequences_test);
	}
	double _compute_log2_Pdataset(std::vector<std::vector<Word*>> &dataset){
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
	double compute_perplexity_train(){
		return _compute_perplexity(_dataset->_word_sequences_train);
	}
	double compute_perplexity_test(){
		return _compute_perplexity(_dataset->_word_sequences_test);
	}
	double _compute_perplexity(std::vector<std::vector<Word*>> &dataset){
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
	void update_hyperparameters(){
		_model->_ithmm->sample_hpylm_hyperparameters();
	}
};