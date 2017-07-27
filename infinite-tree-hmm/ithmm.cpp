#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <string>
#include <set>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "src/hpylm.hpp"
#include "src/ithmm.h"
#include "src/util.h"

#define ID_BOS 0
#define ID_EOS 1
#define ID_UNK 2

class Model{
public:
	unordered_map<id, wstring> _dictionary;
	unordered_map<wstring, id> _dictionary_inv;
	id _autoincrement;
	double** _forward_table;		// 前向き確率計算用
	double** _decode_table;		// viterbiデコーディング用
	iTHMM* _ithmm;
	Model(){
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);

		_ithmm = new iTHMM();
		_forward_table = NULL;
		_decode_table = NULL;
	}
	~Model(){
		delete _ithmm;
	}
};

class Trainer: Model{
public:
	unordered_map<id, int> _word_count;
	vector<vector<Word*>> _dataset_train;
	vector<vector<Word*>> _dataset_test;
	vector<int> _rand_indices;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
	Trainer(){
		_ithmm = new iTHMM();
		_dictionary[ID_BOS] = L"<bos>";
		_dictionary[ID_EOS] = L"<eos>";
		_dictionary[ID_UNK] = L"<unk>";
		_autoincrement = ID_UNK + 1;

		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;

		_forward_table = NULL;
		_decode_table = NULL;
	}
	~Trainer(){
		delete _ithmm;
		for(int n = 0;n < _dataset_train.size();n++){
			vector<Word*> &data = _dataset_train[n];
			for(int m = 0;m < data.size();m++){
				Word* word = data[m];
				delete word;
			}
		}
		for(int n = 0;n < _dataset_test.size();n++){
			vector<Word*> &data = _dataset_test[n];
			for(int m = 0;m < data.size();m++){
				Word* word = data[m];
				delete word;
			}
		}
	}
	double get_alpha(){
		return _ithmm->_alpha;
	}
	double get_gamma(){
		return _ithmm->_gamma;
	}
	double get_lambda_alpha(){
		return _ithmm->_lambda_alpha;
	}
	double get_lambda_gamma(){
		return _ithmm->_lambda_gamma;
	}
	double get_strength(){
		return _ithmm->_strength;
	}
	double get_tau0(){
		return _ithmm->_tau0;
	}
	double get_tau1(){
		return _ithmm->_tau1;
	}
	double get_metropolis_hastings_acceptance_rate(){
		return _ithmm->_num_mh_acceptance / (double)(_ithmm->_num_mh_acceptance + _ithmm->_num_mh_rejection);
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
	python::list get_all_tags(){
		python::list tags;
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			string indices = "[" + node->_dump_indices() + "]";
			tags.append(indices);
		}
		return tags;
	}
	void set_alpha(double alpha){
		_ithmm->_alpha = alpha;
	}
	void set_gamma(double gamma){
		_ithmm->_gamma = gamma;
	}
	void set_lambda_alpha(double lambda_alpha){
		_ithmm->_lambda_alpha = lambda_alpha;
	}
	void set_lambda_gamma(double lambda_gamma){
		_ithmm->_lambda_gamma = lambda_gamma;
	}
	void set_strength(double strength){
		_ithmm->_strength = strength;
	}
	void set_tau0(double tau0){
		_ithmm->_tau0 = tau0;
	}
	void set_tau1(double tau1){
		_ithmm->_tau1 = tau1;
	}
	void set_depth_limit(int limit){
		_ithmm->set_depth_limit(limit);
	}
	void set_metropolis_hastings_enabled(bool enabled){
		_ithmm->_mh_enabled = enabled;
	}
	id add_string(wstring word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			_dictionary[_autoincrement] = word;
			_dictionary_inv[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	id string_to_word_id(wstring word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			return ID_UNK;
		}
		return itr->second;
	}
	void load_textfile(string filename, int train_split){
		c_printf("[*]%s\n", (boost::format("%sを読み込んでいます ...") % filename.c_str()).str().c_str());
		wifstream ifs(filename.c_str());
		wstring line_str;
		if (ifs.fail()){
			c_printf("[R]%s [*]%s", "エラー", (boost::format("%sを開けません.") % filename.c_str()).str().c_str());
			exit(1);
		}
		vector<wstring> lines;
		while (getline(ifs, line_str) && !line_str.empty()){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			lines.push_back(line_str);
		}
		assert(lines.size() > train_split);
		vector<int> rand_indices;
		for(int i = 0;i < lines.size();i++){
			rand_indices.push_back(i);
		}
		shuffle(rand_indices.begin(), rand_indices.end(), Sampler::mt);	// データをシャッフル
		for(int i = 0;i < rand_indices.size();i++){
			wstring &line_str = lines[rand_indices[i]];
			if(i < train_split){
				add_train_data(line_str);
			}else{
				add_test_data(line_str);
			}
		}
		cout << "train: " << _dataset_train.size() << endl;
		cout << "test:  " << _dataset_test.size() << endl;
		c_printf("[*]%s\n", (boost::format("%sを読み込みました.") % filename.c_str()).str().c_str());
	}
	void add_train_data(wstring line_str){
		_add_data_to(line_str, _dataset_train);
	}
	void add_test_data(wstring line_str){
		_add_data_to(line_str, _dataset_test);
	}
	void _add_data_to(wstring &line_str, vector<vector<Word*>> &dataset){
		vector<wstring> word_strs;
		split_word_by(line_str, L' ', word_strs);	// スペースで分割
		if(word_strs.size() > 0){
			vector<Word*> words;

			for(auto word_str: word_strs){
				if(word_str.size() == 0){
					continue;
				}
				Word* word = new Word();
				word->_id = add_string(word_str);
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
		for(int data_index = 0;data_index < _dataset_train.size();data_index++){
			vector<Word*> &data = _dataset_train[data_index];
			for(auto word = data.begin(), end = data.end();word != end;word++){
				id word_id = (*word)->_id;
				int count = get_count_for_word(word_id);
				if(count <= threshold){
					(*word)->_id = ID_UNK;
				}
			}
		}
	}
	void compile(){
		_ithmm->set_word_g0(1.0 / _word_count.size());
		_ithmm->initialize_data(_dataset_train);
	}
	void remove_all_data(){
		_ithmm->remove_all_data(_dataset_train);
		_ithmm->delete_invalid_children();
	}
	bool load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/ithmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _dictionary;
			iarchive >> _dictionary_inv;
			iarchive >> _autoincrement;
			ifs.close();
		}
		return _ithmm->load(dirname);
	}
	bool save(string dirname){
		// 辞書を保存
		std::ofstream ofs(dirname + "/ithmm.dict");
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		oarchive << _dictionary_inv;
		oarchive << _autoincrement;
		ofs.close();
		return _ithmm->save(dirname);
	}
	void perform_gibbs_sampling(){
		if(_rand_indices.size() != _dataset_train.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < _dataset_train.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		_ithmm->_num_mh_acceptance = 0;
		_ithmm->_num_mh_rejection = 0;
		shuffle(_rand_indices.begin(), _rand_indices.end(), Sampler::mt);	// データをシャッフル
		for(int n = 0;n < _dataset_train.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			vector<Word*> &data = _dataset_train[data_index];
			_ithmm->perform_gibbs_sampling_data(data);
		}
		_ithmm->delete_invalid_children();
	}
	void _before_viterbi_decode(vector<Node*> &nodes){
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
	python::list viterbi_decode_test(){
		python::list result_list;
		vector<Node*> nodes;
		_before_viterbi_decode(nodes);
		vector<Node*> series;
		for(int data_index = 0;data_index < _dataset_test.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return result_list;
			}
			vector<Word*> &data = _dataset_test[data_index];
			python::list tuple_list;
			viterbi_decode_data(data, nodes, series);
			for(int i = 0;i < data.size();i++){
				python::list tuple;
				wstring word = _dictionary[data[i]->_id];
				wstring tag = L"[" + series[i]->_wdump_indices() + L"]";
				tuple.append(word);
				tuple.append(tag);
				tuple_list.append(tuple);
			}
			result_list.append(tuple_list);
		}
		_after_viterbi_decode();
		return result_list;
	}
	python::list viterbi_decode_train(){
		python::list result_list;
		vector<Node*> nodes;
		_before_viterbi_decode(nodes);
		vector<Node*> series;
		for(int data_index = 0;data_index < _dataset_train.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return result_list;
			}
			vector<Word*> &data = _dataset_train[data_index];
			python::list tuple_list;
			viterbi_decode_data(data, nodes, series);
			for(int i = 0;i < data.size();i++){
				python::list tuple;
				wstring word = _dictionary[data[i]->_id];
				wstring tag = L"[" + series[i]->_wdump_indices() + L"]";
				tuple.append(word);
				tuple.append(tag);
				tuple_list.append(tuple);
			}
			result_list.append(tuple_list);
		}
		_after_viterbi_decode();
		return result_list;
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void viterbi_decode_data(vector<Word*> &data, vector<Node*> &states, vector<Node*> &series){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < states.size();i++){
			Node* state = states[i];
			Node* state_on_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_on_bos != NULL);
			double Ps = state_on_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			_forward_table[0][i] = Pword_given_s * Ps;
			_decode_table[0][i] = 0;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int j = 0;j < states.size();j++){
				Node* state = states[j];
				_forward_table[t][j] = 0;
				double max_value = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < states.size();i++){
					Node* prev_state = states[i];
					Node* state_on_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_on_prev_htssb != NULL);
					double Ps_given_prev = state_on_prev_htssb->_probability;
					double value = Ps_given_prev * _forward_table[t - 1][i];
					if(value > max_value){
						max_value = value;
						_forward_table[t][j] = value * Pword_given_s;
						_decode_table[t][j] = i;
					}
				}
			}
		}
		// 後ろ向きに系列を復元
		vector<int> series_indices;
		int n = data.size() - 1;
		int k = 0;
		double max_value = 0;
		for(int i = 0;i < states.size();i++){
			if(_forward_table[n][i] > max_value){
				k = i;
				max_value = _forward_table[n][i];
			}
		}
		series_indices.push_back(k);
		for(int t = n - 1;t >= 0;t--){
			k = _decode_table[t + 1][series_indices[n - t - 1]];
			series_indices.push_back(k);
		}
		std::reverse(series_indices.begin(), series_indices.end());
		// ノードをセット
		series.clear();
		for(int t = 0;t <= n;t++){
			int k = series_indices[t];
			series.push_back(states[k]);
		}
	}
	// データの対数尤度を計算
	// 前向きアルゴリズム
	double compute_Pdata(vector<Word*> &data, vector<Node*> &states){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < states.size();i++){
			Node* state = states[i];
			Node* state_on_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_on_bos != NULL);
			double Ps = state_on_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			_forward_table[0][i] = Pword_given_s * Ps;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int j = 0;j < states.size();j++){
				Node* state = states[j];
				_forward_table[t][j] = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < states.size();i++){
					Node* prev_state = states[i];
					Node* state_on_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_on_prev_htssb != NULL);
					double Ps_given_prev = state_on_prev_htssb->_probability;
					_forward_table[t][j] += Pword_given_s * Ps_given_prev * _forward_table[t - 1][i];
				}
			}
		}
		int t = data.size() - 1;
		double Px = 0;
		for(int j = 0;j < states.size();j++){
			Px += _forward_table[t][j];
		}
		return Px;
	}
	void _before_compute_log_Pdataset(vector<Node*> &nodes){
		// あらかじめ全HTSSBの棒の長さを計算しておく
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(auto node: nodes){
			double Peos_given_s = node->compute_transition_probability_to_eos(_ithmm->_tau0, _ithmm->_tau1);
			double total_stick_length = 1.0 - Peos_given_s;	// <eos>以外に遷移する確率をTSSBで分配する
			_ithmm->update_stick_length_of_tssb(node->_transition_tssb, total_stick_length, true);
		}
		_ithmm->update_stick_length_of_tssb(_ithmm->_bos_tssb, 1.0, false);
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
		return _compute_log_Pdataset(_dataset_train);
	}
	double compute_log_Pdataset_test(){
		return _compute_log_Pdataset(_dataset_test);
	}
	double _compute_log_Pdataset(vector<vector<Word*>> &dataset){
		vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			vector<Word*> &data = dataset[data_index];
			double Px = compute_Pdata(data, nodes);
			if(Px > 0){
				log_Pdataset += log(Px);
			}
		}
		_after_compute_log_Pdataset();
		return log_Pdataset;
	}
	double compute_log2_Pdataset_train(){
		return _compute_log2_Pdataset(_dataset_train);
	}
	double compute_log2_Pdataset_test(){
		return _compute_log2_Pdataset(_dataset_test);
	}
	double _compute_log2_Pdataset(vector<vector<Word*>> &dataset){
		vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			vector<Word*> &data = dataset[data_index];
			double Px = compute_Pdata(data, nodes);
			if(Px > 0){
				log_Pdataset += log2(Px);
			}
		}
		_after_compute_log_Pdataset();
		return log_Pdataset;
	}
	double compute_perplexity_train(){
		return _compute_perplexity(_dataset_train);
	}
	double compute_perplexity_test(){
		return _compute_perplexity(_dataset_test);
	}
	double _compute_perplexity(vector<vector<Word*>> &dataset){
		vector<Node*> nodes;
		_before_compute_log_Pdataset(nodes);
		// データごとの対数尤度を足していく
		double log_Pdataset = 0;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return 0;
			}
			vector<Word*> &data = dataset[data_index];
			double Px = compute_Pdata(data, nodes);
			if(Px > 0){
				log_Pdataset += log2(Px) / data.size();
			}
		}
		_after_compute_log_Pdataset();
		return pow(2.0, -log_Pdataset / (double)dataset.size());
	}
	void update_hyperparameters(){
		_ithmm->sample_hpylm_hyperparameters();
	}
	void show_assigned_words_for_each_tag(int number_to_show_for_each_tag, bool show_probability = true){
		auto pair = std::make_pair(0, 0);
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			string indices = node->_dump_indices();
			// linuxでバグるのでstringとwstring両方作る
			string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m" << endl;
			wstring wtab = L"";
			for(int i = 0;i < node->_depth_v;i++){
				wtab += L"	";
			}
			wcout << wtab;
			for(const auto &elem: ranking){
				id word_id = elem.first;
				wstring &word = _dictionary[word_id];
				double p = elem.second;
				int count = node->_num_word_assignment[word_id];
				wcout << "\x1b[1m" << word << "\x1b[0m" << L" (" << count;
				if(show_probability){
					wcout << L";p=" << p;
				} 
				wcout << L") ";
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			wcout << endl;
			ranking.clear();
		}
	}
	void show_assigned_words_and_probability_for_each_tag(int number_to_show_for_each_tag){
		auto pair = std::make_pair(0, 0);
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			string indices = node->_dump_indices();
			cout << "\x1b[32;1m" << "[" << indices << "]" << "\x1b[0m" << endl;
			for(const auto &elem: ranking){
				id word_id = elem.first;
				wstring &word = _dictionary[word_id];
				double p = elem.second;
				int count = node->_num_word_assignment[word_id];
				wcout << "\x1b[1m" << word << "\x1b[0m" << L"	" << count << L"	" << p << endl;
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			wcout << endl;
			ranking.clear();
		}
	}
	void show_hpylm_for_each_tag(){
		auto pair = std::make_pair(0, 0);
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			string indices = node->_dump_indices();
			// linuxでバグるのでstringとwstring両方作る
			string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m" << endl;
			wstring wtab = L"";
			for(int i = 0;i < node->_depth_v;i++){
				wtab += L"	";
			}
			wcout << wtab;
			for(const auto &table: node->_hpylm->_arrangement){
				id word_id = table.first;
				wstring &word = _dictionary[word_id];
				int num_tables = table.second.size();
				int num_customers = std::accumulate(table.second.begin(), table.second.end(), 0);
				wcout << "\x1b[1m" << word << "\x1b[0m" << L" (#t=" << num_tables << ";#c=" << num_customers << L") ";
			}
			wcout << endl;
		}
	}
	void show_sticks(){
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			string indices = node->_dump_indices();
			c_printf("[*]%s\n", ((boost::format("[%s]") % indices.c_str())).str().c_str());
			_show_stick(node);
		}
	}
	void _show_stick(Node* node_on_structure){
		assert(node_on_structure != NULL);
		double p_eos = node_on_structure->compute_transition_probability_to_eos(_ithmm->_tau0, _ithmm->_tau1);
		_ithmm->update_stick_length_of_tssb(node_on_structure->_transition_tssb, 1.0 - p_eos, true);

		vector<Node*> nodes;
		node_on_structure->_transition_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			string indices = node->_dump_indices();
			string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m " << node->_probability << endl;
		}
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<Trainer>("ithmm")
	.def("string_to_word_id", &Trainer::string_to_word_id)
	.def("add_string", &Trainer::add_string)
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling)
	.def("compile", &Trainer::compile)
	.def("load", &Trainer::load)
	.def("save", &Trainer::save)
	.def("load_textfile", &Trainer::load_textfile)
	.def("add_train_data", &Trainer::add_train_data)
	.def("add_test_data", &Trainer::add_test_data)
	.def("mark_low_frequency_words_as_unknown", &Trainer::mark_low_frequency_words_as_unknown)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("get_num_words", &Trainer::get_num_words)
	.def("get_count_for_word", &Trainer::get_count_for_word)
	.def("get_alpha", &Trainer::get_alpha)
	.def("get_gamma", &Trainer::get_gamma)
	.def("get_lambda_alpha", &Trainer::get_lambda_alpha)
	.def("get_lambda_gamma", &Trainer::get_lambda_gamma)
	.def("get_strength", &Trainer::get_strength)
	.def("get_tau0", &Trainer::get_tau0)
	.def("get_tau1", &Trainer::get_tau1)
	.def("get_metropolis_hastings_acceptance_rate", &Trainer::get_metropolis_hastings_acceptance_rate)
	.def("get_all_tags", &Trainer::get_all_tags)
	.def("set_alpha", &Trainer::set_alpha)
	.def("set_gamma", &Trainer::set_gamma)
	.def("set_lambda_alpha", &Trainer::set_lambda_alpha)
	.def("set_lambda_gamma", &Trainer::set_lambda_gamma)
	.def("set_strength", &Trainer::set_strength)
	.def("set_tau0", &Trainer::set_tau0)
	.def("set_tau1", &Trainer::set_tau1)
	.def("set_depth_limit", &Trainer::set_depth_limit)
	.def("set_metropolis_hastings_enabled", &Trainer::set_metropolis_hastings_enabled)
	.def("viterbi_decode_train", &Trainer::viterbi_decode_train)
	.def("viterbi_decode_test", &Trainer::viterbi_decode_test)
	.def("show_assigned_words_for_each_tag", &Trainer::show_assigned_words_for_each_tag)
	.def("show_assigned_words_and_probability_for_each_tag", &Trainer::show_assigned_words_and_probability_for_each_tag)
	.def("show_hpylm_for_each_tag", &Trainer::show_hpylm_for_each_tag)
	.def("show_sticks", &Trainer::show_sticks)
	.def("compute_perplexity_test", &Trainer::compute_perplexity_test)
	.def("compute_perplexity_train", &Trainer::compute_perplexity_train)
	.def("compute_log2_Pdataset_test", &Trainer::compute_log2_Pdataset_test)
	.def("compute_log2_Pdataset_train", &Trainer::compute_log2_Pdataset_train)
	.def("compute_log_Pdataset_test", &Trainer::compute_log_Pdataset_test)
	.def("compute_log_Pdataset_train", &Trainer::compute_log_Pdataset_train)
	.def("remove_all_data", &Trainer::remove_all_data);
}