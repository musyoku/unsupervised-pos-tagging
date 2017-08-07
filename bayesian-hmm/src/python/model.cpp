#include <iostream>
#include "model.h"

namespace bhmm {
	Model::Model(){
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_hmm = new HMM();
	}
	Model::~Model(){
		delete _hmm;
	}
	bool Model::load(std::string filename){
		return _hmm->load(filename);
	}
	bool Model::save(std::string filename){
		return _hmm->save(filename);
	}
	void Model::set_alpha(double alpha){
		_hmm->_alpha = alpha;
	}
	int Model::get_num_tags(){
		return _hmm->_num_tags;
	}
	double Model::get_temperature(){
		return _hmm->_temperature;
	}
	void Model::set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void Model::set_minimum_temperature(double temperature){
		_hmm->_minimum_temperature = temperature;
	}
	void Model::anneal_temperature(double temperature){
		_hmm->anneal_temperature(temperature);
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void Model::_viterbi_decode(std::vector<Word*> &data, std::vector<int> &sampled_state_sequence, double** forward_table, double** decode_table){
		// 初期化
		Word* word = data[0];
		for(int tag = 0;tag < get_num_tags();tag++){
			Node* state = all_states[i];
			Node* state_in_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_in_bos != NULL);
			double Ps = state_in_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			forward_table[0][i] = Pword_given_s * Ps;
			decode_table[0][i] = 0;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int tag = 0;tag < get_num_tags();tag++){
				forward_table[t][tag] = 0;
				double max_value = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < all_states.size();i++){
					Node* prev_state = all_states[i];
					Node* state_in_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_in_prev_htssb != NULL);
					double Ps_given_prev = state_in_prev_htssb->_probability;
					double value = Ps_given_prev * forward_table[t - 1][i];
					if(value > max_value){
						max_value = value;
						forward_table[t][tag] = value * Pword_given_s;
						decode_table[t][tag] = i;
					}
				}
			}
		}
		// 後ろ向きに系列を復元
		std::vector<int> series_indices;
		int n = data.size() - 1;
		int k = 0;
		double max_value = 0;
		for(int i = 0;i < all_states.size();i++){
			if(forward_table[n][i] > max_value){
				k = i;
				max_value = forward_table[n][i];
			}
		}
		series_indices.push_back(k);
		for(int t = n - 1;t >= 0;t--){
			k = decode_table[t + 1][series_indices[n - t - 1]];
			series_indices.push_back(k);
		}
		std::reverse(series_indices.begin(), series_indices.end());
		// ノードをセット
		sampled_state_sequence.clear();
		for(int t = 0;t <= n;t++){
			int k = series_indices[t];
			sampled_state_sequence.push_back(all_states[k]);
		}
	}
}