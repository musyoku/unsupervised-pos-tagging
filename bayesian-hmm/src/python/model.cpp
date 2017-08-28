#include <iostream>
#include "model.h"

namespace bhmm {
	Model::Model(int num_tags){
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_hmm = new HMM(num_tags);
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
	double Model::compute_p_sentence(std::vector<Word*> &sentence, double** forward_table){
		assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		for(int tag = 1;tag <= _hmm->_num_tags;tag++){
			int ti_2 = sentence[0]->_state;
			int ti_1 = sentence[1]->_state;
			int ti = sentence[2]->_state;
			id wi = sentence[2]->_id;
			double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
			double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
			assert(p_s_given_prev > 0);
			assert(p_w_given_s > 0);
			forward_table[2][tag] = p_w_given_s * p_s_given_prev;
		}
		for(int i = 3;i < sentence.size() - 2;i++){
			int ti_2 = sentence[i - 2]->_state;
			int ti_1 = sentence[i - 1]->_state;
			id wi = sentence[i]->_id;
			for(int target_tag = 1;target_tag <= _hmm->_num_tags;target_tag++){
				forward_table[i][target_tag] = 0;
				for(int tag = 1;tag <= _hmm->_num_tags;tag++){
					double p_s_given_prev = _hmm->compute_p_ti_given_t(tag, ti_1, ti_2);
					double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, tag);
					assert(p_s_given_prev > 0);
					assert(p_w_given_s > 0);
					forward_table[i][target_tag] += p_w_given_s * p_s_given_prev * forward_table[i - 1][tag];
				}
			}
		}
		int i = sentence.size() - 2;
		double p_x = 0;
		for(int tag = 1;tag <= _hmm->_num_tags;tag++){
			p_x += forward_table[i][tag];
		}
		return p_x;
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double** forward_table, double** decode_table){
		// 初期化
		// int num_tags = _hmm->_num_tags;
		// int seq_length = sentence.size();
		// assert(seq_length > 4);	// <s><s></s></s>
		// Word* word = sentence[0];
		// for(int t = 1;t < seq_length;t++){
		// 	Word* word = sentence[t];
		// 	for(int tag = 0;tag < get_num_tags();tag++){
		// 		forward_table[t][tag] = 0;
		// 		double max_value = 0;
		// 		double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
		// 		for(int i = 0;i < all_states.size();i++){
		// 			Node* prev_state = all_states[i];
		// 			Node* state_in_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
		// 			assert(state_in_prev_htssb != NULL);
		// 			double Ps_given_prev = state_in_prev_htssb->_probability;
		// 			double value = Ps_given_prev * forward_table[t - 1][i];
		// 			if(value > max_value){
		// 				max_value = value;
		// 				forward_table[t][tag] = value * Pword_given_s;
		// 				decode_table[t][tag] = i;
		// 			}
		// 		}
		// 	}
		// }
		// // 後ろ向きに系列を復元
		// std::vector<int> series_indices;
		// int n = sentence.size() - 1;
		// int k = 0;
		// double max_value = 0;
		// for(int i = 0;i < all_states.size();i++){
		// 	if(forward_table[n][i] > max_value){
		// 		k = i;
		// 		max_value = forward_table[n][i];
		// 	}
		// }
		// series_indices.push_back(k);
		// for(int t = n - 1;t >= 0;t--){
		// 	k = decode_table[t + 1][series_indices[n - t - 1]];
		// 	series_indices.push_back(k);
		// }
		// std::reverse(series_indices.begin(), series_indices.end());
		// // ノードをセット
		// sampled_state_sequence.clear();
		// for(int t = 0;t <= n;t++){
		// 	int k = series_indices[t];
		// 	sampled_state_sequence.push_back(all_states[k]);
		// }
	}
}