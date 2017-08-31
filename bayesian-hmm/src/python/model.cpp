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
	void Model::set_initial_alpha(double alpha){
		_hmm->set_alpha(alpha);
	}
	void Model::set_initial_beta(double beta){
		_hmm->set_beta(beta);
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
	// 文の確率
	// 前向きアルゴリズム
	double Model::compute_p_sentence(std::vector<Word*> &sentence, double*** forward_table){
		assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		int tag_bos = 0;	// <s>
		for(int ti = 1;ti <= _hmm->_num_tags;ti++){
			int ti_2 = tag_bos;	// <s>
			int ti_1 = tag_bos;	// <s>
			id wi = sentence[2]->_id;
			double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
			double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
			assert(p_s_given_prev > 0);
			assert(p_w_given_s > 0);
			forward_table[2][tag_bos][ti] = p_w_given_s * p_s_given_prev;
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				forward_table[2][ti_1][ti] = 0;
			}
		}
		for(int i = 3;i < sentence.size() - 2;i++){
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				for(int ti = 1;ti <= _hmm->_num_tags;ti++){

					id wi = sentence[i]->_id;
					double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
					assert(p_w_given_s > 0);
					forward_table[i][ti_1][ti] = 0;
					if(i == 3){
						forward_table[i][ti_1][ti] += forward_table[i - 1][tag_bos][ti_1] * _hmm->compute_p_ti_given_t(ti, ti_1, tag_bos);
					}else{
						for(int ti_2 = 1;ti_2 <= _hmm->_num_tags;ti_2++){
							forward_table[i][ti_1][ti] += forward_table[i - 1][ti_2][ti_1] * _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
						}
					}
					forward_table[i][ti_1][ti] *= p_w_given_s;


				}
			}
		}
		int i = sentence.size() - 3;
		double p_x = 0;
		for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
			for(int ti = 1;ti <= _hmm->_num_tags;ti++){
				p_x += forward_table[i][ti_1][ti];
			}
		}
		return p_x;
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double*** forward_table, double*** decode_table){
		// assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		// for(int tag = 1;tag <= _hmm->_num_tags;tag++){
		// 	int ti_2 = sentence[0]->_state;
		// 	int ti_1 = sentence[1]->_state;
		// 	int ti = sentence[2]->_state;
		// 	id wi = sentence[2]->_id;
		// 	double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
		// 	double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
		// 	assert(p_s_given_prev > 0);
		// 	assert(p_w_given_s > 0);
		// 	forward_table[2][tag] = p_w_given_s * p_s_given_prev;
		// }
		// for(int i = 3;i < sentence.size() - 2;i++){
		// 	int ti_2 = sentence[i - 2]->_state;
		// 	id wi = sentence[i]->_id;
		// 	for(int target_tag = 1;target_tag <= _hmm->_num_tags;target_tag++){
		// 		forward_table[i][target_tag] = 0;
		// 		double max_value = 0;
		// 		double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, target_tag);
		// 		assert(p_w_given_s > 0);
		// 		for(int tag = 1;tag <= _hmm->_num_tags;tag++){
		// 			double p_s_given_prev = _hmm->compute_p_ti_given_t(target_tag, tag, ti_2);
		// 			assert(p_s_given_prev > 0);
		// 			double value = p_s_given_prev * forward_table[i - 1][tag];
		// 			if(value > max_value){
		// 				max_value = value;
		// 				forward_table[i][target_tag] = value * p_w_given_s;
		// 				decode_table[i][target_tag] = tag;
		// 			}

		// 		}
		// 	}
		// }
		// // 後ろ向きに系列を復元
		// sampled_state_sequence.clear();
		// int n = sentence.size() - 3;
		// int i = n;
		// int k = 0;
		// double max_value = 0;
		// for(int tag = 1;tag <= _hmm->_num_tags;tag++){
		// 	if(forward_table[i][tag] > max_value){
		// 		k = tag;
		// 		max_value = forward_table[i][tag];
		// 	}
		// }
		// sampled_state_sequence.push_back(k);
		// for(int i = n - 1;i >= 0;i--){
		// 	k = decode_table[i + 1][sampled_state_sequence[n - i - 1]];
		// 	sampled_state_sequence.push_back(k);
		// }
		// std::reverse(sampled_state_sequence.begin(), sampled_state_sequence.end());
		// assert(sampled_state_sequence.size() == sentence.size() - 4);
	}
}