#ifndef _bigram_lattice_
#define _bigram_lattice_
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <array>
#include <cfloat>
#include <boost/format.hpp>
#include "hpylm.h"
#include "sampler.h"
#include "cprintf.h"

class Lattice{
public:
	HPYLM** _word_hpylm_for_tag;
	HPYLM* _pos_hpylm;
	double*** _alpha;
	double** _sampling_table;
	vector<int> _pos_context;
	vector<int> _word_context;
	int _max_num_words_in_sentence;
	int _num_tags;
	// max_sentence_lengthは1文に含まれる最大単語数
	Lattice(NPYLM3* npylm, Vocab* vocab, int max_num_words_in_sentence, int num_tags){
		_npylm = npylm;
		_vocab = vocab;
		_max_length = 10;
		_max_num_words_in_sentence = max_num_words_in_sentence;
		_num_tags = num_tags;
		_pos_context = {0, 0}
		_word_context = {0, 0}
		init_sampling_table(max_num_words_in_sentence, num_tags);
	}
	~Lattice(){
		int size = _max_num_words_in_sentence + 1;
		_alpha = new double**[size];
		for(int t = 0;t < size;t++){
			for(int r = 0;r < _num_tags;r++){
				delete[] _alpha[t][r];
			}
			delete[] _alpha[t];
		}
		delete[] _alpha;

		for(int r = 0;r < num_tags;r++){
			delete[] _sampling_table[r];
		}
		delete[] _sampling_table;
	}
	void init_sampling_table(int max_num_words_in_sentence, int num_tags){
		int size = max_num_words_in_sentence;
		_alpha = new double**[size];
		for(int t = 0;t < size;t++){
			_alpha[t] = new double*[num_tags];
			for(int r = 0;r < num_tags;r++){
				_alpha[t][r] = new double[num_tags];
				for(int r = 0;r < num_tags;r++){
					_alpha[t][r][q] = -1;
				}
			}
		}
		_sampling_table = new double*[num_tags];
		for(int r = 0;r < num_tags;r++){
			_sampling_table[r] = new double[num_tags];
		}
	}
	// alpha[t][r][q]の計算
	// word: j -> k -> t
	// pos:  z -> q -> r
	// 位置tの単語が品詞rから生成され、かつtより1つ前の単語が品詞qから生成される確率
	void compute_alpha_t_r_q(vector<Word*> &sentence, int t, int r, int q){
		assert(t >= 2);
		// <bos>2つの場合
		if(t == 2){
			if(q != BEGIN_OF_POS){
				return;
			}
			_pos_context[0] = BEGIN_OF_POS;
			_pos_context[1] = BEGIN_OF_POS;
			double Pz_qr = _pos_hpylm->compute_Pw_h(z, _pos_context);
			HPYLM* word_hpylm = _word_hpylm_for_tag[z];
			_word_context[0] = BEGIN_OF_SENTENSE;
			_word_context[1] = BEGIN_OF_SENTENSE;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			_alpha[t][r][BEGIN_OF_POS] = Pt_h * Pz_qr;
			return;
		}
		// <bos>と何らかの単語
		if(t == 3){
			if(q != BEGIN_OF_POS){
				return;
			}
			_pos_context[0] = q;
			_pos_context[1] = BEGIN_OF_POS;
			double Pz_qr = _pos_hpylm->compute_Pw_h(z, _pos_context);
			HPYLM* word_hpylm = _word_hpylm_for_tag[z];
			_word_context[0] = BEGIN_OF_SENTENSE;
			_word_context[1] = sentence[t - 1]->word_id;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			_alpha[t][r][BEGIN_OF_POS] = Pt_h * Pz_qr * _alpha[t - 1][BEGIN_OF_POS][BEGIN_OF_POS];
			return;
		}
		int token_t_id = sentence[t]->word_id;
		double sum = 0;
		_pos_context[0] = r;
		_pos_context[1] = q;
		for(int z = 0;z < _num_tags;z++){
			HPYLM* word_hpylm = _word_hpylm_for_tag[z];
			double Pz_qr = _pos_hpylm->compute_Pw_h(z, _pos_context);
			_word_context[0] = sentence[t - 2]->word_id;
			_word_context[1] = sentence[t - 1]->word_id;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			sum += Pt_h * Pz_qr * _alpha[t - 1][q][z];
		}
		_alpha[t][r][q] = sum;
	}
	void forward_filtering(vector<Word*> &sentence){
		// <bow>と<eos>の間の部分だけ考える
		for(int t = 2;t < sentence.size() - 1;t++){
			for(int r = 0;r < _num_tags;r++){
				for(int q = 0;q < _num_tags;q++){
					compute_alpha_t_r_q(sentence, t, r, q);
				}
			}
		}
	}
	void backward_sampling(vector<Word*> &sentence, bool argmax = false){
		int r = 0;
		int q = 0;
		// <eop>に繋がる確率からサンプリング
		sample_starting_r_and_q(sentence, r, q);
		int t = sentence.size() - 1;
		sentence[t]->tag_id = r;
		t--;
		sentence[t]->tag_id = q;
		t--;
		// 後ろからサンプリング
		for(;t >= 2;t -= 2){
			if(argmax){
				argmax_backward_k_and_j(sentence, t, k, j);
			}else{
				sample_backward_r_and_q(sentence, t, k, j);
			}
		}
	}
	// EOS, EOPに接続する確率をもとにrとqをサンプリング
	void sample_starting_r_and_q(vector<Word*> &sentence, int &sampled_r, int &sampled_q){
		double sum_p = 0;
		int t = sentence.size() - 1;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				_pos_context[0] = q;
				_pos_context[1] = r;
				double Pend_qr = _pos_hpylm->compute_Pw_h(END_OF_POS, _pos_context);
				HPYLM* word_hpylm = _word_hpylm_for_tag[r];
				_word_context[0] = sentence[t - 2]->word_id;
				_word_context[1] = sentence[t - 1]->word_id;
				double Pend_h = word_hpylm->compute_Pw_h(END_OF_SENTENSE, _word_context);
				double p = Pend_h * Pend_qr * _alpha[t][r][q];
				_sampling_table[r][q] = p;
				sum_p += p;
			}
		}
		double normalizer = 1.0 / sum_p;
		double r = Sampler::uniform(0, 1);
		sum_p = 0;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				sum_p += _sampling_table[r][q] * normalizer;
				if(r < sum_p){
					sampled_r = r;
					sampled_q = q;
					return;
				}
				i++;
			}
		}
		sampled_r = _num_tags - 1;
		sampled_q = _num_tags - 1;
	}
	void sample_backward_r_and_q(wstring &sentence, int t, int &sampled_k, int &sampled_j){
		vector<double> p_k;
		double sum_p = 0;
		int limit_k = min(t, _max_length);
		for(int k = 1;k <= limit_k;k++){
			for(int j = 1;j <= min(t - k, _max_length);j++){
				if(_alpha[t][k][j] == -1){
					c_printf("[r]%s [*]%s\n", "エラー:", (boost::format("前向き確率の計算に不備があります. t = %d, k = %d, j = %d") % t % k % j).str().c_str());
					exit(1);
				}
				double p = _alpha[t][k][j] + DBL_MIN;
				sum_p += p;
				p_k.push_back(p);
			}
			if(t - k == 0){
				double p = _alpha[t][k][0] + DBL_MIN;
				sum_p += p;
				p_k.push_back(p);
			}
		}
		double normalizer = 1.0 / sum_p;
		double r = Sampler::uniform(0, 1);
		int i = 0;
		sum_p = 0;
		for(int k = 1;k <= limit_k;k++){
			for(int j = 1;j <= min(t - k, _max_length);j++){
				sum_p += p_k[i] * normalizer;
				if(r < sum_p){
					sampled_k = k;
					sampled_j = j;
					return;
				}
				i++;
			}
			if(t - k == 0){
				sum_p += p_k[i] * normalizer;
				if(r < sum_p){
					sampled_k = k;
					sampled_j = 0;
					return;
				}
				i++;
			}
		}
	}
	void argmax_backward_k_and_j(wstring &sentence, int t, int &sampled_k, int &sampled_j){
		vector<double> p_k;
		double max_p = 0;
		int max_k, max_j;
		int limit_k = min(t, _max_length);
		for(int k = 1;k <= limit_k;k++){
			for(int j = 1;j <= min(t - k, _max_length);j++){
				if(_alpha[t][k][j] == -1){
					c_printf("[r]%s [*]%s\n", "エラー:", (boost::format("前向き確率の計算に不備があります. t = %d, k = %d, j = %d") % t % k % j).str().c_str());
					exit(1);
				}
				double p = _alpha[t][k][j] + DBL_MIN;
				if(p > max_p){
					max_p = p;
					max_k = k;
					max_j = j;
				}
			}

			if(t - k == 0){
				double p = _alpha[t][k][0] + DBL_MIN;
				if(p > max_p){
					max_p = p;
					max_k = k;
					max_j = 0;
				}
			}
		}
		sampled_k = max_k;
		sampled_j = max_j;
	}
	void perform_blocked_gibbs_sampling(wstring &sentence, vector<int> &segments, bool argmax = false){
		auto t1 = chrono::system_clock::now();
		int size = sentence.size() + 1;
		// for(int k = 0;k < size;k++){
		// 	for(int j = 0;j < size;j++){
		// 		for(int i = 0;i < size;i++){
		// 			_alpha[k][j][i] = -1;
		// 		}
		// 	}
		// }
		for(int i = 0;i < size;i++){
			// memset(_substring_token_id_cache[i], 0, size * sizeof(*_substring_token_id_cache[i]));
			for(int j = 0;j < size;j++){
				// cout << _substring_token_id_cache[i][j] << endl;
				_substring_token_id_cache[i][j] = 0;
			}
		}

		// auto t2 = chrono::system_clock::now();
		// auto t2_1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();

		// wcout << sentence << endl;
		this->forward_filtering(sentence);

		
		// auto t3 = chrono::system_clock::now();
		// auto t3_1 = chrono::duration_cast<chrono::milliseconds>(t3 - t2).count();

		this->backward_sampling(sentence, segments, argmax);

		// auto t4 = chrono::system_clock::now();
		// auto t4_1 = chrono::duration_cast<chrono::milliseconds>(t4 - t3).count();

		// if(sentence.size() == 790){
		// 	cout << t2_1 << ", " << t3_1 << ", " << t4_1 << endl;
		// 	exit(1);
		// }
	}
};


#endif