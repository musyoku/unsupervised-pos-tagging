#ifndef _bigram_lattice_
#define _bigram_lattice_
#include <boost/format.hpp>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <array>
#include <cfloat>
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
	Lattice(int max_num_words_in_sentence, int num_tags, HPYLM* pos_hpylm, HPYLM** word_hpylm_for_tag){
		_num_tags = num_tags;
		_max_num_words_in_sentence = max_num_words_in_sentence;
		_word_hpylm_for_tag = word_hpylm_for_tag;
		_pos_hpylm = pos_hpylm;
		_pos_context.push_back(0);
		_pos_context.push_back(0);
		_word_context.push_back(0);
		_word_context.push_back(0);
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

		for(int r = 0;r < _num_tags;r++){
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
				for(int q = 0;q < num_tags;q++){
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
		int token_t_id = sentence[t]->word_id;
		// <bos>2つの場合
		if(t == 2){
			if(q != BEGIN_OF_POS){
				_alpha[t][r][q] = 0;
				return;
			}
			_pos_context[0] = BEGIN_OF_POS;
			_pos_context[1] = BEGIN_OF_POS;
			double Pz_qr = _pos_hpylm->compute_Pw_h(r, _pos_context);
			HPYLM* word_hpylm = _word_hpylm_for_tag[r];
			_word_context[0] = BEGIN_OF_SENTENSE;
			_word_context[1] = BEGIN_OF_SENTENSE;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			_alpha[t][r][BEGIN_OF_POS] = Pt_h * Pz_qr;
			// cout << (boost::format("_alpha[%d][%d][%d] <- %f * %f") % t % r % BEGIN_OF_POS % Pt_h % Pz_qr).str() << endl;
			return;
		}
		// <bos>と何らかの単語
		if(t == 3){
			_pos_context[0] = BEGIN_OF_POS;
			_pos_context[1] = q;
			double Pz_qr = _pos_hpylm->compute_Pw_h(r, _pos_context);
			HPYLM* word_hpylm = _word_hpylm_for_tag[r];
			_word_context[0] = BEGIN_OF_SENTENSE;
			_word_context[1] = sentence[t - 1]->word_id;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			_alpha[t][r][q] = Pt_h * Pz_qr * _alpha[t - 1][q][BEGIN_OF_POS];
			// cout << (boost::format("_alpha[%d][%d][%d] <- %f * %f * %f") % t % r % q % Pt_h % Pz_qr % _alpha[t - 1][BEGIN_OF_POS][BEGIN_OF_POS]).str() << endl;
			return;
		}
		double sum = 0.1;
		_pos_context[0] = r;
		_pos_context[1] = q;
		for(int z = 0;z < _num_tags;z++){
			HPYLM* word_hpylm = _word_hpylm_for_tag[z];
			double Pz_qr = _pos_hpylm->compute_Pw_h(z, _pos_context);
			double _Pz_qr = _pos_hpylm->_compute_Pw_h(z, _pos_context);
			assert(Pz_qr == _Pz_qr);
			_word_context[0] = sentence[t - 2]->word_id;
			_word_context[1] = sentence[t - 1]->word_id;
			double Pt_h = word_hpylm->compute_Pw_h(token_t_id, _word_context);
			double _Pt_h = word_hpylm->_compute_Pw_h(token_t_id, _word_context);
			assert(Pt_h == _Pt_h);
			sum += Pt_h * Pz_qr * _alpha[t - 1][q][z];
			// cout << (boost::format("sum += %f * %f * %f") % Pt_h % Pz_qr % _alpha[t - 1][q][z]).str() << endl;
		}
		// cout << (boost::format("_alpha[%d][%d][%d] <- %f") % t % r % q % sum).str() << endl;
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
		while(t >= 2){
			if(argmax){
				argmax_backward_r_and_q(sentence, t, r, q);
			}else{
				sample_backward_r_and_q(sentence, t, r, q);
			}
			sentence[t]->tag_id = r;
			t--;
			if(t < 2){	// <bos>を書き換えてはいけない
				assert(q == BEGIN_OF_POS);
				break;
			}
			sentence[t]->tag_id = q;
			t--;
		}
	}
	// <eos>, EOPに接続する確率をもとにrとqをサンプリング
	void sample_starting_r_and_q(vector<Word*> &sentence, int &sampled_r, int &sampled_q){
		double sum_p = 0;
		int t = sentence.size() - 2;	// <eos>の1つ前
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
				assert(p >= 0);
				_sampling_table[r][q] = p;
				sum_p += p;
			}
		}
		double normalizer = 1.0 / sum_p;
		double bernoulli = Sampler::uniform(0, 1);
		sum_p = 0;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				sum_p += _sampling_table[r][q] * normalizer;
				if(bernoulli < sum_p){
					sampled_r = r;
					sampled_q = q;
					return;
				}
			}
		}
		sampled_r = _num_tags - 1;
		sampled_q = _num_tags - 1;
	}
	void sample_backward_r_and_q(vector<Word*> &sentence, int t, int &sampled_r, int &sampled_q){
		double sum_p = 0;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				sum_p += _alpha[t][r][q];
			}
		}
		double normalizer = 1.0 / sum_p;
		double bernoulli = Sampler::uniform(0, 1);
		sum_p = 0;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				sum_p += _alpha[t][r][q] * normalizer;
				if(bernoulli < sum_p){
					sampled_r = r;
					sampled_q = q;
					return;
				}
			}
		}
		sampled_r = _num_tags - 1;
		sampled_q = (t == 3) ? END_OF_POS : _num_tags - 1;
	}
	void argmax_backward_r_and_q(vector<Word*> &sentence, int t, int &sampled_r, int &sampled_q){
		double max_p = 0;
		int max_r = -1, max_q = -1;
		for(int r = 0;r < _num_tags;r++){
			for(int q = 0;q < _num_tags;q++){
				double p = _alpha[t][r][q];
				if(p > max_p){
					max_p = p;
					max_r = r;
					max_q = q;
				}
			}
		}
		sampled_r = max_r;
		sampled_q = max_q;
		assert(max_r != -1);
		assert(max_q != -1);
	}
	void perform_blocked_gibbs_sampling(vector<Word*> &sentence, bool argmax = false){
		this->forward_filtering(sentence);
		this->backward_sampling(sentence, argmax);
	}
};


#endif