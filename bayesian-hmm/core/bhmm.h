#ifndef _bhmm_
#define _bhmm_
#include <boost/format.hpp>
#include <cassert>

class BayesianHMM{
public:
	int _num_phrases;	// 品詞数
	int _num_words;		// 単語数
	int*** _trigram_table;	// 3-gramテーブル
	int** _marginal_counts_tri_bi;	// sum _trigram_table[t_-2][t_-1][.]
	int* _marginal_counts_tri;	// sum _trigram_table[t_-2][.][.]
	HMM(){
		_num_phrases = -1;
		_num_words = -1;
	}
	~HMM(){
		for(int tri_phrase = 0;tri_phrase < _num_phrases;tri_phrase++){
			for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
				free(_trigram_table[tri_phrase][bi_phrase])
			}
			free(_trigram_table[tri_phrase])
		}
		free(_trigram_table);
		for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
			free(_marginal_counts_tri_bi)
		}
		free(_marginal_counts_tri_bi)
		free(_marginal_counts_tri);
	}
	void set_num_phrases(int number){
		_num_phrases = number;
	}
	void set_num_words(int number){
		_num_words = number;
	}
	void init_table(){
		cassert(_num_phrases != -1)
		cassert(_num_words != -1)
		_trigram_table = (int***)malloc(_num_phrases * sizeof(int**));
		for(int tri_phrase = 0;tri_phrase < _num_phrases;tri_phrase++){
			_trigram_table[tri_phrase] = (int**)malloc(_num_phrases * sizeof(int*));
			for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
				_trigram_table[tri_phrase][bi_phrase] = (int*)calloc(_num_phrases, sizeof(int))
			}
		}
		_marginal_counts_tri_bi = (int**)malloc(_num_phrases * sizeof(int*));
		for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
			_marginal_counts_tri_bi = (int*)calloc(_num_phrases, sizeof(int));
		}
		_marginal_counts_tri = (int*)calloc(_num_phrases, sizeof(int));
	}
	int get_marginal_count_unigram(int tri_phrase, int bi_phrase){
		return _marginal_counts_tri_bi[tri_phrase][bi_phrase];
	}
	int get_marginal_count_bigram_unigram(int tri_phrase){
		return _marginal_counts_tri[tri_phrase];
	}
};

#endif