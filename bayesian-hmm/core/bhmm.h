#ifndef _bhmm_
#define _bhmm_
#include <boost/format.hpp>

class BayesianHMM{
public:
	int _num_phrases;	// 品詞数
	int _num_words;		// 単語数
	int*** _trigram_table;	// 3-gramテーブル
	int** _marginal_counts_unigram;	// sum _trigram_table[t_-2][t_-1][.]
	int* _marginal_counts_bigram_unigram;	// sum _trigram_table[t_-2][.][.]
	HMM(int num_phrases, int num_words){
		_num_phrases = num_phrases;
		_num_words = num_words;
		init_table();
	}
	void init_table(){
		_trigram_table = (int***)malloc(_num_phrases * sizeof(int**));
		for(int tri_phrase = 0;tri_phrase < _num_phrases;tri_phrase++){
			_trigram_table[tri_phrase] = (int**)malloc(_num_phrases * sizeof(int*));
			for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
				_trigram_table[tri_phrase][bi_phrase] = (int*)calloc(_num_phrases, sizeof(int))
			}
		}
		_marginal_counts_unigram = (int**)malloc(_num_phrases * sizeof(int*));
		for(int bi_phrase = 0;bi_phrase < _num_phrases;bi_phrase++){
			_marginal_counts_unigram = (int*)calloc(_num_phrases, sizeof(int));
		}
	}
	int get_marginal_count_unigram(int tri_phrase, int bi_phrase){
		int sum_counts = 0;
		int** table = _trigram_table[tri_phrase][bi_phrase];
		for(int uni_phrase = 0;uni_phrase < _num_phrases;uni_phrase++){
			sum_counts += table[uni_phrase];
		}
		return sum_counts;
	}
	int get_marginal_count_bigram_unigram(int tri_phrase){
		int sum_counts = 0;
		int** table = _trigram_table[tri_phrase][bi_phrase];
		for(int uni_phrase = 0;uni_phrase < _num_phrases;uni_phrase++){
			sum_counts += table[uni_phrase];
		}
		return sum_counts;
	}
};

#endif