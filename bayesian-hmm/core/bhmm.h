#ifndef _bhmm_
#define _bhmm_
#include <boost/format.hpp>
#include <cassert>
#include <map>
#include "c_printf.h"
#include "sampler.h"
#include "util.h"
using namespace std;

typedef struct Word {
	int word_id;
	int tag_id;
} Word;

class BayesianHMM{
public:
	int _num_tags;	// 品詞数
	int _num_words;		// 単語数
	int*** _trigram_table;	// 品詞3-gramテーブル
	int** _marginal_counts_tri_bi;	// sum _trigram_table[t_-2][t_-1][.]
	int* _marginal_counts_tri;	// sum _trigram_table[t_-2][.][.]
	map<pair<int, int>, int> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	BayesianHMM(){
		_num_tags = -1;
		_num_words = -1;
	}
	~BayesianHMM(){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				free(_trigram_table[tri_tag][bi_tag]);
			}
			free(_trigram_table[tri_tag]);
		}
		free(_trigram_table);
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			free(_marginal_counts_tri_bi);
		}
		free(_marginal_counts_tri_bi);
		free(_marginal_counts_tri);
	}
	void set_num_tags(int number){
		_num_tags = number;
	}
	void set_num_words(int number){
		_num_words = number;
	}
	void init_table_if_needed(vector<vector<Word*>> &dataset){
		if(_trigram_table == NULL){
			init_table(dataset);
		}
	}
	void init_table(vector<vector<Word*>> &dataset){
		assert(_num_tags != -1);
		assert(_num_words != -1);
		// メモリ確保
		_trigram_table = (int***)malloc(_num_tags * sizeof(int**));
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			_trigram_table[tri_tag] = (int**)malloc(_num_tags * sizeof(int*));
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				_trigram_table[tri_tag][bi_tag] = (int*)calloc(_num_tags, sizeof(int));
			}
		}
		_marginal_counts_tri_bi = (int**)malloc(_num_tags * sizeof(int*));
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			_marginal_counts_tri_bi[bi_tag] = (int*)calloc(_num_tags, sizeof(int));
		}
		_marginal_counts_tri = (int*)calloc(_num_tags, sizeof(int));

		// 最初は品詞をランダムに割り当てる
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			for(int pos = 0;pos < line.size();pos++){
				Word* word = line[pos];
				word->tag_id = Sampler::uniform_int(0, _num_tags - 1);
				increment_tag_word_count(word->tag_id, word->word_id);
			}
		}
	}
	void increment_tag_word_count(int tag_id, int word_id){
		auto pair = std::make_pair(tag_id, word_id);
		auto itr = _tag_word_counts.find(pair);
		if(itr == _tag_word_counts.end()){
			_tag_word_counts[pair] = 1;
			return;
		}
		_tag_word_counts[pair] += 1;
	}
	int get_marginal_count_unigram(int tri_tag, int bi_tag){
		return _marginal_counts_tri_bi[tri_tag][bi_tag];
	}
	int get_marginal_count_bigram_unigram(int tri_tag){
		return _marginal_counts_tri[tri_tag];
	}
};

#endif