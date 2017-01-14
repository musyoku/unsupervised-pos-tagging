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
	int*** _trigram_counts;	// 品詞3-gramのカウント
	int** _bigram_counts;	// 品詞2-gramのカウント
	int* _unigram_counts;	// 品詞1-gramのカウント
	map<pair<int, int>, int> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	BayesianHMM(){
		_num_tags = -1;
		_num_words = -1;
	}
	~BayesianHMM(){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				free(_trigram_counts[tri_tag][bi_tag]);
			}
			free(_trigram_counts[tri_tag]);
		}
		free(_trigram_counts);
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			free(_bigram_counts);
		}
		free(_bigram_counts);
		free(_unigram_counts);
	}
	void set_num_tags(int number){
		_num_tags = number;
	}
	void set_num_words(int number){
		_num_words = number;
	}
	void init_ngram_counts_if_needed(vector<vector<Word*>> &dataset){
		if(_trigram_counts == NULL){
			init_ngram_counts(dataset);
		}
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		assert(_num_tags != -1);
		assert(_num_words != -1);
		// メモリ確保
		_trigram_counts = (int***)malloc(_num_tags * sizeof(int**));
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			_trigram_counts[tri_tag] = (int**)malloc(_num_tags * sizeof(int*));
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				_trigram_counts[tri_tag][bi_tag] = (int*)calloc(_num_tags, sizeof(int));
			}
		}
		_bigram_counts = (int**)malloc(_num_tags * sizeof(int*));
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			_bigram_counts[bi_tag] = (int*)calloc(_num_tags, sizeof(int));
		}
		_unigram_counts = (int*)calloc(_num_tags, sizeof(int));

		// 最初は品詞をランダムに割り当てる
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			for(int pos = 2;pos < line.size() - 1;pos++){	// 3-gramなので3番目から. <eos>は無視
				Word* word = line[pos];
				word->tag_id = Sampler::uniform_int(0, _num_tags - 1);
				increment_tag_word_count(word->tag_id, word->word_id);
				update_ngram_count(line[pos - 2], line[pos - 1], line[pos]);
			}
		}
	}
	void update_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word){
		_trigram_counts[tri_word->tag_id][bi_word->tag_id][uni_word->tag_id] += 1;
		_bigram_counts[bi_word->tag_id][uni_word->tag_id] += 1;
		_unigram_counts[uni_word->tag_id] += 1;
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
	void decrement_tag_word_count(int tag_id, int word_id){
		auto pair = std::make_pair(tag_id, word_id);
		auto itr = _tag_word_counts.find(pair);
		if(itr == _tag_word_counts.end()){
			c_printf("[R]%s [*]%s", "エラー", "品詞-単語ペアのカウントが正しく実装されていません.");
			exit(1);
		}
		_tag_word_counts[pair] -= 1;
	}
	int get_marginal_count_unigram(int tri_tag, int bi_tag){
		return _bigram_counts[tri_tag][bi_tag];
	}
	int get_marginal_count_bigram_unigram(int tri_tag){
		return _unigram_counts[tri_tag];
	}
};

#endif