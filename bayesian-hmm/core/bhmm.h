#ifndef _bhmm_
#define _bhmm_
#include <boost/format.hpp>
#include <cassert>
#include <map>
#include <set>
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
	int* _word_types_for_tag;	// 同じ品詞の単語の種類
	map<pair<int, int>, int> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	double* _sampling_table;	// キャッシュ
	double _alpha;
	double _beta;
	BayesianHMM(){
		_trigram_counts = NULL;
		_bigram_counts = NULL;
		_unigram_counts = NULL;
		_sampling_table = NULL;
		_num_tags = -1;
		_num_words = -1;
		_alpha = 1;
		_beta = 1;
	}
	~BayesianHMM(){
		if(_trigram_counts != NULL){
			for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
				for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
					free(_trigram_counts[tri_tag][bi_tag]);
				}
				free(_trigram_counts[tri_tag]);
			}
			free(_trigram_counts);
		}
		if(_bigram_counts != NULL){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				free(_bigram_counts);
			}
			free(_bigram_counts);
		}
		if(_unigram_counts != NULL){
			free(_unigram_counts);
		}
		if(_sampling_table != NULL){
			free(_sampling_table);
		}
	}
	void set_num_tags(int number){
		_num_tags = number;
	}
	void set_num_words(int number){
		_num_words = number;
	}
	void set_alpha(double alpha){
		_alpha = alpha;
	}
	void set_beta(double beta){
		_beta = beta;
	}
	void init_ngram_counts_if_needed(vector<vector<Word*>> &dataset){
		if(_trigram_counts == NULL){
			init_ngram_counts(dataset);
		}
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		c_printf("[*]%s", "nグラムカウントを初期化してます ...");
		assert(_num_tags != -1);
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
		_word_types_for_tag = (int*)calloc(_num_tags, sizeof(int));
		// 最初は品詞をランダムに割り当てる
		set<pair<int, int>> tag_word_set;
		set<int> word_set;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			for(int pos = 2;pos < line.size();pos++){	// 3-gramなので3番目から.
				Word* word = line[pos];
				word->tag_id = Sampler::uniform_int(0, _num_tags - 1);
				increment_tag_word_count(word->tag_id, word->word_id);
				update_ngram_count(line[pos - 2], line[pos - 1], line[pos]);
				// 単語の種類数をカウント
				word_set.insert(word->word_id);
				// 同じタグの単語集合をカウント
				auto pair = std::make_pair(word->tag_id, word->word_id);
				auto itr = tag_word_set.find(pair);
				if(itr == tag_word_set.end()){
					tag_word_set.insert(pair);
					_word_types_for_tag[word->tag_id] += 1;
				}
			}
		}
		_num_words = word_set.size();
		c_printf("[*]%s", (boost::format("単語数: %d 行数: %d") % _num_words % dataset.size()).str().c_str());
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
		if(itr == _tag_word_counts.end() || itr->second < 1){
			c_printf("[R]%s [*]%s", "エラー", "品詞-単語ペアのカウントが正しく実装されていません.");
			exit(1);
		}
		_tag_word_counts[pair] -= 1;
	}
	int get_count_for_tag_word(int tag_id, int word_id){
		auto pair = std::make_pair(tag_id, word_id);
		auto itr = _tag_word_counts.find(pair);
		if(itr == _tag_word_counts.end()){
			return 0;
		}
		return itr->second;
	}
	// in:  t_{i-2},t_{i-1},t_i,t_{i+1},t_{i+2},w_i
	// out: void
	void add_tag_to_model_parameters(int _t_i_2, int _t_i_1, int t_i, int t_i_1, int t_i_2, int w_i){
		// 1-gram
		_unigram_counts[t_i] += 1;
		// 2-gram
		_bigram_counts[_t_i_1][t_i] += 1;
		_bigram_counts[t_i][t_i_1] += 1;
		// 3-gram
		_trigram_counts[_t_i_2][_t_i_1][t_i] += 1;
		_trigram_counts[_t_i_1][t_i][t_i_1] += 1;
		_trigram_counts[t_i][t_i_1][t_i_2] += 1;
		// 品詞-単語ペア
		auto pair = std::make_pair(t_i, w_i);
		auto itr = _tag_word_counts.find(pair);
		_tag_word_counts[pair] += 1;
	}
	// in:  t_{i-2},t_{i-1},t_i,t_{i+1},t_{i+2},w_i
	// out: void
	void remove_tag_from_model_parameters(int _t_i_2, int _t_i_1, int t_i, int t_i_1, int t_i_2, int w_i){
		// 1-gram
		_unigram_counts[t_i] -= 1;
		assert(_unigram_counts[t_i] >= 0);
		// 2-gram
		_bigram_counts[_t_i_1][t_i] -= 1;
		assert(_bigram_counts[_t_i_1][t_i] >= 0);
		_bigram_counts[t_i][t_i_1] -= 1;
		assert(_bigram_counts[t_i][t_i_1] >= 0);
		// 3-gram
		_trigram_counts[_t_i_2][_t_i_1][t_i] -= 1;
		assert(_trigram_counts[_t_i_2][_t_i_1][t_i] >= 0);
		_trigram_counts[_t_i_1][t_i][t_i_1] -= 1;
		assert(_trigram_counts[_t_i_1][t_i][t_i_1] >= 0);
		_trigram_counts[t_i][t_i_1][t_i_2] -= 1;
		assert(_trigram_counts[t_i][t_i_1][t_i_2] >= 0);
		// 品詞-単語ペア
		auto pair = std::make_pair(t_i, w_i);
		auto itr = _tag_word_counts.find(pair);
		assert(itr != _tag_word_counts.end());
		_tag_word_counts[pair] -= 1;
		assert(itr->second >= 0);
	}
	void perform_gibbs_sampling_with_line(vector<Word*> &line){
		if(_sampling_table == NULL){
			_sampling_table = (double*)malloc(_num_tags * sizeof(double));
		}
		for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
			int t_i = line[pos]->tag_id;
			int w_i = line[pos]->word_id;
			int t_i_1 = line[pos + 1]->tag_id;
			int t_i_2 = line[pos + 2]->tag_id;
			int _t_i_1 = line[pos - 1]->tag_id;
			int _t_i_2 = line[pos - 2]->tag_id;
			// t_iをモデルパラメータから除去
			remove_tag_from_model_parameters(_t_i_2, _t_i_1, t_i, t_i_1, t_i_2, w_i);
			// t_iを再サンプリング
			double sum = 0;
			int new_t_i = 0;
			for(int tag = 0;tag < _num_tags;tag++){
				_sampling_table[tag] = 1;
				double n_t_i_w_i = get_count_for_tag_word(tag, t_i);
				double n_t_i = _unigram_counts[t_i];
				double W_t_i = _word_types_for_tag[t_i];
				double n_t_2_1_i = _trigram_counts[_t_i_2][_t_i_1][t_i];
				double n_t_2_1 = _bigram_counts[_t_i_2][_t_i_1];
				double n_t_1_i_1 = _trigram_counts[_t_i_1][t_i][t_i_1];
				double n_t_1_i = _bigram_counts[_t_i_1][t_i];
				double I_2_1_i_1 = (_t_i_2 == _t_i_1 == t_i == t_i_1) ? 1 : 0;
				double I_2_1_i = (_t_i_2 == _t_i_1 == t_i) ? 1 : 0;
				double n_t_i_1_2 = _trigram_counts[t_i][t_i_1][t_i_2];
				double n_t_i_1 = _bigram_counts[t_i][t_i_1];
				double I_2_i_2_1_1 = (_t_i_2 == t_i == t_i_2 && _t_i_1 == t_i_1) ? 1 : 0;
				double I_1_i_1_2 = (_t_i_1 == t_i == t_i_1 == t_i_2) ? 1 : 0;
				double I_2_i_1_1 = (_t_i_2 == t_i && _t_i_1 == t_i_1) ? 1 : 0;
				double I_1_i_1 = (_t_i_1 == t_i == t_i_1) ? 1 : 0;
				_sampling_table[tag] *= (n_t_i_w_i + _beta) / (n_t_i + W_t_i * _beta);
				_sampling_table[tag] *= (n_t_2_1_i + _alpha) / (n_t_2_1 + _num_tags * _alpha);
				_sampling_table[tag] *= (n_t_1_i_1 + I_2_1_i_1 + _alpha) / (n_t_1_i + I_2_1_i + _num_tags * _alpha);
				_sampling_table[tag] *= (n_t_i_1_2 + I_2_i_2_1_1 + I_1_i_1_2 + _alpha) / (n_t_i_1 + I_2_i_1_1 + I_1_i_1 + _num_tags * _alpha);
				sum += _sampling_table[tag];
			}
			assert(sum > 0);
			double normalizer = 1.0 / sum;
			double r = Sampler::uniform(0, 1);
			for(int tag = 0;tag < _num_tags;tag++){
				sum -= _sampling_table[tag] * normalizer;
				if(sum < r){
					new_t_i = tag;
					break;
				}
			}
			// 新しいt_iをモデルパラメータに追加
			add_tag_to_model_parameters(_t_i_2, _t_i_1, new_t_i, t_i_1, t_i_2, w_i);
		}
	}
};

#endif