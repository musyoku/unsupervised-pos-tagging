#ifndef _bhmm_
#define _bhmm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <cmath>
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
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _num_tags;
		archive & _num_words;
		archive & _alpha;
		archive & _beta;
		archive & _temperature;
		archive & _trigram_counts;
		archive & _bigram_counts;
		archive & _unigram_counts;
		archive & _tag_word_counts;
	}
public:
	int _num_tags;	// 品詞数
	int _num_words;		// 単語数
	int*** _trigram_counts;	// 品詞3-gramのカウント
	int** _bigram_counts;	// 品詞2-gramのカウント
	int* _unigram_counts;	// 品詞1-gramのカウント
	map<int, map<int, int>> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	double* _sampling_table;	// キャッシュ
	double _alpha;
	double _beta;
	double _temperature;
	BayesianHMM(){
		_trigram_counts = NULL;
		_bigram_counts = NULL;
		_unigram_counts = NULL;
		_sampling_table = NULL;
		_num_tags = -1;
		_num_words = -1;
		_alpha = 0.003;
		_beta = 1;
		_temperature = 2;
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
	void set_temperature(double temperature){
		_temperature = temperature;
	}
	void anneal_temperature(double multiplier){
		_temperature *= multiplier;
	}
	void init_ngram_counts_if_needed(vector<vector<Word*>> &dataset){
		if(_trigram_counts == NULL){
			init_ngram_counts(dataset);
		}
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		c_printf("[*]%s\n", "nグラムカウントを初期化してます ...");
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
		// 最初は品詞をランダムに割り当てる
		set<pair<int, int>> tag_word_set;
		set<int> word_set;
		map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			// pos < 2
			// <bos>2つ
			_bigram_counts[line[0]->tag_id][line[1]->tag_id] += 1;
			_unigram_counts[line[0]->tag_id] += 1;
			_unigram_counts[line[1]->tag_id] += 1;
			word_set.insert(line[0]->word_id);
			word_set.insert(line[1]->word_id);
			tag_word_set.insert(std::make_pair(line[0]->tag_id, line[0]->word_id));
			increment_tag_word_count(line[0]->tag_id, line[0]->word_id);
			tag_word_set.insert(std::make_pair(line[1]->tag_id, line[1]->word_id));
			increment_tag_word_count(line[1]->tag_id, line[1]->word_id);
			// line.size() - 2 > pos >= 2
			for(int pos = 2;pos < line.size() - 2;pos++){	// 3-gramなので3番目から.
				Word* word = line[pos];
				auto itr = tag_for_word.find(word->word_id);
				if(itr == tag_for_word.end()){
					word->tag_id = Sampler::uniform_int(0, _num_tags - 1);
					tag_for_word[word->word_id] = word->tag_id;
				}else{
					word->tag_id = itr->second;
				}
				update_ngram_count(line[pos - 2], line[pos - 1], line[pos]);
				// 単語の種類数をカウント
				word_set.insert(word->word_id);
				// 同じタグの単語集合をカウント
				increment_tag_word_count(word->tag_id, word->word_id);
			}
			// pos >= line.size() - 2
			// <eos>2つ
			int end_index = line.size() - 1;
			int t_end = line[end_index]->tag_id;
			int t_end_1 = line[end_index - 1]->tag_id;
			int t_end_2 = line[end_index - 2]->tag_id;
			int t_end_3 = line[end_index - 3]->tag_id;
			_unigram_counts[t_end] += 1;
			_unigram_counts[t_end_1] += 1;
			_bigram_counts[t_end_2][t_end_1] += 1;
			_bigram_counts[t_end_1][t_end] += 1;
			_trigram_counts[t_end_3][t_end_2][t_end_1] += 1;
			_trigram_counts[t_end_2][t_end_1][t_end] += 1;
			int w_end = line[end_index]->word_id;
			int w_end_1 = line[end_index - 1]->word_id;
			word_set.insert(w_end);
			word_set.insert(w_end_1);
			tag_word_set.insert(std::make_pair(t_end, w_end));
			increment_tag_word_count(t_end, w_end);
			tag_word_set.insert(std::make_pair(t_end_1, w_end_1));
			increment_tag_word_count(t_end_1, w_end_1);
		}
		_num_words = word_set.size();
		c_printf("[*]%s\n", (boost::format("単語数: %d - 行数: %d") % _num_words % dataset.size()).str().c_str());
	}
	void update_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word){
		_trigram_counts[tri_word->tag_id][bi_word->tag_id][uni_word->tag_id] += 1;
		_bigram_counts[bi_word->tag_id][uni_word->tag_id] += 1;
		_unigram_counts[uni_word->tag_id] += 1;
	}
	void increment_tag_word_count(int tag_id, int word_id){
		map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			word_counts[word_id] = 1;
			return;
		}
		itr->second += 1;
	}
	void decrement_tag_word_count(int tag_id, int word_id){
		map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			c_printf("[R]%s [*]%s", "エラー", "品詞-単語ペアのカウントが正しく実装されていません.");
			exit(1);
		}
		itr->second -= 1;
		if(itr->second <= 0){
			word_counts.erase(itr);
		}
	}
	int get_count_for_tag_word(int tag_id, int word_id){
		map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			return 0;
		}
		return itr->second;
	}
	int get_word_types_for_tag(int tag_id){
		map<int, int> &word_counts = _tag_word_counts[tag_id];
		return word_counts.size();
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
		increment_tag_word_count(t_i, w_i);
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
		decrement_tag_word_count(t_i, w_i);
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
				double n_t_i_w_i = get_count_for_tag_word(tag, w_i);
				double n_t_i = _unigram_counts[tag];
				double W_t_i = get_word_types_for_tag(tag);
				double n_t_2_1_i = _trigram_counts[_t_i_2][_t_i_1][tag];
				double n_t_2_1 = _bigram_counts[_t_i_2][_t_i_1];
				double n_t_1_i_1 = _trigram_counts[_t_i_1][tag][t_i_1];
				double n_t_1_i = _bigram_counts[_t_i_1][tag];
				double I_2_1_i_1 = (_t_i_2 == _t_i_1 == tag == t_i_1) ? 1 : 0;
				double I_2_1_i = (_t_i_2 == _t_i_1 == tag) ? 1 : 0;
				double n_t_i_1_2 = _trigram_counts[tag][t_i_1][t_i_2];
				double n_t_i_1 = _bigram_counts[tag][t_i_1];
				double I_2_i_2_1_1 = (_t_i_2 == tag == t_i_2 && _t_i_1 == t_i_1) ? 1 : 0;
				double I_1_i_1_2 = (_t_i_1 == tag == t_i_1 == t_i_2) ? 1 : 0;
				double I_2_i_1_1 = (_t_i_2 == tag && _t_i_1 == t_i_1) ? 1 : 0;
				double I_1_i_1 = (_t_i_1 == tag == t_i_1) ? 1 : 0;
				_sampling_table[tag] *= (n_t_i_w_i + _beta) / (n_t_i + W_t_i * _beta);
				_sampling_table[tag] *= (n_t_2_1_i + _alpha) / (n_t_2_1 + _num_tags * _alpha);
				_sampling_table[tag] *= (n_t_1_i_1 + I_2_1_i_1 + _alpha) / (n_t_1_i + I_2_1_i + _num_tags * _alpha);
				_sampling_table[tag] *= (n_t_i_1_2 + I_2_i_2_1_1 + I_1_i_1_2 + _alpha) / (n_t_i_1 + I_2_i_1_1 + I_1_i_1 + _num_tags * _alpha);
				_sampling_table[tag] = pow(_sampling_table[tag], 1.0 / _temperature);
				sum += _sampling_table[tag];
			}
			assert(sum > 0);
			double normalizer = 1.0 / sum;
			double r = Sampler::uniform(0, 1);
			sum = 0;
			for(int tag = 0;tag < _num_tags;tag++){
				sum += _sampling_table[tag] * normalizer;
				if(sum >= r){
					new_t_i = tag;
					break;
				}
			}
			// 新しいt_iをモデルパラメータに追加
			add_tag_to_model_parameters(_t_i_2, _t_i_1, new_t_i, t_i_1, t_i_2, w_i);
			line[pos]->tag_id = new_t_i;
		}
	}
	int get_most_co_occurring_tag(int word_id){
		int max_count = 0;
		int most_co_occurring_tag_id = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			int count = get_count_for_tag_word(tag, word_id);
			if(count > max_count){
				max_count = count;
				most_co_occurring_tag_id = tag;
			}
		}
		return most_co_occurring_tag_id;
	}
	void dump_trigram_counts(){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
					cout << (boost::format("3-gram [%d][%d][%d] = %d") % tri_tag % bi_tag % uni_tag % _trigram_counts[tri_tag][bi_tag][uni_tag]).str() << endl;
				}
			}
		}
	}
	void dump_bigram_counts(){
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
				cout << (boost::format("2-gram [%d][%d] = %d") % bi_tag % uni_tag % _bigram_counts[bi_tag][uni_tag]).str() << endl;
			}
		}
	}
	void dump_unigram_counts(){
		for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
			cout << (boost::format("1-gram [%d] = %d") % uni_tag % _unigram_counts[uni_tag]).str() << endl;
		}
	}
	void dump_word_types(){
		for(int tag = 0;tag < _num_tags;tag++){
			cout << tag << ": " << get_word_types_for_tag(tag) << endl;
		}
	}
	bool save(string filename = "hmm.model"){
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << static_cast<const BayesianHMM&>(*this);
		ofs.close();
		return true;
	}
	bool load(string filename = "hmm.model"){
		std::ifstream ifs(filename);
		if(ifs.good() == false){
			return false;
		}
		boost::archive::binary_iarchive iarchive(ifs);
		iarchive >> *this;
		ifs.close();
		return true;
	}
};

#endif