#ifndef _bhmm_
#define _bhmm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <set>
#include "cprintf.h"
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
		archive & _temperature;
		archive & _minimum_temperature;
		archive & _tag_word_counts;
	}
public:
	int _num_tags;	// 品詞数
	int _num_words;		// 単語数
	int*** _trigram_counts;	// 品詞3-gramのカウント
	int** _bigram_counts;	// 品詞2-gramのカウント
	int* _unigram_counts;	// 品詞1-gramのカウント
	int* _Wt;
	unordered_map<int, unordered_map<int, int>> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	double* _sampling_table;	// キャッシュ
	double _alpha;
	double* _beta;
	double* _new_beta;
	double _temperature;
	double _minimum_temperature;
	BayesianHMM(){
		_trigram_counts = NULL;
		_bigram_counts = NULL;
		_unigram_counts = NULL;
		_sampling_table = NULL;
		_Wt = NULL;
		_num_tags = -1;
		_num_words = -1;
		_alpha = 0.003;
		_beta = NULL;
		_new_beta = NULL;
		_temperature = 1;
		_minimum_temperature = 1;
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
		if(_beta != NULL){
			free(_beta);
		}
		if(_new_beta != NULL){
			free(_new_beta);
		}
		if(_Wt != NULL){
			free(_Wt);
		}
	}
	void anneal_temperature(double multiplier){
		if(_temperature > _minimum_temperature){
			_temperature *= multiplier;
		}
	}
	void initialize(vector<vector<Word*>> &dataset){
		// その他
		alloc_table();
		// nグラムのカウントテーブル
		init_ngram_counts(dataset);
	}
	void alloc_table(){
		assert(_num_tags != -1);
		// 各タグの可能な単語数
		_Wt = (int*)calloc(_num_tags, sizeof(int));
		// Betaの初期化
		// 初期値は1
		_beta = (double*)malloc(_num_tags * sizeof(double));
		for(int tag = 0;tag < _num_tags;tag++){
			_beta[tag] = 1;
		}
		// サンプリングした新しいBetaの一時保存用
		_new_beta = (double*)calloc(_num_tags, sizeof(double));
		// 3-gram		
		_trigram_counts = (int***)malloc(_num_tags * sizeof(int**));
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			_trigram_counts[tri_tag] = (int**)malloc(_num_tags * sizeof(int*));
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				_trigram_counts[tri_tag][bi_tag] = (int*)calloc(_num_tags, sizeof(int));
			}
		}
		// 2-gram
		_bigram_counts = (int**)malloc(_num_tags * sizeof(int*));
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			_bigram_counts[bi_tag] = (int*)calloc(_num_tags, sizeof(int));
		}
		// 1-gram
		_unigram_counts = (int*)calloc(_num_tags, sizeof(int));
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		c_printf("[*]%s\n", "n-gramモデルを構築してます ...");
		assert(_num_tags != -1);
		// 最初は品詞をランダムに割り当てる
		set<int> word_set;
		unordered_map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			// pos < 2
			// <bos>2つ
			_bigram_counts[line[0]->tag_id][line[1]->tag_id] += 1;
			_unigram_counts[line[0]->tag_id] += 1;
			_unigram_counts[line[1]->tag_id] += 1;
			word_set.insert(line[0]->word_id);
			word_set.insert(line[1]->word_id);
			increment_tag_word_count(line[0]->tag_id, line[0]->word_id);
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
			increment_tag_word_count(t_end, w_end);
			increment_tag_word_count(t_end_1, w_end_1);
		}
		_num_words = word_set.size();
		c_printf("[*]%s\n", (boost::format("単語数: %d - 行数: %d") % _num_words % dataset.size()).str().c_str());
	}
	void set_Wt_for_tag(int tag_id, int number){
		assert(_Wt != NULL);
		assert(tag_id < _num_tags);
		_Wt[tag_id] = number;
	}
	void update_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word){
		_trigram_counts[tri_word->tag_id][bi_word->tag_id][uni_word->tag_id] += 1;
		_bigram_counts[bi_word->tag_id][uni_word->tag_id] += 1;
		_unigram_counts[uni_word->tag_id] += 1;
	}
	void increment_tag_word_count(int tag_id, int word_id){
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			word_counts[word_id] = 1;
			return;
		}
		itr->second += 1;
	}
	void decrement_tag_word_count(int tag_id, int word_id){
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
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
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			return 0;
		}
		return itr->second;
	}
	int get_word_types_for_tag(int tag_id){
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		return word_counts.size();
	}
	double compute_log_Pt_alpha(vector<Word*> &line){
		double log_Pt_alpha = 0;
		for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
			int ti_2 = line[pos - 2]->tag_id;
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
			double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
			double Pt_i_alpha = (n_ti_2_ti_1_ti + _alpha) / (n_ti_2_ti_1 + _num_tags * _alpha);
			log_Pt_alpha += log(Pt_i_alpha);
		}
		return log_Pt_alpha;
	}
	double compute_log_Pw_t_alpha(vector<Word*> &line){
		double log_Pt_alpha = 0;
		for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
			int ti_2 = line[pos - 2]->tag_id;
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
			double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
			double Pt_i_alpha = (n_ti_2_ti_1_ti + _alpha) / (n_ti_2_ti_1 + _num_tags * _alpha);
			log_Pt_alpha += log(Pt_i_alpha);
		}
		return log_Pt_alpha;
	}
	// 正規化定数で割る前の値
	double compute_Pti_wi_beta(int ti, int wi, double beta){
		double n_ti_wi = get_count_for_tag_word(ti, wi);
		double n_ti = _unigram_counts[ti];
		double W_ti = _Wt[ti];
		return (n_ti_wi + beta) / (n_ti + W_ti * beta);
	}
	// in:  t_{i-2},t_{i-1},ti,t_{i+1},t_{i+2},w_i
	// out: void
	void add_tag_to_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
		// 1-gram
		_unigram_counts[ti] += 1;
		// 2-gram
		_bigram_counts[ti_1][ti] += 1;
		_bigram_counts[ti][ti1] += 1;
		// 3-gram
		_trigram_counts[ti_2][ti_1][ti] += 1;
		_trigram_counts[ti_1][ti][ti1] += 1;
		_trigram_counts[ti][ti1][ti2] += 1;
		// 品詞-単語ペア
		increment_tag_word_count(ti, wi);
	}
	// in:  t_{i-2},t_{i-1},ti,t_{i+1},t_{i+2},w_i
	// out: void
	void remove_tag_from_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
		// 1-gram
		_unigram_counts[ti] -= 1;
		assert(_unigram_counts[ti] >= 0);
		// 2-gram
		_bigram_counts[ti_1][ti] -= 1;
		assert(_bigram_counts[ti_1][ti] >= 0);
		_bigram_counts[ti][ti1] -= 1;
		assert(_bigram_counts[ti][ti1] >= 0);
		// 3-gram
		_trigram_counts[ti_2][ti_1][ti] -= 1;
		assert(_trigram_counts[ti_2][ti_1][ti] >= 0);
		_trigram_counts[ti_1][ti][ti1] -= 1;
		assert(_trigram_counts[ti_1][ti][ti1] >= 0);
		_trigram_counts[ti][ti1][ti2] -= 1;
		assert(_trigram_counts[ti][ti1][ti2] >= 0);
		// 品詞-単語ペア
		decrement_tag_word_count(ti, wi);
	}
	void perform_gibbs_sampling_with_line(vector<Word*> &line){
		if(_sampling_table == NULL){
			_sampling_table = (double*)malloc(_num_tags * sizeof(double));
		}
		for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
			int ti_2 = line[pos - 2]->tag_id;
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			int ti1 = line[pos + 1]->tag_id;
			int ti2 = line[pos + 2]->tag_id;
			// t_iをモデルパラメータから除去
			remove_tag_from_model_parameters(ti_2, ti_1, ti, ti1, ti2, wi);
			// t_iを再サンプリング
			double sum = 0;
			int new_ti = 0;
			for(int tag = 0;tag < _num_tags;tag++){
				_sampling_table[tag] = 1;
				double n_ti_wi = get_count_for_tag_word(tag, wi);
				double n_ti = _unigram_counts[tag];
				double W_ti = _Wt[tag];
				double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][tag];
				double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
				double n_ti_1_ti_ti1 = _trigram_counts[ti_1][tag][ti1];
				double n_ti_1_ti = _bigram_counts[ti_1][tag];
				double I_ti_2_ti_1_ti_ti1 = (ti_2 == ti_1 == tag == ti1) ? 1 : 0;
				double I_ti_2_ti_1_ti = (ti_2 == ti_1 == tag) ? 1 : 0;
				double n_ti_ti1_ti2 = _trigram_counts[tag][ti1][ti2];
				double n_ti_ti1 = _bigram_counts[tag][ti1];
				double I_ti_2_ti_ti2_and_ti_1_ti1 = (ti_2 == tag == ti2 && ti_1 == ti1) ? 1 : 0;
				double I_ti_1_ti_ti1_ti2 = (ti_1 == tag == ti1 == ti2) ? 1 : 0;
				double I_ti_2_ti_and_ti_1_ti1 = (ti_2 == tag && ti_1 == ti1) ? 1 : 0;
				double I_ti_1_ti_ti1 = (ti_1 == tag == ti1) ? 1 : 0;
				_sampling_table[tag] *= (n_ti_wi + _beta[tag]) / (n_ti + W_ti * _beta[tag]);
				_sampling_table[tag] *= (n_ti_2_ti_1_ti + _alpha) / (n_ti_2_ti_1 + _num_tags * _alpha);
				_sampling_table[tag] *= (n_ti_1_ti_ti1 + I_ti_2_ti_1_ti_ti1 + _alpha) / (n_ti_1_ti + I_ti_2_ti_1_ti + _num_tags * _alpha);
				_sampling_table[tag] *= (n_ti_ti1_ti2 + I_ti_2_ti_ti2_and_ti_1_ti1 + I_ti_1_ti_ti1_ti2 + _alpha) / (n_ti_ti1 + I_ti_2_ti_and_ti_1_ti1 + I_ti_1_ti_ti1 + _num_tags * _alpha);
				_sampling_table[tag] = pow(_sampling_table[tag], 1.0 / _temperature);
				sum += _sampling_table[tag];
			}
			assert(sum > 0);
			double normalizer = 1.0 / sum;
			double bernoulli = Sampler::uniform(0, 1);
			sum = 0;
			for(int tag = 0;tag < _num_tags;tag++){
				sum += _sampling_table[tag] * normalizer;
				if(sum >= bernoulli){
					new_ti = tag;
					break;
				}
			}
			// 新しいt_iをモデルパラメータに追加
			add_tag_to_model_parameters(ti_2, ti_1, new_ti, ti1, ti2, wi);
			line[pos]->tag_id = new_ti;
		}
	}
	// 論文(6)式と(7)式を掛けたものからtiをサンプリング
	int sample_tag_from_Pt_w(int ti_2, int ti_1, int wi){
		double sum_p = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			double Pt_alpha = (_trigram_counts[ti_2][ti_1][tag] + _alpha) / (_bigram_counts[ti_2][ti_1] + _num_tags * _alpha);
			double Pw_t_beta = (get_count_for_tag_word(tag, wi) + _beta[tag]) / (_unigram_counts[tag] + _Wt[tag] * _beta[tag]);
			double Ptw_alpha_beta = Pw_t_beta * Pt_alpha;
			_sampling_table[tag] = Ptw_alpha_beta;
			sum_p += Ptw_alpha_beta;
		}
		double normalizer = 1.0 / sum_p;
		double bernoulli = Sampler::uniform(0, 1);
		sum_p = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			sum_p += _sampling_table[tag] * normalizer;
			if(sum_p > bernoulli){
				return tag;
			}
		}
		return _num_tags - 1;
	}
	// 論文(6)式と(7)式を掛けたものからtiをサンプリング
	int argmax_tag_from_Pt_w(int ti_2, int ti_1, int wi){
		double max_p = 0;
		double max_tag = 0;
		// cout << (boost::format("argmax(%d, %d, %d)") % ti_2 % ti_1 % wi).str() << endl;
		for(int tag = 0;tag < _num_tags;tag++){
			double Pt_alpha = (_trigram_counts[ti_2][ti_1][tag] + _alpha) / (_bigram_counts[ti_2][ti_1] + _num_tags * _alpha);
			double Pw_t_beta = (get_count_for_tag_word(tag, wi) + _beta[tag]) / (_unigram_counts[tag] + _Wt[tag] * _beta[tag]);
			double Ptw_alpha_beta = Pw_t_beta * Pt_alpha;
			// cout << (boost::format("%f = %f * %f") % Ptw_alpha_beta % Pw_t_beta % Pt_alpha).str() << endl;
			if(Ptw_alpha_beta > max_p){
				max_p = Ptw_alpha_beta;
				max_tag = tag;
			}
		}
		// cout << "return " << max_tag << endl;
		return max_tag;
	}
	Word* _get_random_word_with_tag(int tag, vector<vector<Word*>> &dataset){
		int random_index = Sampler::uniform_int(0, dataset.size() - 1);
		vector<Word*> &line = dataset[random_index];
		for(int pos = 0;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			if(ti == tag){
				return line[pos];
			}
		}
		return NULL;
	}
	// 新しいBetaをサンプリング
	void sample_new_beta(vector<vector<Word*>> &dataset){
		for(int tag = 0;tag < _num_tags;tag++){
			_new_beta[tag] = Sampler::normal(_beta[tag], 0.1 * _beta[tag]);
			Word* random_word = NULL;
			int limit = 100;
			while(random_word == NULL){
				if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
					return;
				}
				random_word = _get_random_word_with_tag(tag, dataset);
				limit--;
				if(limit < 0){
					break;
				}
			}
			if(random_word == NULL){
				continue;
			}
			double Pti_wi_beta = compute_Pti_wi_beta(random_word->tag_id, random_word->word_id, _beta[tag]);
			double Pti_wi_new_beta = compute_Pti_wi_beta(random_word->tag_id, random_word->word_id, _new_beta[tag]);
			if(Pti_wi_new_beta == 0){
				continue;
			}
			double adoption_rate = std::min(1.0, Pti_wi_new_beta / Pti_wi_beta);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli < adoption_rate){
				_beta[tag] = _new_beta[tag];
			}
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
	bool save(string dir = "out"){
		ofstream ofs(dir + "/hmm.obj");
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << static_cast<const BayesianHMM&>(*this);
		ofs.close();
		ofstream ofs_bin;
		// 3-gram
		ofs_bin.open(dir + "/hmm.trigram", ios::binary);
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				ofs_bin.write((char*)(_trigram_counts[tri_tag][bi_tag]), _num_tags * sizeof(int));
			}
		}
		ofs_bin.close();
		// 2-gram
		ofs_bin.open(dir + "/hmm.bigram", ios::binary);
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			ofs_bin.write((char*)(_bigram_counts[bi_tag]), _num_tags * sizeof(int));
		}
		ofs_bin.close();
		// 1-gram
		ofs_bin.open(dir + "/hmm.unigram", ios::binary);
		ofs_bin.write((char*)(_unigram_counts), _num_tags * sizeof(int));
		ofs_bin.close();
		// beta
		ofs_bin.open(dir + "/hmm.beta", ios::binary);
		ofs_bin.write((char*)(_beta), _num_tags * sizeof(double));
		ofs_bin.close();
		// Wt
		ofs_bin.open(dir + "/hmm.wt", ios::binary);
		ofs_bin.write((char*)(_Wt), _num_tags * sizeof(int));
		ofs_bin.close();
		return true;
	}
	bool load(string dir = "out"){
		bool complete = true;
		ifstream ifs(dir + "/hmm.obj");
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> *this;
		}else{
			complete = false;
		}
		ifs.close();

		// ポインタの読み込み
		if(_trigram_counts == NULL){
			alloc_table();
		}
		ifstream ifs_bin;
		// 3-gram
		ifs_bin.open(dir + "/hmm.trigram", ios::binary);
		if(ifs_bin.good()){
			for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
				for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
					ifs_bin.read((char*)(_trigram_counts[tri_tag][bi_tag]), _num_tags * sizeof(int));
				}
			}
		}else{
			complete = false;
		}
		ifs_bin.close();
		// 2-gram
		ifs_bin.open(dir + "/hmm.bigram", ios::binary);
		if(ifs_bin.good()){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				ifs_bin.read((char*)(_bigram_counts[bi_tag]), _num_tags * sizeof(int));
			}
		}else{
			complete = false;
		}
		ifs_bin.close();
		// 1-gram
		ifs_bin.open(dir + "/hmm.unigram", ios::binary);
		if(ifs_bin.good()){
			ifs_bin.read((char*)(_unigram_counts), _num_tags * sizeof(int));
		}else{
			complete = false;
		}
		ifs_bin.close();
		// beta
		ifs_bin.open(dir + "/hmm.beta", ios::binary);
		if(ifs_bin.good()){
			ifs_bin.read((char*)(_beta), _num_tags * sizeof(double));
		}else{
			complete = false;
		}
		ifs_bin.close();
		// Wt
		ifs_bin.open(dir + "/hmm.wt", ios::binary);
		if(ifs_bin.good()){
			ifs_bin.read((char*)(_Wt), _num_tags * sizeof(int));
		}else{
			complete = false;
		}
		ifs_bin.close();
		return complete;
	}
};

#endif