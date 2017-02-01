#ifndef _ihmm_
#define _ihmm_
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

class InfiniteHMM{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _tag_count;
		archive & _bigram_tag_counts;
		archive & _tag_word_counts;
		archive & _oracle_word_counts;
		archive & _oracle_tag_counts;
		archive & _alpha;
		archive & _beta;
		archive & _gamma;
		archive & _beta_emission;
		archive & _gamma_emission;
	}
public:
	vector<int> _tag_count;	// 全ての状態とそのカウント
	unordered_map<int, unordered_map<int, int>> _bigram_tag_counts;	// 品詞と単語のペアの出現頻度
	unordered_map<int, unordered_map<int, int>> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	unordered_map<int, int> _oracle_word_counts;	// 品詞と単語のペアの出現頻度
	unordered_map<int, int> _oracle_tag_counts;	// 品詞と単語のペアの出現頻度
	double _alpha;
	double _beta;
	double _gamma;
	double _beta_emission;
	double _gamma_emission;
	int _initial_num_tags;
	int _sum_oracle_tags_count;
	int _sum_oracle_words_count;
	unordered_map<int, int> _sum_bigram_destination;
	unordered_map<int, int> _sum_word_count_for_tag;
	InfiniteHMM(int initial_num_tags){
		_alpha = 1;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_initial_num_tags = initial_num_tags;
		_sum_oracle_words_count = 0;
		_sum_oracle_tags_count = 0;
	}
	void initialize(vector<vector<Word*>> &dataset){
		for(int tag = 0;tag < _initial_num_tags;tag++){
			_tag_count.push_back(0);
		}
		// nグラムのカウントテーブル
		init_ngram_counts(dataset);
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		c_printf("[*]%s\n", "n-gramモデルを構築してます ...");
		// 最初は品詞をランダムに割り当てる
		set<int> word_set;
		unordered_map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			word_set.insert(line[0]->word_id);
			increment_tag_word_count(line[0]);
			increment_tag_count(line[0]->tag_id);
			
			for(int pos = 1;pos < line.size();pos++){	// 2-gramなので3番目から.
				Word* word = line[pos];

				int ti_1 = line[pos - 1]->tag_id;
				int ti = (pos == line.size() - 1) ? 0 :  Sampler::uniform_int(0, _initial_num_tags - 1);
				int wi = word->word_id;
				double empirical_p, coeff_oracle_p;

				_compute_Ptag_context(ti, ti_1, empirical_p, coeff_oracle_p);
				double normalizer = 1 / (empirical_p + coeff_oracle_p);
				double bernoulli = Sampler::uniform(0, 1);
				if(bernoulli < empirical_p * normalizer){
					increment_tag_bigram_count(ti_1, ti);
					increment_tag_count(ti);
				}else{	// oracleから生成された場合
					increment_tag_bigram_count(ti_1, ti);
					increment_tag_count(ti);
					increment_oracle_tag_count(ti);
				}

				_compute_Pword_tag(wi, ti, empirical_p, coeff_oracle_p);
				normalizer = 1 / (empirical_p + coeff_oracle_p);
				bernoulli = Sampler::uniform(0, 1);
				if(bernoulli < empirical_p * normalizer){
					increment_tag_word_count(ti, wi);
				}else{	// oracleから生成された場合
					increment_tag_word_count(ti, wi);
					increment_oracle_word_count(wi);
				}

				word_set.insert(wi);
				word->tag_id = ti;
			}
		}
		c_printf("[*]%s\n", (boost::format("単語数: %d - 行数: %d") % word_set.size() % dataset.size()).str().c_str());
	}
	void increment_tag_count(int tag_id){
		while(tag_id >= _tag_count.size()){
			_tag_count.push_back(0);
		}
		_tag_count[tag_id] += 1;
	}
	void decrement_tag_count(int tag_id){
		while(tag_id >= _tag_count.size()){
			_tag_count.push_back(0);
		}
		_tag_count[tag_id] -= 1;
		assert(_tag_count[tag_id] >= 0);
	}
	void increment_tag_bigram_count(int context_tag_id, int tag_id){
		_bigram_tag_counts[context_tag_id][tag_id] += 1;
		_sum_bigram_destination[context_tag_id] += 1;
	}
	void increment_tag_word_count(Word* word){
		increment_tag_word_count(word->tag_id, word->word_id);
	}
	void increment_tag_word_count(int tag_id, int word_id){
		_sum_word_count_for_tag[tag_id] += 1;
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			word_counts[word_id] = 1;
			return;
		}
		itr->second += 1;
	}
	void increment_oracle_tag_count(int tag_id){
		_oracle_tag_counts[tag_id] += 1;
		_sum_oracle_tags_count += 1;
	}
	void increment_oracle_word_count(int word_id){
		_oracle_word_counts[word_id] += 1;
		_sum_oracle_words_count += 1;
	}
	void decrement_oracle_word_count(int word_id){
		_oracle_word_counts[word_id] -= 1;
		assert(_oracle_word_counts[word_id] >= 0);
		_sum_oracle_words_count -= 1;
		assert(_sum_oracle_words_count >= 0);
	}
	void decrement_oracle_tag_count(int tag_id){
		_oracle_tag_counts[tag_id] -= 1;
		assert(_oracle_tag_counts[tag_id] >= 0);
		_sum_oracle_tags_count -= 1;
		assert(_sum_oracle_tags_count >= 0);
	}
	void decrement_tag_bigram_count(int context_tag_id, int tag_id){
		_bigram_tag_counts[context_tag_id][tag_id] -= 1;
		assert(_bigram_tag_counts[context_tag_id][tag_id] >= 0);
		auto itr = _sum_bigram_destination.find(context_tag_id);
		assert(itr != _sum_bigram_destination.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_sum_bigram_destination.erase(itr);
		}
	}
	void decrement_tag_word_count(int tag_id, int word_id){
		unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			c_printf("[R]%s [*]%s", "エラー", "品詞-単語ペアのカウントが正しく実装されていません.");
			exit(1);
		}
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			word_counts.erase(itr);
		}
		auto itr_sum = _sum_word_count_for_tag.find(tag_id);
		assert(itr_sum != _sum_word_count_for_tag.end());
		itr_sum->second -= 1;
		assert(itr_sum->second >= 0);
		if(itr_sum->second == 0){
			_sum_word_count_for_tag.erase(itr_sum);
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
	bool is_tag_new(int tag_id){
		if(tag_id >= _tag_count.size()){
			return true;
		}
		if(_tag_count[tag_id] == 0){
			return true;
		}
		return false;
	}
	int get_new_tag_id(){
		for(int tag = 0;tag < _tag_count.size();tag++){
			if(_tag_count[tag] == 0){
				return tag;
			}
		}
		return _tag_count.size();
	}
	int sum_oracle_words_count(){
		return _sum_oracle_words_count;
		// int sum = 0;
		// for(const auto &count: _oracle_word_counts){
		// 	sum += count.second;
		// }
		// return sum;
	}
	int sum_word_count_for_tag(int tag_id){
		return _sum_word_count_for_tag[tag_id];
		// int sum = 0;
		// unordered_map<int, int> &counts = _tag_word_counts[tag_id];
		// for(const auto &count: counts){
		// 	sum += count.second;
		// }
		// return sum;
	}
	int sum_bigram_destination(int tag_id){
		return _sum_bigram_destination[tag_id];
		// auto itr = _sum_bigram_destination.find(tag_id);
		// if(itr == _sum_bigram_destination.end()){
		// 	int sum = 0;
		// 	unordered_map<int, int> &unigram_table = _bigram_tag_counts[tag_id];
		// 	for(const auto &unigram: unigram_table){
		// 		sum += unigram.second;
		// 	}
		// 	_sum_bigram_destination[tag_id] = sum;
		// 	return sum;
		// }
		// return itr->second;
	}
	int sum_oracle_tags_count(){
		return _sum_oracle_tags_count;
		// int sum = 0;
		// for(const auto &unigram: _oracle_tag_counts){
		// 	sum += unigram.second;
		// }
		// return sum;
	}
	void _compute_Ptag_context(int tag_id, int context_tag_id, double &empirical_p, double &coeff_oracle_p){
		empirical_p = 0;
		double n_i = sum_bigram_destination(context_tag_id);
		if(is_tag_new(tag_id) == false){
			double n_ij = _bigram_tag_counts[context_tag_id][tag_id];
			empirical_p = n_ij / (n_i + _beta);
		}
		coeff_oracle_p = _beta / (n_i + _beta);
	}
	// P(s_{t+1}|s_t)
	double compute_Ptag_context(int tag_id, int context_tag_id){
		double empirical_p = 0;		// 自分が生成したサンプルからなる経験確率.新しい品詞の場合は存在しないので0.
		double n_i = sum_bigram_destination(context_tag_id);
		if(is_tag_new(tag_id) == false){
			double n_ij = _bigram_tag_counts[context_tag_id][tag_id];
			empirical_p = n_ij / (n_i + _beta);
		}
		double coeff_oracle_p = _beta / (n_i + _beta);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		double oracle_p;
		double n_o = sum_oracle_tags_count();
		if(is_tag_new(tag_id)){
			oracle_p = _gamma / (n_o + _gamma);
		}else{
			double n_oj = _oracle_tag_counts[tag_id];
			oracle_p = n_oj / (n_o + _gamma);
		}
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	void _compute_Pword_tag(int word_id, int tag_id, double &empirical_p, double &coeff_oracle_p){
		empirical_p = 0;
		double m_i = sum_word_count_for_tag(tag_id);
		if(is_tag_new(tag_id) == false){
			double m_iq = _tag_word_counts[tag_id][word_id];
			empirical_p = m_iq / (m_i + _beta_emission);
		}
		coeff_oracle_p = _beta_emission / (m_i + _beta_emission);
	}
	// P(y_t|s_t)
	double compute_Pword_tag(int word_id, int tag_id){
		double empirical_p = 0;		// 自分が生成したサンプルからなる経験確率.新しい品詞の場合は存在しないので0.
		double m_i = sum_word_count_for_tag(tag_id);
		if(is_tag_new(tag_id) == false){
			double m_iq = _tag_word_counts[tag_id][word_id];
			empirical_p = m_iq / (m_i + _beta_emission);
		}
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);	// 親の分布から生成される確率. 親からword_idが生成される確率とは別物.
		double oracle_p;
		double m_o = sum_oracle_words_count();
		if(is_tag_new(tag_id)){
			oracle_p = _gamma_emission / (m_o + _gamma_emission);
		}else{
			double m_oq = _oracle_word_counts[word_id];
			oracle_p = m_oq / (m_o + _gamma_emission);
		}
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	void remove_tag_from_model(int context_tag_id, int tag_id){
		assert(is_tag_new(tag_id) == false);
		double n_ij = _bigram_tag_counts[context_tag_id][tag_id];
		double n_i = sum_bigram_destination(context_tag_id);
		double empirical_p = n_ij / (n_i + _beta);
		double coeff_oracle_p = _beta / (n_i + _beta);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		double n_o = sum_oracle_tags_count();
		double n_oj = _oracle_tag_counts[tag_id];
		if(n_oj == 0){
			decrement_tag_bigram_count(context_tag_id, tag_id);
			return;
		}
		// double oracle_p = n_oj / (n_o + _gamma);
		double normalizer = 1 / (empirical_p + coeff_oracle_p);
		double bernoulli = Sampler::uniform(0, 1);
		if(bernoulli < empirical_p * normalizer){
			decrement_tag_bigram_count(context_tag_id, tag_id);
		}else{
			decrement_oracle_tag_count(tag_id);
		}
	}
	void remove_word_from_model(int word_id, int tag_id){
		assert(is_tag_new(tag_id) == false);
		double m_iq = _tag_word_counts[tag_id][word_id];
		double m_i = sum_word_count_for_tag(tag_id);
		double empirical_p = m_iq / (m_i + _beta_emission);
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);
		double m_oq = _oracle_word_counts[word_id];
		double m_o = sum_oracle_words_count();
		// double oracle_p = m_oq / (m_o + _gamma_emission);
		double normalizer = 1 / (empirical_p + coeff_oracle_p);
		double bernoulli = Sampler::uniform(0, 1);
		if(bernoulli < empirical_p * normalizer){
			decrement_tag_word_count(tag_id, word_id);
		}else{
			decrement_oracle_word_count(word_id);
		}
	}
	// t_{i-1} -> t_i -> t_{i+1}
	int sampling_new_tag(int ti_1, int ti, int ti1, int wi){
		cout << "sampling_new_tag" << endl;
		// ギブスサンプリング
		vector<double> sampling_table;
		bool new_tag_included = false;
		double sum = 0;
		for(int tag = 0;tag < _tag_count.size();tag++){
			if(is_tag_new(tag)){
				new_tag_included = true;
			}
			double p_emission = compute_Pword_tag(wi, ti);
			double p_generation = compute_Ptag_context(tag, ti_1);
			double p_likelihood = compute_Ptag_context(ti1, tag);
			double p_conditional = p_emission * p_generation * p_likelihood;
			cout << p_emission << ",";
			cout << p_generation << ",";
			cout << p_likelihood << ",";
			cout << endl;
			sampling_table.push_back(p_conditional);
			sum += p_conditional;
		}
		int new_tag = _tag_count.size();
		if(new_tag_included == false){
			double p_emission = compute_Pword_tag(wi, ti);
			double p_generation = compute_Ptag_context(new_tag, ti_1);
			double p_likelihood = compute_Ptag_context(ti1, new_tag);
			double p_conditional = p_emission * p_generation * p_likelihood;
			sampling_table.push_back(p_conditional);
			sum += p_conditional;
		}
		assert(sum > 0);
		double normalizer = 1 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int tag = 0;tag < _tag_count.size();tag++){
			sum += sampling_table[tag] * normalizer;
			if(bernoulli < sum){
				return tag;
			}
		}
		cout << sum << endl;
		cout << "new_tag generated" << endl;
		return new_tag;
	}
	void perform_gibbs_sampling_with_line(vector<Word*> &line){
		for(int pos = 1;pos < line.size() - 1;pos++){
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			int ti1 = line[pos + 1]->tag_id;

			// 現在のtiをモデルから除去
			// remove_word_from_model(wi, ti);	// 先に品詞-単語ペアから除去するとassrtで引っかからない
			remove_tag_from_model(ti_1, ti);
			remove_tag_from_model(ti, ti1);
			decrement_tag_count(ti);
			decrement_tag_word_count(ti, wi);

			int new_tag = sampling_new_tag(ti_1, ti, ti1, wi);
			// モデルに追加
			increment_tag_word_count(new_tag, wi);
			increment_tag_bigram_count(ti_1, new_tag);
			increment_tag_bigram_count(new_tag, ti1);
			increment_tag_count(new_tag);
			if(is_tag_new(new_tag)){
				increment_oracle_tag_count(new_tag);
				// increment_oracle_word_count(wi);
			}
			line[pos]->tag_id = new_tag;
		}
	}
	void dump_tags(){
		for(int tag = 0;tag < _tag_count.size();tag++){
			cout << _tag_count[tag] << ", ";
		}
		cout << endl;
	}
	bool load(string dir = "out"){
		return true;
	}
	bool save(string dir = "out"){
		return true;
	}

};

#endif