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
		archive & _all_states_count;
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
	vector<int> _all_states_count;	// 全ての状態とそのカウント
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
	InfiniteHMM(int initial_num_tags){
		_alpha = 1;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_initial_num_tags = initial_num_tags;
	}
	void initialize(vector<vector<Word*>> &dataset){
		for(int tag = 0;tag < _initial_num_tags;tag++){
			_all_states_count.push_back(0);
		}
		// nグラムのカウントテーブル
		init_ngram_counts(dataset);
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset, bool assign_different_tags_to_the_same_word = true){
		c_printf("[*]%s\n", "n-gramモデルを構築してます ...");
		// 最初は品詞をランダムに割り当てる
		set<int> word_set;
		unordered_map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			word_set.insert(line[0]->word_id);
			increment_tag_word_count(line[0]);
			for(int pos = 1;pos < line.size();pos++){	// 2-gramなので3番目から.
				Word* word = line[pos];
				auto itr = tag_for_word.find(word->word_id);
				if(itr == tag_for_word.end() || assign_different_tags_to_the_same_word){
					word->tag_id = Sampler::uniform_int(0, _initial_num_tags - 1);
					tag_for_word[word->word_id] = word->tag_id;
				}else{
					word->tag_id = itr->second;
				}
				word_set.insert(word->word_id);
				increment_tag_bigram_count(line[pos - 1], word);
				increment_tag_word_count(word);
			}
		}
		c_printf("[*]%s\n", (boost::format("単語数: %d - 行数: %d") % word_set.size() % dataset.size()).str().c_str());
	}
	void increment_tag_bigram_count(Word* bi_word, Word* uni_word){
		_bigram_tag_counts[bi_word->tag_id][uni_word->tag_id] += 1;
	}
	void increment_tag_word_count(Word* word){
		increment_tag_word_count(word->tag_id, word->word_id);
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
	bool is_tag_new(int tag_id){
		if(tag_id >= _all_states_count.size()){
			return true;
		}
		if(_all_states_count[tag_id] == 0){
			return true;
		}
		return false;
	}
	int get_new_tag_id(){
		for(int tag = 0;tag < _all_states_count.size();tag++){
			if(_all_states_count[tag] == 0){
				return tag;
			}
		}
		return _all_states_count.size();
	}
	int sum_oracle_words_count(){
		int sum = 0;
		for(const auto &count: _oracle_word_counts){
			sum += count.second;
		}
		return sum;
	}
	int sum_word_count_for_tag(int tag_id){
		int sum = 0;
		unordered_map<int, int> &counts = _tag_word_counts[tag_id];
		for(const auto &count: counts){
			sum += count.second;
		}
		return sum;
	}
	int sum_bigram_destination(int tag_id){
		int sum = 0;
		unordered_map<int, int> &unigram_table = _bigram_tag_counts[tag_id];
		for(const auto &unigram: unigram_table){
			sum += unigram.second;
		}
		return sum;
	}
	int sum_oracle_tags_count(){
		int sum = 0;
		for(const auto &unigram: _oracle_tag_counts){
			sum += unigram.second;
		}
		return sum;
	}
	// P(s_{t+1}|s_t)
	double compute_Ptag_context(int tag_id, int context_id){
		double empirical_p = 0;		// 自分が生成したサンプルからなる経験確率.新しい品詞の場合は存在しないので0.
		double n_i = sum_bigram_destination(context_id);
		if(is_tag_new(tag_id) == false){
			double n_ij = _bigram_tag_counts[context_id][tag_id];
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
	bool load(string dir = "out"){
		return true;
	}
	bool save(string dir = "out"){
		return true;
	}

};

#endif