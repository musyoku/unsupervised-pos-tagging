#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <unordered_map>
#include <vector>
#include "common.h"
#include "sampler.h"
#include "util.h"
#include "table.h"

// <s>と</s>のIDは0

namespace ihmm {
	class InfiniteHMM {
	public:
		int _initial_num_tags;
		int _num_words;
		int _prev_num_tags;
		std::vector<std::vector<Table*>> _n_ij_tables;	// 品詞bigramの出現頻度
		std::vector<Table**> _m_iq_tables;	// 品詞と単語のペアの出現頻度
		int* _oracle_m_q_counts;
		int _oracle_sum_n_over_j;			// \sum_j{n_j^o}
		int _oracle_sum_m_over_q;			// \sum_j{m_q^o}
		std::vector<int> _sum_n_i_over_j;	// \sum_j{n_ij}の計算用
		std::vector<int> _oracle_n_j_counts;
		std::vector<int> _sum_m_i_over_q;
		double* _gibbs_sampling_table;
		double _alpha;
		double _beta;
		double _gamma;
		double _beta_emission;
		double _gamma_emission;
		InfiniteHMM();
		InfiniteHMM(int initial_num_tags, int num_words);
		~InfiniteHMM();
		void initialize_with_training_dataset(std::vector<std::vector<Word*>> &dataset);
		void _remove_all_training_dataset(std::vector<std::vector<Word*>> &dataset);
		int get_num_tags();
		int get_num_valid_tags();
		int get_num_words();
		int get_sum_n_i_over_j(int tag);
		int get_n_ij(int context_tag, int tag);
		int get_oracle_sum_n_over_j();
		int get_oracle_n_j(int tag);
		int get_sum_m_i_over_q(int tag);
		int get_m_iq(int tag, int word_id);
		int get_oracle_sum_m_over_q();
		int get_oracle_m_q(int word_id);
		int _get_new_tag();
		bool is_tag_new(int tag);
		double compute_p_tag_given_context(int tag, int context_tag);
		double compute_p_word_given_tag(int word_id, int tag);
		bool save(std::string filename);
		bool load(std::string filename);
		// void _add_special_tag();
		int _add_new_tag();
		void _delete_tag(int tag);
		void _increment_tag_bigram_count(int context_tag, int tag);
		void _decrement_tag_bigram_count(int context_tag, int tag);
		void _increment_tag_word_count(int tag, int word_id);
		void _decrement_tag_word_count(int tag, int word_id);
		void _increment_oracle_tag_count(int tag);
		void _decrement_oracle_tag_count(int tag);
		void _increment_oracle_word_count(int word_id);
		void _decrement_oracle_word_count(int word_id);
		int _perform_gibbs_sampling_on_markov_blanket(int ti_1, int ti1, int wi);
		void perform_gibbs_sampling_with_sequence(std::vector<Word*> &word_vec);
	};
}

// class InfiniteHMM{
// 	bool is_tag_new(int tag_id){
// 		assert(tag_id != EOS);
// 		if(tag_id >= _tag_unigram_count.size()){
// 			return true;
// 		}
// 		if(_tag_unigram_count[tag_id] == 0){
// 			return true;
// 		}
// 		return false;
// 	}
// 	int get_new_tag_id(){
// 		for(int tag = EOS + 1;tag < _tag_unigram_count.size();tag++){
// 			if(_tag_unigram_count[tag] == 0){
// 				return tag;
// 			}
// 		}
// 		return _tag_unigram_count.size();
// 	}
// 	bool is_word_new(int word_id){
// 		auto itr = _oracle_m_q_counts.find(word_id);
// 		if(itr == _oracle_m_q_counts.end()){
// 			return false;
// 		}
// 		return itr->second == 0;
// 	}
// 	double compute_gamma_distribution(double v, double a, double b){
// 		return pow(b, a) / tgamma(a) * pow(v, a - 1) * exp(-b * v);
// 	}
// 	double compute_log_Pdata(vector<Word*> &line){
// 		double p = 0;
// 		int ti_1 = BOS;
// 		for(int pos = 0;pos < line.size();pos++){
// 			int ti = line[pos]->_tag;
// 			int wi = line[pos]->_id;
// 			double log_p_tag = log(compute_Ptag_context(ti, ti_1));
// 			double log_p_word = log(compute_Pword_tag(wi, ti));
// 			p += log_p_tag + log_p_word;
// 			ti_1 = ti;
// 		}
// 		return p;
// 	}
// 	// t_{i-1} -> t_i -> t_{i+1}
// 	int gibbs_sample_new_tag(int ti_1, int ti1, int wi){
// 		// ギブスサンプリング
// 		double sum = 0;
// 		for(int tag = EOS + 1;tag < _tag_unigram_count.size();tag++){
// 			if(is_tag_new(tag)){
// 				_gibbs_sampling_table[tag] = 0;
// 				continue;
// 			}
// 			double p_emission = compute_Pword_tag(wi, tag);
// 			double p_generation = compute_Ptag_context(tag, ti_1);
// 			int correcting_count_for_bigram = (ti_1 == tag == ti1) ? 1 : 0;
// 			int correcting_count_for_destination = (ti_1 == tag) ? 1 : 0;
// 			double p_likelihood = compute_Ptag_context(ti1, tag, correcting_count_for_bigram, correcting_count_for_destination);
// 			double p_conditional = p_emission * p_generation * p_likelihood;
// 			_gibbs_sampling_table[tag] = p_conditional;
// 			sum += p_conditional;
// 		}
// 		int new_tag = get_new_tag_id();
// 		double p_emission = compute_Pword_tag(wi, new_tag);
// 		double p_generation = compute_Ptag_context(new_tag, ti_1);
// 		double p_likelihood = compute_Ptag_context(ti1, new_tag);
// 		double p_conditional = p_emission * p_generation * p_likelihood;
// 		sum += p_conditional;
// 		// new_tag > _tag_unigram_count.size()ならサンプリングテーブルに入れなくてもよい.
// 		if(new_tag < _tag_unigram_count.size()){
// 			_gibbs_sampling_table[new_tag] = p_conditional;
// 		}
// 		assert(sum > 0);
// 		double normalizer = 1.0 / sum;
// 		double bernoulli = Sampler::uniform(0, 1);
// 		sum = 0;
// 		for(int tag = EOS + 1;tag < _tag_unigram_count.size();tag++){
// 			sum += _gibbs_sampling_table[tag] * normalizer;
// 			if(bernoulli <= sum){
// 				return tag;
// 			}
// 		}
// 		return new_tag;
// 	}
// 	int argmax_Ptag_given_context_word(int context_tag_id, int word_id){
// 		double max_p = 0;
// 		double max_tag = 0;
// 		for(int tag = EOS + 1;tag < _tag_unigram_count.size();tag++){
// 			if(is_tag_new(tag)){
// 				continue;
// 			}
// 			double Ptag = compute_Ptag_context(tag, context_tag_id);
// 			double Pword = compute_Pword_tag(word_id, tag);
// 			double p = Ptag * Pword;
// 			if(p > max_p){
// 				max_p = p;
// 				max_tag = tag;
// 			}
// 		}
// 		return max_tag;
// 	}
// 	void perform_gibbs_sampling_line(vector<Word*> &line){
// 		if(line.size() < 2){
// 			return;
// 		}
// 		int pos = 0;
// 		int ti_1 = BOS;
// 		for(;pos < line.size();pos++){
// 			int ti = line[pos]->_tag;
// 			int wi = line[pos]->_id;
// 			int ti1 = (pos == line.size() - 1) ? EOS : line[pos + 1]->_tag;

// 			// 現在のtiをモデルから除去
// 			decrement_tag_bigram_count(ti_1, ti);
// 			decrement_tag_bigram_count(ti, ti1);
// 			decrement_tag_unigram_count(ti);
// 			decrement_tag_word_count(ti, wi);
// 			increment_tag_bigram_count(ti_1, ti1);
// 			// 新しい状態をサンプリング
// 			int new_tag = gibbs_sample_new_tag(ti_1, ti1, wi);
// 			// モデルに追加
// 			increment_tag_word_count(new_tag, wi);
// 			increment_tag_bigram_count(ti_1, new_tag);
// 			increment_tag_bigram_count(new_tag, ti1);
// 			increment_tag_unigram_count(new_tag);
// 			decrement_tag_bigram_count(ti_1, ti1);
// 			line[pos]->_tag = new_tag;
// 			ti_1 = new_tag;
// 		}
// 	}
// 	void check_oracle_tag_count(){
// 		unordered_map<int, int> counts;
// 		for(const auto &contexts: _n_ij_tables){
// 			for(const auto &tags: contexts.second){
// 				vector<int> &table = tags.second->_arrangement;
// 				int num_tables = table.size();
// 				counts[tags.first] += num_tables;
// 			}
// 		}
// 		for(const auto &elem: counts){
// 			int tag_id = elem.first;
// 			int num_tables = elem.second;
// 			int count = get_oracle_count_for_tag(tag_id);
// 			assert(num_tables == count);
// 		}
// 	}
// 	void check_oracle_word_count(){
// 		unordered_map<int, int> counts;
// 		for(const auto &tags: _m_iq_tables){
// 			for(const auto &words: tags.second){
// 				vector<int> &table = words.second->_arrangement;
// 				int num_tables = table.size();
// 				counts[words.first] += num_tables;
// 			}
// 		}
// 		for(const auto &elem: counts){
// 			int word_id = elem.first;
// 			int num_tables = elem.second;
// 			int count = get_oracle_count_for_word(word_id);
// 			assert(num_tables == count);
// 		}
// 	}
// 	void check_sum_bigram_destination(){
// 		for(const auto &tag: _sum_n_i_over_j){
// 			int tag_id = tag.first;
// 			int count = tag.second;
// 			int sum = 0;
// 			for(const auto &table: _n_ij_tables[tag_id]){
// 				sum += table.second->_num_customers;
// 			}
// 			assert(count == sum);
// 		}
// 	}
// 	void check_tag_count(){
// 		int num_non_zero = 0;
// 		for(int i = 0;i < _tag_unigram_count.size();i++){
// 			if(_tag_unigram_count[i] > 0){
// 				num_non_zero += 1;
// 			}
// 		}
// 		assert(num_non_zero == _oracle_n_j_counts.size());
// 	}
// 	void check_sum_word_customers(){
// 		int num_customers = 0;
// 		for(const auto &tag: _m_iq_tables){
// 			for(const auto &word: tag.second){
// 				Table* table = word.second;
// 				num_customers += table->_num_customers;
// 			}
// 		}
// 		assert(num_customers == _num_words);
// 	}
// 	void check_sum_tag_customers(){
// 		int num_customers_in_bigram = 0;
// 		for(const auto &bigram: _n_ij_tables){
// 			for(const auto &unigram: bigram.second){
// 				Table* table = unigram.second;
// 				num_customers_in_bigram += table->_num_customers;
// 			}
// 		}
// 		int num_customers_in_unigram = 0;
// 		for(auto itr = _tag_unigram_count.begin();itr != _tag_unigram_count.end();itr++){
// 			num_customers_in_unigram += *itr;
// 		}
// 		assert(num_customers_in_bigram == num_customers_in_unigram);
// 	}
// 	bool save(string dir = "out"){
// 		bool success = false;
// 		ofstream ofs(dir + "/ihmm.model");
// 		if(ofs.good()){
// 			boost::archive::binary_oarchive oarchive(ofs);
// 			oarchive << static_cast<const InfiniteHMM&>(*this);
// 			success = true;
// 		}
// 		ofs.close();
// 		return success;
// 	}
// 	bool load(string dir = "out"){
// 		bool success = false;
// 		ifstream ifs(dir + "/ihmm.model");
// 		if(ifs.good()){
// 			boost::archive::binary_iarchive iarchive(ifs);
// 			iarchive >> *this;
// 			success = true;
// 		}
// 		ifs.close();
// 		return success;
// 	}

// };