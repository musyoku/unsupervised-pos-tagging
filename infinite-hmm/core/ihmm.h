#ifndef _ihmm_
#define _ihmm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "cprintf.h"
#include "sampler.h"
#include "util.h"
using namespace std;

typedef struct Word {
	int word_id;
	int tag_id;
} Word;

class Table{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _arrangement;
		archive & _num_customers;
		archive & _token_id;
	}
public:
	vector<int> _arrangement;
	int _num_customers;
	int _token_id;
	Table(){
		_num_customers = 0;
		_token_id = 0;
	}
	Table(int token_id){
		_num_customers = 0;
		_token_id = token_id;
	}
	bool is_empty(){
		return _arrangement.size() == 0;
	}
	void add_customer(double concentration_parameter, bool &new_table_generated){
		_num_customers += 1;
		if(_arrangement.size() == 0){
			_arrangement.push_back(1);
			new_table_generated = true;
			return;
		}
		new_table_generated = false;
		double sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0) + concentration_parameter;
		double normalizer = 1 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int i = 0;i < _arrangement.size();i++){
			sum += _arrangement[i] * normalizer;
			if(bernoulli <= sum){
				_arrangement[i] += 1;
				return;
			}
		}
		_arrangement.push_back(1);
		new_table_generated = true;
	}
	void remove_customer(bool &empty_table_deleted){
		assert(_arrangement.size() > 0);
		empty_table_deleted = false;
		_num_customers -= 1;
		int sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
		int bernoulli = Sampler::uniform_int(0, sum);
		sum = 0;
		int target_index = _arrangement.size() - 1;
		for(int i = 0;i < _arrangement.size();i++){
			sum += _arrangement[i];
			if(bernoulli <= sum){
				target_index = i;
				break;
			}
		}
		_arrangement[target_index] -= 1;
		if(_arrangement[target_index] == 0){
			_arrangement.erase(_arrangement.begin() + target_index);
			empty_table_deleted = true;
		}
	}
};

class InfiniteHMM{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _tag_unigram_count;
		archive & _bigram_tag_table;
		archive & _tag_word_table;
		archive & _oracle_word_counts;
		archive & _oracle_tag_counts;
		archive & _alpha;
		archive & _beta;
		archive & _gamma;
		archive & _beta_emission;
		archive & _gamma_emission;
		archive & _initial_num_tags;
		archive & _sum_oracle_tags_count;
		archive & _sum_oracle_words_count;
		archive & _num_words;
		archive & _sum_bigram_destination;
		archive & _sum_word_count_for_tag;
	}
public:
	vector<int> _tag_unigram_count;	// 全ての状態とそのカウント
	unordered_map<int, unordered_map<int, Table*>> _bigram_tag_table;	// 品詞と単語のペアの出現頻度
	unordered_map<int, unordered_map<int, Table*>> _tag_word_table;	// 品詞と単語のペアの出現頻度
	unordered_map<int, int> _oracle_word_counts;	// 品詞と単語のペアの出現頻度
	unordered_map<int, int> _oracle_tag_counts;	// 品詞と単語のペアの出現頻度
	double _alpha;
	double _beta;
	double _gamma;
	double _beta_emission;
	double _gamma_emission;
	int _initial_num_tags;
	int _num_words;
	int _sum_oracle_tags_count;
	int _sum_oracle_words_count;
	double _temperature;
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
		_num_words = 0;
		_temperature = 1;
	}
	void initialize(vector<vector<Word*>> &dataset){
		for(int tag = 0;tag < _initial_num_tags;tag++){
			_tag_unigram_count.push_back(0);
		}
		// nグラムのカウントテーブル
		init_ngram_counts(dataset);
	}
	void init_ngram_counts(vector<vector<Word*>> &dataset){
		c_printf("[*]%s\n", "n-gramモデルを構築してます ...");
		// 最初は品詞をランダムに割り当てる
		set<int> word_set;
		int num_words = 0;
		unordered_map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			word_set.insert(line[0]->word_id);
			// increment_tag_word_count(line[0]);
			// increment_tag_unigram_count(line[0]->tag_id);
			
			for(int pos = 1;pos < line.size();pos++){	// 2-gramなので3番目から.
				Word* word = line[pos];

				int ti_1 = line[pos - 1]->tag_id;
				int ti = (pos == line.size() - 1) ? 0 :  Sampler::uniform_int(0, _initial_num_tags - 1);
				// int ti = Sampler::uniform_int(0, _initial_num_tags - 1);
				int wi = word->word_id;
				increment_tag_bigram_count(ti_1, ti);
				increment_tag_unigram_count(ti);
				increment_tag_word_count(ti, wi);

				word_set.insert(wi);
				num_words += 1;
				word->tag_id = ti;
			}
		}
		c_printf("[*]%s\n", (boost::format("単語数: %d - 単語種: %d - 行数: %d") % num_words % word_set.size() % dataset.size()).str().c_str());
		_num_words += num_words;
	}
	void increment_tag_unigram_count(int tag_id){
		while(tag_id >= _tag_unigram_count.size()){
			_tag_unigram_count.push_back(0);
		}
		_tag_unigram_count[tag_id] += 1;
	}
	void decrement_tag_unigram_count(int tag_id){
		assert(tag_id < _tag_unigram_count.size());
		_tag_unigram_count[tag_id] -= 1;
		assert(_tag_unigram_count[tag_id] >= 0);
	}
	void increment_tag_bigram_count(int context_tag_id, int tag_id){
		_sum_bigram_destination[context_tag_id] += 1;

		Table* table = NULL;
		auto itr_context = _bigram_tag_table.find(context_tag_id);
		if(itr_context == _bigram_tag_table.end()){
			table = new Table(tag_id);
			_bigram_tag_table[context_tag_id][tag_id] = table;
		}else{
			unordered_map<int, Table*> &tables = itr_context->second;
			auto itr_table = tables.find(tag_id);
			if(itr_table == tables.end()){
				table = new Table(tag_id);
				tables[tag_id] = table;
			}else{
				table = itr_table->second;
			}
		}
		bool new_table_generated = false;
		table->add_customer(_beta, new_table_generated);
		if(new_table_generated){
			increment_oracle_tag_count(tag_id);
		}
	}
	void increment_tag_word_count(Word* word){
		increment_tag_word_count(word->tag_id, word->word_id);
	}
	void increment_tag_word_count(int tag_id, int word_id){
		_sum_word_count_for_tag[tag_id] += 1;

		Table* table = NULL;
		auto itr_tag = _tag_word_table.find(tag_id);
		if(itr_tag == _tag_word_table.end()){
			table = new Table(word_id);
			_tag_word_table[tag_id][word_id] = table;
		}else{
			unordered_map<int, Table*> &tables = itr_tag->second;
			auto itr_table = tables.find(word_id);
			if(itr_table == tables.end()){
				table = new Table(word_id);
				tables[word_id] = table;
			}else{
				table = itr_table->second;
			}
		}
		bool new_table_generated = false;
		table->add_customer(_beta_emission, new_table_generated);
		if(new_table_generated){
			increment_oracle_word_count(word_id);
		}
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
		auto itr = _oracle_tag_counts.find(tag_id);
		assert(itr != _oracle_tag_counts.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_oracle_tag_counts.erase(itr);
		}
		_sum_oracle_tags_count -= 1;
		assert(_sum_oracle_tags_count >= 0);
	}
	void decrement_tag_bigram_count(int context_tag_id, int tag_id){
		auto itr = _sum_bigram_destination.find(context_tag_id);
		assert(itr != _sum_bigram_destination.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_sum_bigram_destination.erase(itr);
		}

		auto itr_context = _bigram_tag_table.find(context_tag_id);
		assert(itr_context != _bigram_tag_table.end());
		unordered_map<int, Table*> &tables = itr_context->second;
		auto itr_table = tables.find(tag_id);
		assert(itr_table != tables.end());
		Table* table = itr_table->second;
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			decrement_oracle_tag_count(tag_id);
		}
		if(table->is_empty()){
			tables.erase(itr_table);
		}
		if(tables.size() == 0){
			_bigram_tag_table.erase(itr_context);
		}
	}
	void decrement_tag_word_count(int tag_id, int word_id){
		auto itr_sum = _sum_word_count_for_tag.find(tag_id);
		assert(itr_sum != _sum_word_count_for_tag.end());
		itr_sum->second -= 1;
		assert(itr_sum->second >= 0);
		if(itr_sum->second == 0){
			_sum_word_count_for_tag.erase(itr_sum);
		}

		auto itr_tag = _tag_word_table.find(tag_id);
		assert(itr_tag != _tag_word_table.end());
		unordered_map<int, Table*> &tables = itr_tag->second;
		auto itr_table = tables.find(word_id);
		assert(itr_table != tables.end());
		Table* table = itr_table->second;
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			decrement_oracle_word_count(word_id);
		}
		if(table->is_empty()){
			tables.erase(itr_table);
		}
		if(tables.size() == 0){
			_tag_word_table.erase(itr_tag);
		}
	}
	int get_bigram_tag_count(int context_tag_id, int tag_id){
		auto itr_context = _bigram_tag_table.find(context_tag_id);
		if(itr_context == _bigram_tag_table.end()){
			return 0;
		}
		unordered_map<int, Table*> &tables = itr_context->second;
		auto itr_table = tables.find(tag_id);
		if(itr_table == tables.end()){
			return 0;
		}
		Table* table = itr_table->second;
		return table->_num_customers;
	}
	int get_oracle_count_for_tag(int tag_id){
		auto itr = _oracle_tag_counts.find(tag_id);
		if(itr == _oracle_tag_counts.end()){
			return 0;
		}
		return itr->second;
	}
	int get_oracle_count_for_word(int word_id){
		auto itr = _oracle_word_counts.find(word_id);
		if(itr == _oracle_word_counts.end()){
			return 0;
		}
		return itr->second;
	}
	int get_tag_word_count(int tag_id, int word_id){
		auto itr_tag = _tag_word_table.find(tag_id);
		if(itr_tag == _tag_word_table.end()){
			return 0;
		}
		unordered_map<int, Table*> &tables = itr_tag->second;
		auto itr_table = tables.find(word_id);
		if(itr_table == tables.end()){
			return 0;
		}
		Table* table = itr_table->second;
		return table->_num_customers;
	}
	int get_num_times_oracle_tag_used(){
		int count = 0;
		for(const auto &tags: _bigram_tag_table){
			for(const auto &words: tags.second){
				count += words.second->_arrangement.size();
			}
		}
		return count;
	}
	int get_num_times_oracle_word_used(){
		int count = 0;
		for(const auto &tags: _tag_word_table){
			for(const auto &words: tags.second){
				count += words.second->_arrangement.size();
			}
		}
		return count;
	}
	int get_num_tags(){
		return _oracle_tag_counts.size();
	}
	int get_num_words(){
		return _num_words;
	}
	bool is_tag_new(int tag_id){
		if(tag_id >= _tag_unigram_count.size()){
			return true;
		}
		if(_tag_unigram_count[tag_id] == 0){
			return true;
		}
		return false;
	}
	int get_new_tag_id(){
		for(int tag = 0;tag < _tag_unigram_count.size();tag++){
			if(_tag_unigram_count[tag] == 0){
				return tag;
			}
		}
		return _tag_unigram_count.size();
	}
	bool is_word_new(int word_id){
		auto itr = _oracle_word_counts.find(word_id);
		if(itr == _oracle_word_counts.end()){
			return false;
		}
		return itr->second == 0;
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
		// 	unordered_map<int, int> &unigram_table = _bigram_tag_table[tag_id];
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
	// P(s_{t+1}|s_t)
	double compute_Ptag_context(int tag_id, int context_tag_id, int correcting_count_for_bigram = 0, int correcting_count_for_destination = 0){
		double n_i = sum_bigram_destination(context_tag_id);
		double n_ij = get_bigram_tag_count(context_tag_id, tag_id);
		// if(is_tag_new(tag_id) == false){
		// 	assert(n_ij > 0);
		// }
		double alpha = (tag_id == context_tag_id) ? _alpha : 0;
		double empirical_p = (n_ij + alpha) / (n_i + _beta + _alpha);
		double coeff_oracle_p = _beta / (n_i + _beta);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		double n_o = sum_oracle_tags_count();
		double n_oj = get_oracle_count_for_tag(tag_id);
		double T = get_num_tags();
		double g0 = 1.0 / (T + 1.0);
		double oracle_p = (n_oj + _gamma) / (n_o + _gamma);
		// if(is_tag_new(tag_id)){
		// 	oracle_p = _gamma / (n_o + _gamma);
		// }else{
		// 	double n_oj = get_oracle_count_for_tag(tag_id);
		// 	// double T = get_num_tags();
		// 	// return (n_ij + _beta * g0) / (n_i + _beta);
		// 	// double numerator = is_tag_new(tag_id) ? _gamma : n_oj;
		// 	oracle_p = n_oj / (n_o + _gamma);
		// }
		// double oracle_p = numerator / (n_o + _gamma);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	// P(y_t|s_t)
	double compute_Pword_tag(int word_id, int tag_id){
		double m_i = sum_word_count_for_tag(tag_id);
		double m_iq = get_tag_word_count(tag_id, word_id);
		// if(is_word_new(word_id) == false){
		// 	assert(m_iq > 0);
		// }
		double empirical_p = m_iq / (m_i + _beta_emission);
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);	// 親の分布から生成される確率. 親からword_idが生成される確率とは別物.
		double m_o = sum_oracle_words_count();
		double m_oq = get_oracle_count_for_word(word_id);
		double W = get_num_words();
		double g0 = 1.0 / (W + 1);
		double oracle_p = (m_oq + _gamma_emission) / (m_o + _gamma_emission);
		// if(is_word_new(word_id)){
		// 	oracle_p = _gamma_emission / (m_o + _gamma_emission);
		// }else{
		// 	oracle_p = m_oq / (m_o + _gamma_emission);
		// }
		// double numerator = is_word_new(word_id) ? _gamma_emission : m_oq;
		// double oracle_p = numerator / (m_o + _gamma_emission);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	double compute_gamma_distribution(double v, double a, double b){
		return pow(b, a) / tgamma(a) * pow(v, a - 1) * exp(-b * v);
	}
	double compute_log_posterior_alpha_and_beta(double alpha, double beta){
		assert(beta > 0);
		double a_alpha = 1;
		double b_alpha = 1;
		double gamma_dist_alpha = compute_gamma_distribution(_alpha, a_alpha, b_alpha);
		double a_beta = 1;
		double b_beta = 1;
		double gamma_dist_beta = compute_gamma_distribution(_beta, a_beta, b_beta);
		double log_sum = 0;
		for(int tag = 0;tag < _tag_unigram_count.size();tag++){
			if(is_tag_new(tag)){
				continue;
			}
			double Ki = sum_bigram_destination(tag);
			double log_first_term = (Ki - 1) * log(beta) + lgamma(alpha + beta) - lgamma(alpha);
			double n_ii = get_bigram_tag_count(tag, tag);
			double n_i = sum_bigram_destination(tag);
			double log_second_term = lgamma(n_ii + alpha) - lgamma(n_i + alpha + beta);
			log_sum += log_first_term + log_second_term;
		}
		return log(gamma_dist_alpha) + log(gamma_dist_beta) + log_sum;
	}
	double compute_log_posterior_gamma(double gamma){
		double a = 1;
		double b = 1;
		double gamma_dist = compute_gamma_distribution(gamma, a, b);
		double K = get_num_tags();
		double To = get_num_times_oracle_tag_used();
		double log_term = K * log(gamma) + lgamma(gamma) - lgamma(To + gamma);
		return log(gamma_dist) + log_term;
	}
	double compute_log_posterior_gamma_emission(double gamma){
		double a = 1;
		double b = 1;
		double gamma_dist = compute_gamma_distribution(gamma, a, b);
		double K = _num_words;
		double To = get_num_times_oracle_word_used();
		double log_term = K * log(gamma) + lgamma(gamma) - lgamma(To + gamma);
		// cout << gamma << endl;
		// cout << gamma_dist << endl;
		// cout << K << endl;
		// cout << To << endl;
		// cout << log(gamma) << endl;
		// cout << lgamma(gamma) << endl;
		// cout << lgamma(To + gamma) << endl;
		// cout << log_term << endl;
		// cout << endl;
		return log(gamma_dist) + log_term;
	}
	double compute_log_Pdata(vector<Word*> &line){
		double p = 0;
		for(int pos = 1;pos < line.size();pos++){
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			double log_p_tag = log(compute_Ptag_context(ti, ti_1));
			double log_p_word = log(compute_Pword_tag(wi, ti));
			p += log_p_tag + log_p_word;
		}
		return p;
	}
	void sample_alpha_and_beta(double lr = 1e-2){
		double eps = 1e-6;
		{
			double log_p_old = compute_log_posterior_alpha_and_beta(_alpha, _beta);
			double log_p_new = compute_log_posterior_alpha_and_beta(_alpha + eps, _beta);
			double grad = (log_p_new - log_p_old) / eps;
			_alpha -= grad * lr;
			_alpha = std::max(eps, _alpha);
		}
		{
			double log_p_old = compute_log_posterior_alpha_and_beta(_alpha, _beta);
			double log_p_new = compute_log_posterior_alpha_and_beta(_alpha, _beta + eps);
			double grad = (log_p_new - log_p_old) / eps;
			_beta -= grad * lr;
			_beta = std::max(eps, _beta);
		}
	}
	void sample_beta_emission(double lr = 1e-2){

	}
	void sample_gamma(double lr = 1e-2){
		double eps = 1e-6;
		double log_p_old = compute_log_posterior_gamma(_gamma);
		double log_p_new = compute_log_posterior_gamma(_gamma + eps);
		double grad = (log_p_new - log_p_old) / eps;
		_gamma -= grad * lr;
		_gamma = std::max(eps, _gamma);
	}
	void sample_gamma_emission(double lr = 1e-2){
		double eps = 1e-6;
		double log_p_old = compute_log_posterior_gamma_emission(_gamma_emission);
		double log_p_new = compute_log_posterior_gamma_emission(_gamma_emission + eps);
		double grad = (log_p_new - log_p_old) / eps;
		_gamma_emission -= grad * lr;
		_gamma_emission = std::max(eps, _gamma_emission);
	}
	// t_{i-1} -> t_i -> t_{i+1}
	int sample_new_tag(int ti_1, int ti1, int wi){
		// ギブスサンプリング
		vector<double> sampling_table;
		double sum = 0;
		bool new_tag_included = false;
		for(int tag = 0;tag < _tag_unigram_count.size();tag++){
			if(is_tag_new(tag)){
				if(new_tag_included){
					sampling_table.push_back(0);
					continue;
				}
				new_tag_included = true;
			}
			double p_emission = compute_Pword_tag(wi, tag);
			double p_generation = compute_Ptag_context(tag, ti_1);
			int correcting_count_for_bigram = (ti_1 == tag == ti1) ? 1 : 0;
			int correcting_count_for_destination = (ti_1 == tag) ? 1 : 0;
			double p_likelihood = compute_Ptag_context(ti1, tag, correcting_count_for_bigram, correcting_count_for_destination);
			double p_conditional = p_emission * p_generation * p_likelihood;
				// if(p_conditional == 0){
				// 	cout << "ti_1: " << ti_1 << endl;
				// 	cout << "tag: " << tag << endl;
				// 	cout << "ti1: " << ti1 << endl;
				// 	cout << p_emission << ",";
				// 	cout << p_generation << ",";
				// 	cout << p_likelihood << ",";
				// 	cout << p_conditional << ",";
				// 	cout << endl;
				// 	cout << get_oracle_count_for_tag(ti1) << endl;
				// 	cout << get_bigram_tag_count(tag, ti1) << endl;
				// 	cout << new_tag_included << endl;
				// 	exit(0);
				// }
			p_conditional = pow(p_conditional, 1.0 / _temperature);
			sampling_table.push_back(p_conditional);
			sum += p_conditional;
		}
		assert(sampling_table.size() == _tag_unigram_count.size());
		int new_tag = _tag_unigram_count.size();
		if(new_tag_included == false){
			double p_emission = compute_Pword_tag(wi, new_tag);
			double p_generation = compute_Ptag_context(new_tag, ti_1);
			double p_likelihood = compute_Ptag_context(ti1, new_tag);
			double p_conditional = p_emission * p_generation * p_likelihood;
				// cout << "tag: " << new_tag << endl;
				// cout << p_emission << ",";
				// cout << p_generation << ",";
				// cout << p_likelihood << ",";
				// cout << p_conditional << ",";
				// cout << endl;
			p_conditional = pow(p_conditional, 1.0 / _temperature);
			sampling_table.push_back(p_conditional);
			sum += p_conditional;
		}
		assert(sum > 0);
		double normalizer = 1.0 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int tag = 0;tag < _tag_unigram_count.size();tag++){
			sum += sampling_table[tag] * normalizer;
			if(bernoulli <= sum){
				return tag;
			}
		}
			// cout << sum << endl;
			// cout << "new_tag generated" << endl;
		return new_tag;
	}
	void perform_gibbs_sampling_with_line(vector<Word*> &line){
			// c_printf("[*]%s\n", "perform_gibbs_sampling_with_line");
			// for(int pos = 0;pos < line.size();pos++){
			// 	cout << line[pos]->tag_id << " -> ";
			// }
			// cout << endl;

		for(int pos = 1;pos < line.size() - 1;pos++){
				// c_printf("[*]%s\n", (boost::format("pos = %d") % pos).str().c_str());
			int ti_1 = line[pos - 1]->tag_id;
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			int ti1 = line[pos + 1]->tag_id;

			// 現在のtiをモデルから除去
			decrement_tag_bigram_count(ti_1, ti);
				// dump_oracle_tags();
				// dump_bigram_table();
			decrement_tag_bigram_count(ti, ti1);
				// dump_oracle_tags();
				// dump_bigram_table();
			decrement_tag_unigram_count(ti);
			decrement_tag_word_count(ti, wi);
			increment_tag_bigram_count(ti_1, ti1);

			int new_tag = sample_new_tag(ti_1, ti1, wi);
			// モデルに追加
			increment_tag_word_count(new_tag, wi);
			increment_tag_bigram_count(ti_1, new_tag);
				// dump_oracle_tags();
				// dump_bigram_table();
			increment_tag_bigram_count(new_tag, ti1);
				// dump_oracle_tags();
				// dump_bigram_table();
			increment_tag_unigram_count(new_tag);
			decrement_tag_bigram_count(ti_1, ti1);
			line[pos]->tag_id = new_tag;
		}
	}
	void dump_tags(){
		for(int tag = 0;tag < _tag_unigram_count.size();tag++){
			cout << _tag_unigram_count[tag] << ", ";
		}
		cout << endl;
	}
	void dump_oracle_tags(){
		c_printf("[*]%s\n", "dump_oracle_tags");
		for(const auto &elem: _oracle_tag_counts){
			cout << elem.first << ": " << elem.second << endl;
		}
	}
	void dump_oracle_words(){
		c_printf("[*]%s\n", "dump_oracle_words");
		for(const auto &elem: _oracle_word_counts){
			cout << elem.first << ": " << elem.second << endl;
		}
	}
	void dump_bigram_table(){
		c_printf("[*]%s\n", "dump_bigram_table");
		for(const auto &contexts: _bigram_tag_table){
			cout << contexts.first << ":" << endl;
			for(const auto &tables: contexts.second){
				cout << "	" << tables.first << ":" << endl;
				cout << "		";
				Table* table = tables.second;
				for(int i = 0;i < table->_arrangement.size();i++){
					cout << table->_arrangement[i] << ", ";
				}
				cout << endl;
			}
		}
	}
	void dump_hyperparameters(){
		c_printf("[*]%s\n", "dump_hyperparameters");
		cout << "alpha <- " << _alpha << endl;
		cout << "beta <- " << _beta << endl;
		cout << "beta_e <- " << _beta_emission << endl;
		cout << "gamma <- " << _gamma << endl;
		cout << "gamma_e <- " << _gamma_emission << endl;
	}
	void check_oracle_tag_count(){
		unordered_map<int, int> counts;
		for(const auto &contexts: _bigram_tag_table){
			for(const auto &tags: contexts.second){
				vector<int> &table = tags.second->_arrangement;
				int num_tables = table.size();
				counts[tags.first] += num_tables;
			}
		}
		for(const auto &elem: counts){
			int tag_id = elem.first;
			int num_tables = elem.second;
			int count = get_oracle_count_for_tag(tag_id);
			assert(num_tables == count);
		}
	}
	void check_oracle_word_count(){
		unordered_map<int, int> counts;
		for(const auto &tags: _tag_word_table){
			for(const auto &words: tags.second){
				vector<int> &table = words.second->_arrangement;
				int num_tables = table.size();
				counts[words.first] += num_tables;
			}
		}
		for(const auto &elem: counts){
			int word_id = elem.first;
			int num_tables = elem.second;
			int count = get_oracle_count_for_word(word_id);
			assert(num_tables == count);
		}
	}
	void check_sum_bigram_destination(){
		for(const auto &tag: _sum_bigram_destination){
			int tag_id = tag.first;
			int count = tag.second;
			int sum = 0;
			for(const auto &table: _bigram_tag_table[tag_id]){
				sum += table.second->_num_customers;
			}
			assert(count == sum);
		}
	}
	void check_tag_count(){
		int num_non_zero = 0;
		for(int i = 0;i < _tag_unigram_count.size();i++){
			if(_tag_unigram_count[i] > 0){
				num_non_zero += 1;
			}
		}
		assert(num_non_zero == _oracle_tag_counts.size());
	}
	void check_sum_word_customers(){
		int num_customers = 0;
		for(const auto &tag: _tag_word_table){
			for(const auto &word: tag.second){
				Table* table = word.second;
				num_customers += table->_num_customers;
			}
		}
		assert(num_customers == _num_words);
	}
	void check_sum_tag_customers(){
		int num_customers_in_bigram = 0;
		for(const auto &bigram: _bigram_tag_table){
			for(const auto &unigram: bigram.second){
				Table* table = unigram.second;
				num_customers_in_bigram += table->_num_customers;
			}
		}
		int num_customers_in_unigram = 0;
		for(auto itr = _tag_unigram_count.begin();itr != _tag_unigram_count.end();itr++){
			num_customers_in_unigram += *itr;
		}
		assert(num_customers_in_bigram == num_customers_in_unigram);
	}
	bool save(string dir = "out"){
		bool success = false;
		ofstream ofs(dir + "/ihmm.model");
		if(ofs.good()){
			boost::archive::binary_oarchive oarchive(ofs);
			oarchive << static_cast<const InfiniteHMM&>(*this);
			success = true;
		}
		ofs.close();
		return success;
	}
	bool load(string dir = "out"){
		bool success = false;
		ifstream ifs(dir + "/ihmm.model");
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> *this;
			success = true;
		}
		ifs.close();
		return success;
	}

};

#endif