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

#define BOP 0
#define EOP 0

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
	int _prev_tag_unigram_count_size;
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
	int _max_sequence_length;
	double _temperature;
	unordered_map<int, int> _sum_bigram_destination;
	unordered_map<int, int> _sum_word_count_for_tag;
	double* _gibbs_sampling_table;
	double* _beam_sampling_table_u;
	double** _beam_sampling_table_s;
	InfiniteHMM(int initial_num_tags){
		_alpha = 0.1;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_initial_num_tags = initial_num_tags;
		_sum_oracle_words_count = 0;
		_sum_oracle_tags_count = 0;
		_num_words = 0;
		_temperature = 1;
		_max_sequence_length = 0;
		_gibbs_sampling_table = NULL;
		_beam_sampling_table_u = NULL;
		_beam_sampling_table_s = NULL;
	}
	void initialize(vector<vector<Word*>> &dataset){
		// サンプリングテーブル
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			if(line.size() > _max_sequence_length){
				_max_sequence_length = line.size();
			}
		}
		for(int tag = 0;tag < _initial_num_tags;tag++){
			_tag_unigram_count.push_back(0);
		}
		_prev_tag_unigram_count_size = _tag_unigram_count.size();
		_gibbs_sampling_table = (double*)malloc(_prev_tag_unigram_count_size * sizeof(double));
		assert(_max_sequence_length > 0);
		_beam_sampling_table_u = (double*)malloc((_max_sequence_length + 1) * sizeof(double));
		_beam_sampling_table_s = (double**)malloc(_max_sequence_length * sizeof(double));
		for(int pos = 0;pos < _max_sequence_length;pos++){
			_beam_sampling_table_s[pos] = (double*)malloc((_initial_num_tags + 1) * sizeof(double));
		}
		// nグラムカウントテーブル
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
			// Word* word = line[0];
			// word->tag_id = BOP;
			// word_set.insert(word->word_id);
			// increment_tag_word_count(line[0]);
			// increment_tag_unigram_count(line[0]->tag_id);
			
			int ti_1 = BOP;
			for(int pos = 0;pos < line.size();pos++){	// 2-gramなので2番目から.
				Word* word = line[pos];
				int ti = Sampler::uniform_int(EOP + 1, _initial_num_tags - 1);
				increment_tag_bigram_count(ti_1, ti);
				increment_tag_unigram_count(ti);
				int wi = word->word_id;
				increment_tag_word_count(ti, wi);
				word_set.insert(wi);
				num_words += 1;
				word->tag_id = ti;
				ti_1 = ti;
			}
			increment_tag_bigram_count(ti_1, EOP);
		}
		c_printf("[*]%s\n", (boost::format("単語数: %d - 単語種: %d - 行数: %d") % num_words % word_set.size() % dataset.size()).str().c_str());
		_num_words += num_words;
	}
	void increment_tag_unigram_count(int tag_id){
		assert(tag_id != BOP);
		assert(tag_id != EOP);
		while(tag_id >= _tag_unigram_count.size()){
			_tag_unigram_count.push_back(0);
		}
		_tag_unigram_count[tag_id] += 1;
		if(_tag_unigram_count.size() != _prev_tag_unigram_count_size){
			assert(_gibbs_sampling_table != NULL);
			free(_gibbs_sampling_table);
			assert(_beam_sampling_table_s != NULL);
			for(int pos = 0;pos < _max_sequence_length;pos++){
				assert(_beam_sampling_table_s[pos] != NULL);
				free(_beam_sampling_table_s[pos]);
			}
			// 更新
			_prev_tag_unigram_count_size = _tag_unigram_count.size();
			// re-alloc
			_gibbs_sampling_table = (double*)malloc(_prev_tag_unigram_count_size * sizeof(double));
			for(int pos = 0;pos < _max_sequence_length;pos++){
				// ここだけは既存タグ数+1
				_beam_sampling_table_s[pos] = (double*)malloc((_prev_tag_unigram_count_size + 1) * sizeof(double));
			}
		}
	}
	void decrement_tag_unigram_count(int tag_id){
		assert(tag_id != BOP);
		assert(tag_id != EOP);
		assert(tag_id < _tag_unigram_count.size());
		_tag_unigram_count[tag_id] -= 1;
		assert(_tag_unigram_count[tag_id] >= 0);
		int num_pop = 0;
		for(int i = _tag_unigram_count.size() - 1;i >= 0;i--){
			if(_tag_unigram_count[i]){
				break;
			}
			num_pop += 1;
		}
		for(int n = 0;n < num_pop;n++){
			_tag_unigram_count.pop_back();
		}
		if(_tag_unigram_count.size() != _prev_tag_unigram_count_size){
			assert(_gibbs_sampling_table != NULL);
			free(_gibbs_sampling_table);
			assert(_beam_sampling_table_s != NULL);
			for(int pos = 0;pos < _max_sequence_length;pos++){
				assert(_beam_sampling_table_s[pos] != NULL);
				free(_beam_sampling_table_s[pos]);
			}
			// 更新
			_prev_tag_unigram_count_size = _tag_unigram_count.size();
			// re-alloc
			_gibbs_sampling_table = (double*)malloc(_prev_tag_unigram_count_size * sizeof(double));
			for(int pos = 0;pos < _max_sequence_length;pos++){
				// ここだけは既存タグ数+1
				_beam_sampling_table_s[pos] = (double*)malloc((_prev_tag_unigram_count_size + 1) * sizeof(double));
			}
		}
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
		assert(tag_id != EOP);
		if(tag_id >= _tag_unigram_count.size()){
			return true;
		}
		if(_tag_unigram_count[tag_id] == 0){
			return true;
		}
		return false;
	}
	int get_new_tag_id(){
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
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
	}
	int sum_word_count_for_tag(int tag_id){
		return _sum_word_count_for_tag[tag_id];
	}
	int sum_bigram_destination(int tag_id){
		return _sum_bigram_destination[tag_id];
	}
	int sum_oracle_tags_count(){
		return _sum_oracle_tags_count;
	}
	// P(s_{t+1}|s_t)
	double compute_Ptag_context(int tag_id, int context_tag_id, int correcting_count_for_bigram = 0, int correcting_count_for_destination = 0){
		double n_i = sum_bigram_destination(context_tag_id);
		double n_ij = get_bigram_tag_count(context_tag_id, tag_id);
		double alpha = (tag_id == context_tag_id) ? _alpha : 0;
		double empirical_p = (n_ij + alpha) / (n_i + _beta + _alpha);
		double coeff_oracle_p = _beta / (n_i + _beta + _alpha);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		double n_o = sum_oracle_tags_count();
		double n_oj = get_oracle_count_for_tag(tag_id);
		double T = get_num_tags();
		double g0 = 1.0 / (T + 1.0);
		double oracle_p = (n_oj + _gamma * g0) / (n_o + _gamma);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	// P(y_t|s_t)
	double compute_Pword_tag(int word_id, int tag_id){
		double m_i = sum_word_count_for_tag(tag_id);
		double m_iq = get_tag_word_count(tag_id, word_id);
		double empirical_p = m_iq / (m_i + _beta_emission);
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);	// 親の分布から生成される確率. 親からword_idが生成される確率とは別物.
		double m_o = sum_oracle_words_count();
		double m_oq = get_oracle_count_for_word(word_id);
		double W = get_num_words();
		double g0 = 1.0 / (W + 1);
		// cout << "tag = " << tag_id <<  ", m_i = " << m_i << ", m_iq = " << m_iq << endl;
		// cout << "m_o = " << m_o << ", m_oq = " << m_oq << endl;
		double oracle_p = (m_oq + _gamma_emission * g0) / (m_o + _gamma_emission);
		// cout << "empirical_p = " << empirical_p << ", coeff_oracle_p = " << coeff_oracle_p << ", oracle_p = " << oracle_p << endl;
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	double compute_gamma_distribution(double v, double a, double b){
		return pow(b, a) / tgamma(a) * pow(v, a - 1) * exp(-b * v);
	}
	double compute_log_Pdata(vector<Word*> &line){
		double p = 0;
		int ti_1 = BOP;
		for(int pos = 0;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			double log_p_tag = log(compute_Ptag_context(ti, ti_1));
			double log_p_word = log(compute_Pword_tag(wi, ti));
			p += log_p_tag + log_p_word;
			ti_1 = ti;
		}
		return p;
	}
	// t_{i-1} -> t_i -> t_{i+1}
	int gibbs_sample_new_tag(int ti_1, int ti1, int wi){
		// ギブスサンプリング
		double sum = 0;
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
			if(is_tag_new(tag)){
				_gibbs_sampling_table[tag] = 0;
				continue;
			}
			double p_emission = compute_Pword_tag(wi, tag);
			double p_generation = compute_Ptag_context(tag, ti_1);
			int correcting_count_for_bigram = (ti_1 == tag == ti1) ? 1 : 0;
			int correcting_count_for_destination = (ti_1 == tag) ? 1 : 0;
			double p_likelihood = compute_Ptag_context(ti1, tag, correcting_count_for_bigram, correcting_count_for_destination);
			double p_conditional = p_emission * p_generation * p_likelihood;
			// p_conditional = pow(p_conditional, 1.0 / _temperature);
			_gibbs_sampling_table[tag] = p_conditional;
			sum += p_conditional;
		}
		int new_tag = get_new_tag_id();
		double p_emission = compute_Pword_tag(wi, new_tag);
		double p_generation = compute_Ptag_context(new_tag, ti_1);
		double p_likelihood = compute_Ptag_context(ti1, new_tag);
		double p_conditional = p_emission * p_generation * p_likelihood;
		// p_conditional = pow(p_conditional, 1.0 / _temperature);
		sum += p_conditional;
		// new_tag > _tag_unigram_count.size()ならサンプリングテーブルに入れなくてもよい.
		if(new_tag < _tag_unigram_count.size()){
			_gibbs_sampling_table[new_tag] = p_conditional;
		}
		assert(sum > 0);
		double normalizer = 1.0 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
			sum += _gibbs_sampling_table[tag] * normalizer;
			if(bernoulli <= sum){
				return tag;
			}
		}
		return new_tag;
	}
	int argmax_Ptag_context_word(int context_tag_id, int word_id){
		double max_p = 0;
		double max_tag = 0;
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
			if(is_tag_new(tag)){
				continue;
			}
			double Ptag = compute_Ptag_context(tag, context_tag_id);
			double Pword = compute_Pword_tag(word_id, tag);
			double p = Ptag * Pword;
			if(p > max_p){
				max_p = p;
				max_tag = tag;
			}
		}
		return max_tag;
	}
	void perform_gibbs_sampling_with_line(vector<Word*> &line){
		if(line.size() < 2){
			return;
		}
		int pos = 0;
		int ti_1 = BOP;
		for(;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			int ti1 = (pos == line.size() - 1) ? EOP : line[pos + 1]->tag_id;

			// 現在のtiをモデルから除去
			decrement_tag_bigram_count(ti_1, ti);
			decrement_tag_bigram_count(ti, ti1);
			decrement_tag_unigram_count(ti);
			decrement_tag_word_count(ti, wi);
			increment_tag_bigram_count(ti_1, ti1);
			// 新しい状態をサンプリング
			int new_tag = gibbs_sample_new_tag(ti_1, ti1, wi);
			// モデルに追加
			increment_tag_word_count(new_tag, wi);
			increment_tag_bigram_count(ti_1, new_tag);
			increment_tag_bigram_count(new_tag, ti1);
			increment_tag_unigram_count(new_tag);
			decrement_tag_bigram_count(ti_1, ti1);
			line[pos]->tag_id = new_tag;
			ti_1 = new_tag;
		}
	}
	void perform_beam_sampling_with_line(vector<Word*> &line){
		if(line.size() < 1){
			return;
		}
		// 品詞をモデルから除去
		int ti_1 = BOP;
		for(int pos = 0;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			int wi = line[pos]->word_id;
			decrement_tag_bigram_count(ti_1, ti);
			decrement_tag_unigram_count(ti);
			decrement_tag_word_count(ti, wi);
			ti_1 = ti;
			// cout << wi << " -> ";
			// cout << "pos = " << pos << endl;
		}
		// cout << endl;
		decrement_tag_bigram_count(ti_1, EOP);

		// dump_tags();

		// for(int ti_1 = BOP;ti_1 < _tag_unigram_count.size();ti_1++){
		// 	for(int ti = BOP;ti < _tag_unigram_count.size();ti++){
		// 		cout << (boost::format("p(%d|%d) = %f") % ti % ti_1 % compute_Ptag_context(ti, ti_1)).str() << endl;
		// 	}
		// }

		// for(int wi = 0;wi < _num_words;wi++){
		// 	for(int ti = BOP;ti < _tag_unigram_count.size();ti++){
		// 		cout << (boost::format("p(wi=%d|%d) = %f") % wi % ti % compute_Pword_tag(wi, ti)).str() << endl;
		// 	}
		// }

		// cout << "starting u" << endl;
		// uのサンプリング
		ti_1 = BOP;
		for(int pos = 0;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			double p = compute_Ptag_context(ti, ti_1);
			_beam_sampling_table_u[pos] = Sampler::uniform(0, p);
			// cout << (boost::format("u[%d] <- %f; p = %f; %d -> %d") % pos % _beam_sampling_table_u[pos] % p % ti_1 % ti).str() << endl;
			ti_1 = ti;
		}
		double p = compute_Ptag_context(EOP, ti_1);
		_beam_sampling_table_u[line.size()] = Sampler::uniform(0, p);
		// cout << (boost::format("u[%d] <- %f; p = %f; %d -> %d") % line.size() % _beam_sampling_table_u[line.size()] % p % ti_1 % EOP).str() << endl;
		// sのサンプリング
		// cout << "starting s" << endl;
		int new_tag = get_new_tag_id();
		// cout << "new_tag: " << new_tag << endl;
		//// forwardパス
		////// pos == 1
		int wi = line[0]->word_id;
		double ui = _beam_sampling_table_u[0];
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
			if(is_tag_new(tag)){
				_beam_sampling_table_s[0][tag] = 0;
				continue;
			}
			double Pti_ti_1 = compute_Ptag_context(tag, BOP);
			double Pwi_ti = compute_Pword_tag(wi, tag);
			_beam_sampling_table_s[0][tag] = (Pti_ti_1 >= ui) ? Pwi_ti : 0; 

			// if(line.size() == 200){
				// cout << (boost::format("s[%d][%d] <- %e; %f > %f; p(wi=%d|%d) = %f") % 0 % tag % _beam_sampling_table_s[0][tag] % Pti_ti_1 % ui % wi % tag % Pwi_ti).str() << endl;
			// }
		}
		double Pti_ti_1 = compute_Ptag_context(new_tag, BOP);
		double Pwi_ti = compute_Pword_tag(wi, new_tag);
		_beam_sampling_table_s[0][new_tag] = (Pti_ti_1 > ui) ? Pwi_ti : 0; 

		// if(line.size() == 200){
			// cout << (boost::format("s[%d][%d] <- %e; %f > %f; p(wi=%d|%d) = %f") % 0 % new_tag % _beam_sampling_table_s[0][new_tag] % Pti_ti_1 % ui % wi % new_tag % Pwi_ti).str() << endl;
		// }
		////// pos > 1
		for(int pos = 1;pos < line.size();pos++){
			int wi = line[pos]->word_id;
			double ui = _beam_sampling_table_u[pos];
			double sum_over_tag = 0;
			// 新しい品詞以外
			for(int ti = EOP + 1;ti < _tag_unigram_count.size();ti++){
				if(is_tag_new(ti)){
					_beam_sampling_table_s[pos][ti] = 0;
					continue;
				}
				double Pwi_ti = compute_Pword_tag(wi, ti);
				// 周辺化
				double sum = 0;
				for(int ti_1 = EOP + 1;ti_1 < _tag_unigram_count.size();ti_1++){
					if(is_tag_new(ti_1)){
						continue;
					}
					if(_beam_sampling_table_s[pos - 1][ti_1] > 0){
						double Pti_ti_1 = compute_Ptag_context(ti, ti_1);
						if(Pti_ti_1 >= ui){
							sum += _beam_sampling_table_s[pos - 1][ti_1];

							// if(line.size() == 200){
								// cout << (boost::format("sum += %e; %f > %f") % _beam_sampling_table_s[pos - 1][ti_1] % Pti_ti_1 % ui).str() << endl;
							// }
						}
					}
				}
				if(_beam_sampling_table_s[pos - 1][new_tag] > 0){
					double Pti_ti_1 = compute_Ptag_context(ti, new_tag);
					if(Pti_ti_1 >= ui){
						sum += _beam_sampling_table_s[pos - 1][new_tag];

						// if(line.size() == 200){
							// cout << (boost::format("sum += %e; %f > %f") % _beam_sampling_table_s[pos - 1][new_tag] % Pti_ti_1 % ui).str() << endl;
						// }
					}
				}
				_beam_sampling_table_s[pos][ti] = Pwi_ti * sum;
				sum_over_tag += Pwi_ti * sum;

				// if(line.size() == 200){
					// cout << (boost::format("s[%d][%d] <- %e; p = %e * %e") % pos % ti % (Pwi_ti * sum) % Pwi_ti % sum).str() << endl;
				// }
			}
			// 新しい品詞
			{
				double Pwi_ti = compute_Pword_tag(wi, new_tag);
				// 周辺化
				double sum = 0;
				for(int ti_1 = EOP + 1;ti_1 < _tag_unigram_count.size();ti_1++){
					if(is_tag_new(ti_1)){
						continue;
					}
					if(_beam_sampling_table_s[pos - 1][ti_1] > 0){
						double Pti_ti_1 = compute_Ptag_context(new_tag, ti_1);
						if(Pti_ti_1 >= ui){
							sum += _beam_sampling_table_s[pos - 1][ti_1];

							// if(line.size() == 200){
								// cout << (boost::format("sum += %e; %f > %f") % _beam_sampling_table_s[pos - 1][ti_1] % Pti_ti_1 % ui).str() << endl;
							// }
						}
					}
				}
				if(_beam_sampling_table_s[pos - 1][new_tag] > 0){
					double Pti_ti_1 = compute_Ptag_context(new_tag, new_tag);
					if(Pti_ti_1 >= ui){
						sum += _beam_sampling_table_s[pos - 1][new_tag];

						// if(line.size() == 200){
							// cout << (boost::format("sum += %e; %f > %f") % _beam_sampling_table_s[pos - 1][new_tag] % Pti_ti_1 % ui).str() << endl;
						// }
					}
				}
				_beam_sampling_table_s[pos][new_tag] = Pwi_ti * sum;
				sum_over_tag += Pwi_ti * sum;

				// if(line.size() == 200){
					// cout << (boost::format("s[%d][%d] <- %e; p = %e * %e") % pos % new_tag % (Pwi_ti * sum) % Pwi_ti % sum).str() << endl;
				// }
			}
			assert(sum_over_tag > 0);
			for(int ti = EOP + 1;ti < _tag_unigram_count.size();ti++){
				if(is_tag_new(ti)){
					continue;
				}
				_beam_sampling_table_s[pos][ti] /= sum_over_tag;
			}
			_beam_sampling_table_s[pos][new_tag] /= sum_over_tag;
		}
		//// backwardパス
		// cout << "backward" << endl;
		int sampled_tag = EOP;
		for(int pos = line.size() - 1;pos >= 0;pos--){
			double sum = 0;
			double ui1 = _beam_sampling_table_u[pos + 1];
			for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
				if(is_tag_new(tag)){
					_gibbs_sampling_table[tag] = 0;
					continue;
				}
				double Pti_yi_ui = _beam_sampling_table_s[pos][tag];
				double Pti1_ti = compute_Ptag_context(sampled_tag, tag);
				if(Pti1_ti < ui1){
					Pti1_ti = 0;
				}else{
					Pti1_ti = 1;
				}
				// if(line.size() == 200){
				// }
				double Pti_ti1_yi_ui = Pti1_ti * Pti_yi_ui;
				sum += Pti_ti1_yi_ui;
				// cout << (boost::format("sum += %e; sum = %e; s[%d][%d] = %e;") % p % sum % pos % tag % p).str() << endl;
				_gibbs_sampling_table[tag] = Pti_ti1_yi_ui;
			}
			double Pti_yi_ui = _beam_sampling_table_s[pos][new_tag];
			double Pti1_ti = compute_Ptag_context(sampled_tag, new_tag);
			if(Pti1_ti < ui1){
				Pti1_ti = 0;
			}else{
				Pti1_ti = 1;
			}
			double Pti_ti1_yi_ui = Pti1_ti * Pti_yi_ui;
			// cout << (boost::format("sum += %e; sum = %e; s[%d][%d] = %e;") % p % sum % pos % new_tag % p).str() << endl;
			sum += Pti_ti1_yi_ui;
			// new_tag > _tag_unigram_count.size()ならサンプリングテーブルに入れない.
			if(new_tag < _tag_unigram_count.size()){
				_gibbs_sampling_table[new_tag] = Pti_ti1_yi_ui;
			}
			// if(sum <= 0){
			// 	cout << line.size() << endl;
			// }
			assert(sum > 0);
			double normalizer = 1.0 / sum;
			double bernoulli = Sampler::uniform(0, 1);
			sum = 0;
			sampled_tag = -1;
			for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
				sum += _gibbs_sampling_table[tag] * normalizer;
				if(bernoulli <= sum){
					sampled_tag = tag;
					break;
				}
			}
			sampled_tag = (sampled_tag == -1) ? new_tag : sampled_tag;
			// cout << "sampled: " << sampled_tag << endl;
			line[pos]->tag_id = sampled_tag;
		}

		// サンプリングした品詞をモデルに追加
		ti_1 = BOP;
		for(int pos = 0;pos < line.size();pos++){
			int ti = line[pos]->tag_id;
			increment_tag_bigram_count(ti_1, ti);
			int wi = line[pos]->word_id;
			increment_tag_unigram_count(ti);
			increment_tag_word_count(ti, wi);
			ti_1 = ti;
		}
		increment_tag_bigram_count(ti_1, EOP);
		// exit(0);
		// exit(0);
	}
	void dump_tags(){
		for(int tag = EOP + 1;tag < _tag_unigram_count.size();tag++){
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