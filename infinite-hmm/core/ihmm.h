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

class Table{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _arrangement;
		archive & _num_customers;
	}
public:
	vector<int> _arrangement;
	int _num_customers;
	int _token_id;
	Table(int token_id){
		_num_customers = 0;
		_token_id = token_id;
	}
	bool is_empty(){
		return _arrangement.size() == 0;
	}
	void add_customer(double self_transition_alpha, double concentration_parameter, bool &new_table_generated){
		_num_customers += 1;
		if(_arrangement.size() == 0){
			_arrangement.push_back(1);
			new_table_generated = true;
			return;
		}
		new_table_generated = false;
		double sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0) + self_transition_alpha + concentration_parameter;
		double normalizer = 1 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int i = 0;i < _arrangement.size();i++){
			sum += (_arrangement[i] + self_transition_alpha / _arrangement.size()) * normalizer;
			if(bernoulli < sum){
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
		archive & _tag_count;
		archive & _bigram_tag_table;
		archive & _tag_word_table;
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
				increment_tag_bigram_count(ti_1, ti);
				increment_tag_count(ti);

				_compute_Pword_tag(wi, ti, empirical_p, coeff_oracle_p);
				double normalizer = 1 / (empirical_p + coeff_oracle_p);
				double bernoulli = Sampler::uniform(0, 1);
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
		double alpha = (context_tag_id == tag_id) ? _alpha : 0;
		table->add_customer(alpha, _beta, new_table_generated);
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
		table->add_customer(0, _beta_emission, new_table_generated);
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
	int get_oracle_tag_count(int tag_id){
		auto itr = _oracle_tag_counts.find(tag_id);
		if(itr == _oracle_tag_counts.end()){
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
			double n_ij = get_bigram_tag_count(context_tag_id, tag_id);
			empirical_p = n_ij / (n_i + _beta);
		}
		coeff_oracle_p = _beta / (n_i + _beta);
	}
	// P(s_{t+1}|s_t)
	double compute_Ptag_context(int tag_id, int context_tag_id){
		double empirical_p = 0;		// 自分が生成したサンプルからなる経験確率.新しい品詞の場合は存在しないので0.
		double n_i = sum_bigram_destination(context_tag_id);
		if(is_tag_new(tag_id) == false){
			double n_ij = get_bigram_tag_count(context_tag_id, tag_id);
			empirical_p = n_ij / (n_i + _beta);
		}
		double coeff_oracle_p = _beta / (n_i + _beta);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		double n_o = sum_oracle_tags_count();
		double n_oj = get_oracle_tag_count(tag_id);
		double oracle_p = (n_oj + _gamma) / (n_o + _gamma);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	void _compute_Pword_tag(int word_id, int tag_id, double &empirical_p, double &coeff_oracle_p){
		empirical_p = 0;
		double m_i = sum_word_count_for_tag(tag_id);
		if(is_tag_new(tag_id) == false){
			double m_iq = get_tag_word_count(tag_id, word_id);
			empirical_p = m_iq / (m_i + _beta_emission);
		}
		coeff_oracle_p = _beta_emission / (m_i + _beta_emission);
	}
	// P(y_t|s_t)
	double compute_Pword_tag(int word_id, int tag_id){
		double empirical_p = 0;		// 自分が生成したサンプルからなる経験確率.新しい品詞の場合は存在しないので0.
		double m_i = sum_word_count_for_tag(tag_id);
		if(is_tag_new(tag_id) == false){
			double m_iq = get_tag_word_count(tag_id, word_id);
			empirical_p = m_iq / (m_i + _beta_emission);
		}
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);	// 親の分布から生成される確率. 親からword_idが生成される確率とは別物.
		double m_o = sum_oracle_words_count();
		double m_oq = _oracle_word_counts[word_id];
		double oracle_p = (m_oq + _gamma_emission) / (m_o + _gamma_emission);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	void remove_tag_from_bigram_count(int context_tag_id, int tag_id){
		// assert(is_tag_new(tag_id) == false);
		// double n_ij = get_bigram_tag_count(context_tag_id, tag_id);
		// double n_i = sum_bigram_destination(context_tag_id);
		// double empirical_p = n_ij / (n_i + _beta);
		// double coeff_oracle_p = _beta / (n_i + _beta);	// 親の分布から生成される確率. 親からtag_idが生成される確率とは別物.
		// double n_o = sum_oracle_tags_count();
		// double n_oj = _oracle_tag_counts[tag_id];
		// if(n_oj == 0){
		// 	decrement_tag_bigram_count(context_tag_id, tag_id);
		// 	return;
		// }
		// // double oracle_p = n_oj / (n_o + _gamma);
		// double normalizer = 1 / (empirical_p + coeff_oracle_p);
		// double bernoulli = Sampler::uniform(0, 1);
		// if(bernoulli < empirical_p * normalizer){
		// 	decrement_tag_bigram_count(context_tag_id, tag_id);
		// }else{
		// 	decrement_oracle_tag_count(tag_id);
		// }
	}
	void remove_word_from_model(int word_id, int tag_id){
		assert(is_tag_new(tag_id) == false);
		double m_iq = get_tag_word_count(tag_id, word_id);
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
			// cout << "sampling_new_tag" << endl;
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
				// cout << "tag: " << tag << endl;
				// cout << p_emission << ",";
				// cout << p_generation << ",";
				// cout << p_likelihood << ",";
				// cout << endl;
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
			// remove_word_from_model(wi, ti);	// 先に品詞-単語ペアから除去するとassrtで引っかからない
			decrement_tag_bigram_count(ti_1, ti);
				// dump_oracle_tag();
				// dump_bigram_table();
			decrement_tag_bigram_count(ti, ti1);
				// dump_oracle_tag();
				// dump_bigram_table();
			decrement_tag_count(ti);
			decrement_tag_word_count(ti, wi);

			int new_tag = sampling_new_tag(ti_1, ti, ti1, wi);
			// モデルに追加
			increment_tag_word_count(new_tag, wi);
			increment_tag_bigram_count(ti_1, new_tag);
				// dump_oracle_tag();
				// dump_bigram_table();
			increment_tag_bigram_count(new_tag, ti1);
				// dump_oracle_tag();
				// dump_bigram_table();
			increment_tag_count(new_tag);
			line[pos]->tag_id = new_tag;
		}
	}
	void dump_tags(){
		for(int tag = 0;tag < _tag_count.size();tag++){
			cout << _tag_count[tag] << ", ";
		}
		cout << endl;
	}
	void dump_oracle_tag(){
		c_printf("[*]%s\n", "dump_oracle_tag");
		for(const auto &elem: _oracle_tag_counts){
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
	bool load(string dir = "out"){
		return true;
	}
	bool save(string dir = "out"){
		return true;
	}

};

#endif