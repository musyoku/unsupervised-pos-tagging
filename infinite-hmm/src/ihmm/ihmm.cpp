#include <cassert>
#include <iostream>
#include <set>
#include "ihmm.h"

// num_tagsは全て特殊タグを除いた個数を表す
// 内部的には<s></s>を考慮するためnum_tags+1個のタグがある
// 例えばnum_tags=5のとき、ユーザーから見ればタグは[1, 2, 3, 4, 5]であり、内部的には[0, 1, 2, 3, 4, 5]となっている

namespace ihmm {
	InfiniteHMM::InfiniteHMM(){
		_initial_num_tags = 0;
		_num_words = 0;
		_alpha = 1;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_oracle_sum_n_over_j = 0;
		_oracle_sum_m_over_q = 0;
	}
	InfiniteHMM::InfiniteHMM(int initial_num_tags, int num_words): InfiniteHMM(){
		_initial_num_tags = initial_num_tags;
		_num_words = num_words;
		for(int tag = 0;tag <= _initial_num_tags;tag++){	// <s>と</s>のID=0
			_add_new_tag();
		}
		_prev_num_tags = _initial_num_tags;
		_oracle_m_q_counts = new int[num_words];
	}
	InfiniteHMM::~InfiniteHMM(){
		for(int tag = get_num_tags();tag > 0;tag--){
			_delete_tag(tag);
		}
		// <s></s>を削除
		std::cout << _m_iq_tables.size() << endl;
		Table** table_array = _m_iq_tables[0];
		delete[] _oracle_m_q_counts;
	}
	// <s></s>を除いたタグの数を返す
	int InfiniteHMM::get_num_tags(){
		return _oracle_n_j_counts.size() - 1;
	}
	int InfiniteHMM::get_num_words(){
		return _num_words;
	}
	void InfiniteHMM::initialize_with_training_corpus(std::vector<std::vector<Word*>> &dataset){
		// 最初は品詞をランダムに割り当てる
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &word_ids = dataset[data_index];
			int ti_1 = 0;
			for(int i = 0;i < word_ids.size();i++){
				Word* word = word_ids[i];
				int ti = sampler::uniform_int(1, _initial_num_tags);
				int wi = word->_id;
				_increment_tag_bigram_count(ti_1, ti);
				_increment_tag_word_count(ti, wi);
				word->_tag = ti;
				ti_1 = ti;
			}
			_increment_tag_bigram_count(ti_1, 0);	// </s>への遷移
		}
	}
	int InfiniteHMM::_add_new_tag(){
		int new_num_tags = get_num_tags() + 1;	// <s></s>は除く
		// bigramのカウント
		// 0: *, *				0: *, *
		// 1: *, *		->   	1: *, *
		// 						2: *, *, *
		std::vector<Table*> table_vec;
		table_vec.push_back(NULL);
		for(int tag = 1;tag <= new_num_tags;tag++){
			table_vec.push_back(new Table());
		}
		_n_ij_tables.push_back(table_vec);
		// 0: *, *				0: *, *, *
		// 1: *, *		->   	1: *, *, *
		// 2: *, *, *			2: *, *, *
		for(int context_tag = 1;context_tag < new_num_tags;context_tag++){
			std::vector<Table*> &table_vec = _n_ij_tables[context_tag];
			assert(table_vec.size() == new_num_tags);
			table_vec.push_back(new Table());
		}
		// tagと単語のペアの出現頻度
		Table** table_array = new Table*[_num_words];
		for(id word_id = 0;word_id < _num_words;word_id++){
			table_array[word_id] = new Table();
		}
		_m_iq_tables.push_back(table_array);
		// その他カウント
		_oracle_n_j_counts.push_back(0);
		_sum_m_i_over_q.push_back(0);
		_sum_n_i_over_j.push_back(0);

		assert(_n_ij_tables.size() == new_num_tags + 1);
		assert(_m_iq_tables.size() == new_num_tags + 1);
		assert(_oracle_n_j_counts.size() == new_num_tags + 1);
		assert(_sum_m_i_over_q.size() == new_num_tags + 1);
		assert(_sum_n_i_over_j.size() == new_num_tags + 1);
		return new_num_tags;
	}
	// タグが無くなったら詰める
	void InfiniteHMM::_delete_tag(int tag){
		assert(1 <= tag && tag <= get_num_tags());
		int new_num_tags = get_num_tags() - 1;	// <s></s>は除く
		std::vector<Table*> &tables_vec = _n_ij_tables[tag];
		for(Table* table: tables_vec){
			delete table;
		}
		Table** tables_array = _m_iq_tables[tag];
		for(id word_id = 0;word_id < _num_words;word_id++){
			delete tables_array[word_id];
		}
		delete[] tables_array;
		for(int t = tag;t <= new_num_tags;t++){
			_n_ij_tables[t] = _n_ij_tables[t + 1];
			_m_iq_tables[t] = _m_iq_tables[t + 1];
			_oracle_n_j_counts[t] = _oracle_n_j_counts[t + 1];
			_sum_m_i_over_q[t] = _sum_m_i_over_q[t + 1];
			_sum_n_i_over_j[t] = _sum_n_i_over_j[t + 1];
		}
		_n_ij_tables.pop_back();
		_m_iq_tables.pop_back();
		_oracle_n_j_counts.pop_back();
		_sum_m_i_over_q.pop_back();
		_sum_n_i_over_j.pop_back();
		assert(_n_ij_tables.size() == new_num_tags + 1);
		assert(_m_iq_tables.size() == new_num_tags + 1);
		assert(_oracle_n_j_counts.size() == new_num_tags + 1);
		assert(_sum_m_i_over_q.size() == new_num_tags + 1);
		assert(_sum_n_i_over_j.size() == new_num_tags + 1);
	}
	void InfiniteHMM::_increment_tag_bigram_count(int context_tag, int tag){
		assert(1 <= context_tag && context_tag <= get_num_tags());
		assert(1 <= tag && tag <= get_num_tags());
		_sum_n_i_over_j[context_tag] += 1;
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		bool new_table_generated = false;
		table->add_customer(_beta, new_table_generated);
		if(new_table_generated){
			_increment_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_decrement_tag_bigram_count(int context_tag, int tag){
		assert(1 <= context_tag && context_tag <= get_num_tags());
		assert(1 <= tag && tag <= get_num_tags());
		_sum_n_i_over_j[context_tag] -= 1;
		assert(_sum_n_i_over_j[context_tag] >= 0);
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_increment_tag_word_count(int tag, id word_id){
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_m_i_over_q[tag] += 1;
		Table* table = _m_iq_tables[tag][word_id];
		bool new_table_generated = false;
		table->add_customer(_beta_emission, new_table_generated);
		if(new_table_generated){
			_increment_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::_decrement_tag_word_count(int tag, int word_id){
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_m_i_over_q[tag] += 1;
		Table* table = _m_iq_tables[tag][word_id];
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::_increment_oracle_tag_count(int tag){
		assert(1 <= tag && tag <= get_num_tags());
		_oracle_n_j_counts[tag] += 1;
		_oracle_sum_n_over_j += 1;
	}
	void InfiniteHMM::_decrement_oracle_tag_count(int tag){
		assert(1 <= tag && tag <= get_num_tags());
		_oracle_n_j_counts[tag] -= 1;
		_oracle_sum_n_over_j -= 1;
		assert(_oracle_n_j_counts[tag] >= 0);
		assert(_oracle_sum_n_over_j >= 0);
	}
	void InfiniteHMM::_increment_oracle_word_count(int word_id){
		assert(0 <= word_id < _num_words);
		_oracle_m_q_counts[word_id] += 1;
		_oracle_sum_m_over_q += 1;
	}
	void InfiniteHMM::_decrement_oracle_word_count(int word_id){
		assert(0 <= word_id < _num_words);
		_oracle_m_q_counts[word_id] -= 1;
		_oracle_sum_m_over_q -= 1;
		assert(_oracle_m_q_counts[word_id] >= 0);
		assert(_oracle_sum_m_over_q >= 0);
	}
	int InfiniteHMM::get_sum_n_i_over_j(int tag){
		assert(1 <= tag && tag <= get_num_tags());
		return _sum_n_i_over_j[tag];
	}
	int InfiniteHMM::get_n_ij(int context_tag, int tag){
		if(tag > get_num_tags()){
			return 0;
		}
		assert(1 <= tag && tag <= get_num_tags());
		assert(1 <= context_tag && context_tag <= get_num_tags());
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		return table->get_num_customers();
	}
	int InfiniteHMM::get_oracle_sum_n_over_j(){
		return _oracle_sum_n_over_j;
	}
	int InfiniteHMM::get_oracle_n_j(int tag){
		if(tag > get_num_tags()){
			return 0;
		}
		assert(1 <= tag && tag <= get_num_tags());
		return _oracle_n_j_counts[tag];
	}
	int InfiniteHMM::get_sum_m_i_over_q(int tag){
		assert(1 <= tag && tag <= get_num_tags());
		return _sum_m_i_over_q[tag];
	}
	int InfiniteHMM::get_m_iq(int tag, id word_id){
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id && word_id < get_num_words());
		Table* table = _m_iq_tables[tag][word_id];
		assert(table != NULL);
		return table->get_num_customers();
	}
	int InfiniteHMM::get_oracle_sum_m_over_q(){
		return _oracle_sum_m_over_q;
	}
	int InfiniteHMM::get_oracle_m_q(id word_id){
		assert(0 <= word_id && word_id < get_num_words());
		return _oracle_m_q_counts[word_id];
	}
	// P(s_{t+1}|s_t)
	double InfiniteHMM::compute_p_tag_given_context(int tag, int context_tag){
		double n_i = get_sum_n_i_over_j(context_tag);
		double n_ij = get_n_ij(context_tag, tag);
		assert(n_i >= n_ij);
		double alpha = (tag == context_tag) ? _alpha : 0;
		double empirical_p = (n_ij + alpha) / (n_i + _beta + _alpha);
		double coeff_oracle_p = _beta / (n_i + _beta + _alpha);	// 親の分布を選択する確率. 親からtagが生成される確率とは別物.
		double n_o = get_oracle_sum_n_over_j();
		double n_oj = get_oracle_n_j(tag);
		assert(n_o >= n_oj);
		double T = get_num_tags();
		double g0 = 1.0 / T;
		double oracle_p = (n_oj + _gamma * g0) / (n_o + _gamma);
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	// P(y_t|s_t)
	double InfiniteHMM::compute_p_word_given_tag(id word_id, int tag){
		double m_i = get_sum_m_i_over_q(tag);
		double m_iq = get_m_iq(tag, word_id);
		assert(m_i >= m_iq);
		double empirical_p = m_iq / (m_i + _beta_emission);
		double coeff_oracle_p = _beta_emission / (m_i + _beta_emission);	// 親の分布を選択する確率. 親からword_idが生成される確率とは別物.
		double m_o = get_oracle_sum_m_over_q();
		double m_oq = get_oracle_m_q(word_id);
		assert(m_o >= m_oq);
		double W = get_num_words();
		double g0 = 1.0 / W;
		// cout << "tag = " << tag <<  ", m_i = " << m_i << ", m_iq = " << m_iq << endl;
		// cout << "m_o = " << m_o << ", m_oq = " << m_oq << endl;
		double oracle_p = (m_oq + _gamma_emission * g0) / (m_o + _gamma_emission);
		// cout << "empirical_p = " << empirical_p << ", coeff_oracle_p = " << coeff_oracle_p << ", oracle_p = " << oracle_p << endl;
		return empirical_p + coeff_oracle_p * oracle_p;
	}

}