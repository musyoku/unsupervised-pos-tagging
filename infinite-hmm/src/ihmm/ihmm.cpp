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
		_alpha = 0;
		_beta = 1;
		_gamma = 1;
		_beta_emission = 1;
		_gamma_emission = 1;
		_oracle_sum_n_over_j = 0;
		_oracle_sum_m_over_q = 0;
		_gibbs_sampling_table = NULL;
	}
	InfiniteHMM::InfiniteHMM(int initial_num_tags, int num_words): InfiniteHMM(){
		_initial_num_tags = initial_num_tags;
		_num_words = num_words;
		for(int tag = 0;tag <= _initial_num_tags;tag++){
			_add_new_tag();
		}
		_prev_num_tags = _initial_num_tags;
		_oracle_m_q_counts = new int[num_words];
		for(int word_id = 0;word_id < num_words;word_id++){
			_oracle_m_q_counts[word_id] = 0;
		}
	}
	InfiniteHMM::~InfiniteHMM(){
		for(int tag = get_num_tags();tag > 0;tag--){
			_delete_tag(tag);
		}
		delete[] _oracle_m_q_counts;
		if(_gibbs_sampling_table != NULL){
			delete[] _gibbs_sampling_table;
		}
	}
	// <s></s>を除いたタグの数を返す
	int InfiniteHMM::get_num_tags(){
		return _oracle_n_j_counts.size() - 1;
	}
	int InfiniteHMM::get_num_words(){
		return _num_words;
	}
	int InfiniteHMM::get_num_valid_tags(){
		int num_tags = 0;
		for(int tag = 1;tag <= get_num_tags();tag++){
			int count = _oracle_n_j_counts[tag];
			if(count > 0){
				num_tags++;
			}
		}
		return num_tags;
	}
	void InfiniteHMM::initialize_with_training_dataset(std::vector<std::vector<Word*>> &dataset){
		// 最初は品詞をランダムに割り当てる
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &word_ids = dataset[data_index];
			int ti_1 = 0;
			for(int i = 1;i < word_ids.size() - 1;i++){
				Word* word = word_ids[i];












































				int ti = sampler::uniform_int(1, _initial_num_tags);


				ti = word->_tag;




				int wi = word->_id;
				_increment_tag_bigram_count(ti_1, ti);
				_increment_tag_word_count(ti, wi);
				// word->_tag = ti;
				ti_1 = ti;
			}
			_increment_tag_bigram_count(ti_1, 0);	// </s>への遷移
		}
	}
	void InfiniteHMM::_remove_all_training_dataset(std::vector<std::vector<Word*>> &dataset){
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &word_ids = dataset[data_index];
			int ti_1 = 0;
			for(int i = 1;i < word_ids.size() - 1;i++){
				Word* word = word_ids[i];
				int ti = word->_tag;
				int wi = word->_id;
				_decrement_tag_bigram_count(ti_1, ti);
				_decrement_tag_word_count(ti, wi);
				ti_1 = ti;
			}
			_decrement_tag_bigram_count(ti_1, 0);	// </s>への遷移
		}
	}
	// void InfiniteHMM::_add_special_tag(){
	// 	assert(get_num_tags() == -1);
	// 	std::vector<Table*> table_vec;
	// 	table_vec.push_back(NULL);
	// 	_n_ij_tables.push_back(table_vec);
	// 	_m_iq_tables.push_back(NULL);
	// 	_oracle_n_j_counts.push_back(0);
	// 	_sum_m_i_over_q.push_back(0);
	// 	_sum_n_i_over_j.push_back(0);
	// 	assert(_n_ij_tables.size() == 1);
	// 	assert(_m_iq_tables.size() == 1);
	// 	assert(_oracle_n_j_counts.size() == 1);
	// 	assert(_sum_m_i_over_q.size() == 1);
	// 	assert(_sum_n_i_over_j.size() == 1);
	// }
	int InfiniteHMM::_add_new_tag(){
		int new_num_tags = get_num_tags() + 1;	// <s></s>は除く
		assert(new_num_tags >= 0);
		// bigramのカウント
		// 0: *, *				0: *, *
		// 1: *, *		->   	1: *, *
		// 						2: *, *, *
		std::vector<Table*> table_vec;
		for(int tag = 0;tag <= new_num_tags;tag++){
			table_vec.push_back(new Table(tag));
		}
		_n_ij_tables.push_back(table_vec);
		// 0: *, *				0: *, *, *
		// 1: *, *		->   	1: *, *, *
		// 2: *, *, *			2: *, *, *
		for(int context_tag = 0;context_tag < new_num_tags;context_tag++){
			std::vector<Table*> &table_vec = _n_ij_tables[context_tag];
			assert(table_vec.size() == new_num_tags);
			table_vec.push_back(new Table(new_num_tags));
		}
		// tagと単語のペアの出現頻度
		Table** table_array = new Table*[_num_words];
		for(int word_id = 0;word_id < _num_words;word_id++){
			table_array[word_id] = new Table(word_id);
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

		if(_gibbs_sampling_table != NULL){
			delete[] _gibbs_sampling_table;
		}
		_gibbs_sampling_table = new double[new_num_tags + 2];	// <s>と新しいタグの分も確保
		return new_num_tags;
	}
	// タグが無くなったら詰める
	void InfiniteHMM::_delete_tag(int delete_tag){
		assert(1 <= delete_tag && delete_tag <= get_num_tags());
		int current_num_tags = get_num_tags();
		int new_num_tags = get_num_tags() - 1;	// <s></s>は除く
		// 対応する行を削除
		std::vector<Table*> &tables_vec = _n_ij_tables[delete_tag];
		for(Table* &table: tables_vec){
			assert(table != NULL);
			delete table;
			table = NULL;
		}
		// 対応する列を削除
		for(int context_tag = 0;context_tag <= current_num_tags;context_tag++){
			if(context_tag == delete_tag){
				continue;	// 交差する点は二重に解放されるのでスキップ
			}
			Table* &table = _n_ij_tables[context_tag][delete_tag];
			assert(table != NULL);
			delete table;
			table = NULL;
			for(int t = delete_tag;t <= new_num_tags;t++){
				Table* table = _n_ij_tables[context_tag][t + 1];
				table->_identifier -= 1;
				_n_ij_tables[context_tag][t] = table;
			}
			_n_ij_tables[context_tag].pop_back();
		}
		Table** tables_array = _m_iq_tables[delete_tag];
		for(int word_id = 0;word_id < _num_words;word_id++){
			delete tables_array[word_id];
		}
		delete[] tables_array;
		for(int t = delete_tag;t <= new_num_tags;t++){
			_n_ij_tables[t] = _n_ij_tables[t + 1];
			_m_iq_tables[t] = _m_iq_tables[t + 1];
			_oracle_n_j_counts[t] = _oracle_n_j_counts[t + 1];
			_sum_m_i_over_q[t] = _sum_m_i_over_q[t + 1];
			_sum_n_i_over_j[t] = _sum_n_i_over_j[t + 1];
			assert(_n_ij_tables[t].size() == new_num_tags + 1);
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
		assert(0 <= context_tag && context_tag <= get_num_tags());
		assert(0 <= tag && tag <= get_num_tags());
		_sum_n_i_over_j[context_tag] += 1;
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		assert(table->_identifier == tag);
		bool new_table_generated = false;
		table->add_customer(_beta, new_table_generated);
		if(new_table_generated){
			_increment_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_decrement_tag_bigram_count(int context_tag, int tag){
		assert(0 <= context_tag && context_tag <= get_num_tags());
		assert(0 <= tag && tag <= get_num_tags());
		_sum_n_i_over_j[context_tag] -= 1;
		assert(_sum_n_i_over_j[context_tag] >= 0);
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		assert(table->_identifier == tag);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_increment_tag_word_count(int tag, int word_id){
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_m_i_over_q[tag] += 1;
		Table* table = _m_iq_tables[tag][word_id];
		assert(table != NULL);
		assert(table->_identifier == word_id);
		bool new_table_generated = false;
		table->add_customer(_beta_emission, new_table_generated);
		if(new_table_generated){
			_increment_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::_decrement_tag_word_count(int tag, int word_id){
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_m_i_over_q[tag] -= 1;
		Table* table = _m_iq_tables[tag][word_id];
		assert(table != NULL);
		assert(table->_identifier == word_id);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::_increment_oracle_tag_count(int tag){
		assert(0 <= tag && tag <= get_num_tags());
		_oracle_n_j_counts[tag] += 1;
		_oracle_sum_n_over_j += 1;
	}
	void InfiniteHMM::_decrement_oracle_tag_count(int tag){
		assert(0 <= tag && tag <= get_num_tags());
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
		if(tag > get_num_tags()){
			return 0;
		}
		assert(0 <= tag && tag <= get_num_tags());
		return _sum_n_i_over_j[tag];
	}
	int InfiniteHMM::get_n_ij(int context_tag, int tag){
		if(tag > get_num_tags()){
			return 0;
		}
		if(context_tag > get_num_tags()){
			return 0;
		}
		assert(0 <= tag && tag <= get_num_tags());
		assert(0 <= context_tag && context_tag <= get_num_tags());
		Table* table = _n_ij_tables[context_tag][tag];
		assert(table != NULL);
		assert(table->_identifier == tag);
		return table->get_num_customers();
	}
	int InfiniteHMM::get_oracle_sum_n_over_j(){
		return _oracle_sum_n_over_j;
	}
	int InfiniteHMM::get_oracle_n_j(int tag){
		if(tag > get_num_tags()){
			return 0;
		}
		assert(0 <= tag && tag <= get_num_tags());
		return _oracle_n_j_counts[tag];
	}
	int InfiniteHMM::get_sum_m_i_over_q(int tag){
		if(tag > get_num_tags()){
			return 0;
		}
		assert(1 <= tag && tag <= get_num_tags());
		return _sum_m_i_over_q[tag];
	}
	int InfiniteHMM::get_m_iq(int tag, int word_id){
		if(tag > get_num_tags()){
			return 0;
		}
		assert(1 <= tag && tag <= get_num_tags());
		assert(0 <= word_id && word_id < get_num_words());
		Table* table = _m_iq_tables[tag][word_id];
		assert(table != NULL);
		assert(table->_identifier == word_id);
		return table->get_num_customers();
	}
	int InfiniteHMM::get_oracle_sum_m_over_q(){
		return _oracle_sum_m_over_q;
	}
	int InfiniteHMM::get_oracle_m_q(int word_id){
		assert(0 <= word_id && word_id < get_num_words());
		return _oracle_m_q_counts[word_id];
	}
	int InfiniteHMM::_get_new_tag(){
		for(int tag = 1;tag < get_num_tags();tag++){
			if(_oracle_n_j_counts[tag] == 0){
				return tag;
			}
		}
		return get_num_tags() + 1;
	}
	bool InfiniteHMM::is_tag_new(int tag){
		assert(tag != 0);
		if(tag > get_num_tags()){
			return true;
		}
		if(_oracle_n_j_counts[tag] == 0){
			return true;
		}
		return false;
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
		double g0 = 1.0 / (T + 1);
		double oracle_p = (n_oj + _gamma * g0) / (n_o + _gamma);
		std::cout << "		compute_p_tag_given_context(" << tag << ", " << context_tag << std::endl;
		std::cout << "		n_i=" << n_i << ", n_ij=" << n_ij << ", empirical_p=" << empirical_p << ", coeff_oracle_p=" << coeff_oracle_p << std::endl;
		std::cout << "		n_o=" << n_o << ", n_oj=" << n_oj << ", T=" << T << ", oracle_p=" << oracle_p << std::endl;
		return empirical_p + coeff_oracle_p * oracle_p;
	}
	// P(y_t|s_t)
	double InfiniteHMM::compute_p_word_given_tag(int word_id, int tag){
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

		std::cout << "		compute_p_word_given_tag(" << word_id << ", " << tag << std::endl;
		std::cout << "		m_i=" << m_i << ", m_iq=" << m_iq << ", empirical_p=" << empirical_p << ", coeff_oracle_p=" << coeff_oracle_p << std::endl;
		std::cout << "		m_o=" << m_o << ", m_oq=" << m_oq << ", W=" << W << ", oracle_p=" << oracle_p << std::endl;

		return empirical_p + coeff_oracle_p * oracle_p;
	}
	// ギブスサンプリング
	int InfiniteHMM::_perform_gibbs_sampling_on_markov_blanket(int ti_1, int ti1, int wi){
		std::cout << "_perform_gibbs_sampling_on_markov_blanket(" << ti_1 << ", " << ti1 << ", " << wi << std::endl;
		double sum = 0;
		for(int tag = 1;tag <= get_num_tags();tag++){
			double p_emission = compute_p_word_given_tag(wi, tag);
			double p_generation = compute_p_tag_given_context(tag, ti_1);
			// int correcting_count_for_bigram = (ti_1 == tag == ti1) ? 1 : 0;
			// int correcting_count_for_destination = (ti_1 == tag) ? 1 : 0;
			double p_likelihood = compute_p_tag_given_context(ti1, tag);
			assert(p_emission > 0);
			assert(p_generation > 0);
			assert(p_likelihood > 0);
			double conditional_p = p_emission * p_generation * p_likelihood;
			_gibbs_sampling_table[tag] = conditional_p;
			std::cout << "	_gibbs_sampling_table[" << tag << "] = " << conditional_p << std::endl;
			sum += conditional_p;
		}
		int stack_size = get_num_tags();
		int new_tag = _get_new_tag();
		assert(1 <= new_tag && new_tag <= get_num_tags() + 1);
		if(new_tag == get_num_tags() + 1){
			double p_emission = compute_p_word_given_tag(wi, new_tag);
			double p_generation = compute_p_tag_given_context(new_tag, ti_1);
			double p_likelihood = compute_p_tag_given_context(ti1, new_tag);
			assert(p_emission > 0);
			assert(p_generation > 0);
			assert(p_likelihood > 0);
			double conditional_p = p_emission * p_generation * p_likelihood;
			_gibbs_sampling_table[new_tag] = conditional_p;
			std::cout << "	_gibbs_sampling_table[" << new_tag << "] = " << conditional_p << std::endl;
			sum += conditional_p;
			stack_size += 1;
		}
		assert(sum > 0);
		double normalizer = 1.0 / sum;
		double bernoulli = sampler::uniform(0, 1);
		double stack = 0;
		std::cout << "bernoulli = " << bernoulli << std::endl;
		for(int tag = 1;tag <= stack_size;tag++){
			stack += _gibbs_sampling_table[tag] * normalizer;
			if(bernoulli <= stack){
				return tag;
			}
		}
		assert(false);
		return stack_size;
	}
	void InfiniteHMM::perform_gibbs_sampling_with_sequence(std::vector<Word*> &word_vec){
		for(int i = 1;i < word_vec.size() - 1;i++){	// <s>と</s>の内側だけ考える
			int ti_1 = word_vec[i - 1]->_tag;
			int ti = word_vec[i]->_tag;
			int wi = word_vec[i]->_id;
			int ti1 = word_vec[i + 1]->_tag;
			// 現在のtiをモデルから除去
			_decrement_tag_bigram_count(ti_1, ti);
			_decrement_tag_bigram_count(ti, ti1);
			_decrement_tag_word_count(ti, wi);
			// t_iを再サンプリング
			int new_ti = _perform_gibbs_sampling_on_markov_blanket(ti_1, ti1, wi);
			std::cout << "new_ti = " << new_ti << std::endl;
			assert(1 <= new_ti && new_ti <= get_num_tags() + 1);	// 新しいタグも許可
			if(new_ti == get_num_tags() + 1){
				_add_new_tag();
			}
			// 新しいt_iをモデルパラメータに追加
			_increment_tag_bigram_count(ti_1, new_ti);
			_increment_tag_bigram_count(new_ti, ti1);
			_increment_tag_word_count(new_ti, wi);
			word_vec[i]->_tag = new_ti;
		}
	}
	bool InfiniteHMM::save(std::string filename){
		return true;
	}
	bool InfiniteHMM::load(std::string filename){
		return true;
	}
}