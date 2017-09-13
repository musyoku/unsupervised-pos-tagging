#include <cassert>
#include <set>
#include <iostream>
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
		_sum_oracle_words_count = 0;
		_sum_oracle_tags_count = 0;
	}
	InfiniteHMM::InfiniteHMM(int initial_num_tags, int num_words){
		InfiniteHMM();
		_initial_num_tags = initial_num_tags;
		_num_words = num_words;
		for(int tag = 0;tag <= _initial_num_tags;tag++){	// <s>と</s>のID=0
			_add_new_tag();
		}
		_prev_num_tags = _initial_num_tags;
		_oracle_word_counts = new int[num_words];
	}
	InfiniteHMM::~InfiniteHMM(){
		for(std::vector<Table*> &tables: _bigram_tag_table){
			for(Table* table: tables){
				delete table;
			}
		}
		for(Table** tables: _tag_word_table){
			for(id word_id = 0;word_id < _num_words;word_id++){
				delete tables[word_id];
			}
			delete[] tables;
		}
		delete[] _oracle_word_counts;
	}
	// <s></s>を除いたタグの数を返す
	int InfiniteHMM::get_num_tags(){
		return _tag_unigram_count.size() - 1;
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
				_increment_tag_unigram_count(ti);
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
		_bigram_tag_table.push_back(table_vec);
		// 0: *, *				0: *, *, *
		// 1: *, *		->   	1: *, *, *
		// 2: *, *, *			2: *, *, *
		for(int context_tag = 1;context_tag < new_num_tags;context_tag++){
			std::vector<Table*> &table_vec = _bigram_tag_table[context_tag];
			assert(table_vec.size() == new_num_tags);
			table_vec.push_back(new Table());
		}
		// tagと単語のペアの出現頻度
		Table** table_array = new Table*[_num_words];
		for(id word_id = 0;word_id < _num_words;word_id++){
			table_array[word_id] = new Table();
		}
		_tag_word_table.push_back(table_array);
		_oracle_tag_counts.push_back(0);
		_tag_unigram_count.push_back(0);
		_sum_word_count_of_tag.push_back(0);
		return new_num_tags;
	}
	// タグが無くなったら詰める
	void InfiniteHMM::_delete_tag(int tag){
		std::vector<Table*> &tables_vec = _bigram_tag_table[tag];
		for(Table* table: tables_vec){
			delete table;
		}
		Table** tables_array = _tag_word_table[tag];
		for(id word_id = 0;word_id < _num_words;word_id++){
			delete tables_array[word_id];
		}
		delete[] tables_array;
		for(int t = tag;t < get_num_tags();t++){
			_bigram_tag_table[t] = _bigram_tag_table[t + 1];
			_tag_word_table[t] = _tag_word_table[t + 1];
			_oracle_tag_counts[t] = _oracle_tag_counts[t + 1];
			_tag_unigram_count[t] = _tag_unigram_count[t + 1];
			_sum_word_count_of_tag[t] = _sum_word_count_of_tag[t + 1];
		}
		_bigram_tag_table.pop_back();
		_tag_word_table.pop_back();
		_oracle_tag_counts.pop_back();
		_tag_unigram_count.pop_back();
		_sum_word_count_of_tag.pop_back();
		assert(_bigram_tag_table.size() == _tag_word_table.size());
		assert(_tag_word_table.size() == _oracle_tag_counts.size());
		assert(_oracle_tag_counts.size() == _tag_unigram_count.size());
		assert(_tag_unigram_count.size() == _sum_word_count_of_tag.size());
	}
	void InfiniteHMM::_increment_tag_unigram_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_tag_unigram_count[tag] += 1;
	}
	void InfiniteHMM::_decrement_tag_unigram_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_tag_unigram_count[tag] -= 1;
		assert(_tag_unigram_count[tag] >= 0);
		if(_tag_unigram_count[tag] == 0){
			_delete_tag(tag);
		}
	}
	void InfiniteHMM::_increment_tag_bigram_count(int context_tag, int tag){
		assert(1 <= context_tag <= get_num_tags());
		assert(1 <= tag <= get_num_tags());
		Table* table = _bigram_tag_table[context_tag][tag];
		assert(table != NULL);
		bool new_table_generated = false;
		table->add_customer(_beta, new_table_generated);
		if(new_table_generated){
			_increment_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_increment_tag_word_count(int tag, id word_id){
		assert(1 <= tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_word_count_of_tag[tag] += 1;
		Table* table = _tag_word_table[tag][word_id];
		bool new_table_generated = false;
		table->add_customer(_beta_emission, new_table_generated);
		if(new_table_generated){
			_increment_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::_increment_oracle_tag_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_oracle_tag_counts[tag] += 1;
		_sum_oracle_tags_count += 1;
	}
	void InfiniteHMM::_increment_oracle_word_count(int word_id){
		assert(0 <= word_id < _num_words);
		_oracle_word_counts[word_id] += 1;
		_sum_oracle_words_count += 1;
	}
	void InfiniteHMM::_decrement_oracle_word_count(int word_id){
		assert(0 <= word_id < _num_words);
		_oracle_word_counts[word_id] -= 1;
		assert(_oracle_word_counts[word_id] >= 0);
		_sum_oracle_words_count -= 1;
		assert(_sum_oracle_words_count >= 0);
	}
	void InfiniteHMM::_decrement_oracle_tag_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_oracle_tag_counts[tag] -= 1;
		assert(_oracle_tag_counts[tag] >= 0);
	}
	void InfiniteHMM::_decrement_tag_bigram_count(int context_tag, int tag){
		assert(1 <= context_tag <= get_num_tags());
		assert(1 <= tag <= get_num_tags());
		Table* table = _bigram_tag_table[context_tag][tag];
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::_decrement_tag_word_count(int tag, int word_id){
		assert(1 <= tag <= get_num_tags());
		assert(0 <= word_id < _num_words);
		_sum_word_count_of_tag[tag] += 1;
		Table* table = _tag_word_table[tag][word_id];
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			_decrement_oracle_word_count(word_id);
		}
	}

}