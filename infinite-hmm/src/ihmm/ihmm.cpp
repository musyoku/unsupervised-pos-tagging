#include <cassert>
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
		_sum_oracle_words_count = 0;
		_sum_oracle_tags_count = 0;
	}
	InfiniteHMM::InfiniteHMM(int initial_num_tags, int num_words){
		InfiniteHMM();
		_initial_num_tags = initial_num_tags;
		_num_words = num_words;
		for(int tag = 0;tag <= _initial_num_tags;tag++){	// <s>と</s>のID=0
			_tag_unigram_count.push_back(0);
		}
		_prev_num_tags = _initial_num_tags;
	}
	InfiniteHMM::~InfiniteHMM(){
		for(auto &tables: _bigram_tag_table){
			for(auto &table: tables){
				delete table.second;
			}
		}
		for(auto &tables: _tag_word_table){
			for(auto &table: tables){
				delete table.second;
			}
		}
	}
	// <s></s>を除いたタグの数を返す
	int InfiniteHMM::get_num_tags(){
		return _tag_unigram_count.size() - 1;
	}
	void InfiniteHMM::initialize_with_training_corpus(std::vector<std::vector<Word*>> &dataset){
		// 最初は品詞をランダムに割り当てる
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &word_ids = dataset[data_index];
			int ti_1 = 0;
			for(int i = 0;i < word_ids.size();i++){
				Word* word = word_ids[i];
				int ti = sampler::uniform_int(1, _initial_num_tags);
				int wi = word->_id;
				increment_tag_bigram_count(ti_1, ti);
				increment_tag_unigram_count(ti);
				increment_tag_word_count(ti, wi);
				word->_tag = ti;
				ti_1 = ti;
			}
			increment_tag_bigram_count(ti_1, 0);	// </s>への遷移
		}
	}
	// タグが無くなったら詰める
	void InfiniteHMM::_delete_tag(int tag){
		auto &tables = _bigram_tag_table[t];
		for(auto &table: tables){
			delete table.second;
		}
		auto &tables = _tag_word_table[t];
		for(auto &table: tables){
			delete table.second;
		}
		for(int t = tag;t < get_num_tags();t++){
			_bigram_tag_table[t] = _bigram_tag_table[t + 1];
			_tag_word_table[t] = _tag_word_table[t + 1];
			_oracle_word_counts[t] = _oracle_word_counts[t + 1];
			_oracle_tag_counts[t] = _oracle_tag_counts[t + 1];
			_tag_unigram_count[t] = _tag_unigram_count[t + 1];
			_sum_word_count_of_tag[t] = _sum_word_count_of_tag[t + 1];
		}
		_bigram_tag_table.pop_back();
		_tag_word_table.pop_back();
		_oracle_word_counts.pop_back();
		_oracle_tag_counts.pop_back();
		_tag_unigram_count.pop_back();
		_sum_word_count_of_tag.pop_back();
		assert(_bigram_tag_table.size() == _tag_word_table);
		assert(_tag_word_table.size() == _oracle_word_counts);
		assert(_oracle_word_counts.size() == _oracle_tag_counts);
		assert(_oracle_tag_counts.size() == _tag_unigram_count);
		assert(_tag_unigram_count.size() == _sum_word_count_of_tag);
	}
	void InfiniteHMM::increment_tag_unigram_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_tag_unigram_count[tag] += 1;
	}
	void InfiniteHMM::decrement_tag_unigram_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_tag_unigram_count[tag] -= 1;
		assert(_tag_unigram_count[tag] >= 0);
		if(_tag_unigram_count[tag] == 0){
			_delete_tag(tag);
		}
	}
	void InfiniteHMM::increment_tag_bigram_count(int context_tag, int tag){
		assert(1 <= context_tag <= get_num_tags());
		assert(1 <= tag <= get_num_tags());
		Table* table = _bigram_tag_table[context_tag][tag];
		assert(table != NULL);
		bool new_table_generated = false;
		table->add_customer(_beta, new_table_generated);
		if(new_table_generated){
			increment_oracle_tag_count(tag);
		}
	}
	void InfiniteHMM::increment_tag_word_count(Word* word){
		increment_tag_word_count(word->_tag, word->_id);
	}
	void InfiniteHMM::increment_tag_word_count(int tag, id word_id){
		assert(1 <= tag <= get_num_tags());
		_sum_word_count_of_tag[tag] += 1;
		unordered_map<int, Table*> &tables = _tag_word_table[tag];
		auto itr_table = tables.find(word_id);
		if(itr_table == tables.end()){
			table = new Table(word_id);
			tables[word_id] = table;
		}else{
			table = itr_table->second;
		}
		bool new_table_generated = false;
		table->add_customer(_beta_emission, new_table_generated);
		if(new_table_generated){
			increment_oracle_word_count(word_id);
		}
	}
	void InfiniteHMM::increment_oracle_tag_count(int tag){
		assert(1 <= tag <= get_num_tags());
		_oracle_tag_counts[tag] += 1;
		_sum_oracle_tags_count += 1;
	}
	void InfiniteHMM::increment_oracle_word_count(int word_id){
		_oracle_word_counts[word_id] += 1;
		_sum_oracle_words_count += 1;
	}
	void InfiniteHMM::decrement_oracle_word_count(int word_id){
		_oracle_word_counts[word_id] -= 1;
		assert(_oracle_word_counts[word_id] >= 0);
		_sum_oracle_words_count -= 1;
		assert(_sum_oracle_words_count >= 0);
	}
	void InfiniteHMM::decrement_oracle_tag_count(int tag){
		assert(1 <= tag <= get_num_tags());
		auto itr = _oracle_tag_counts.find(tag);
		assert(itr != _oracle_tag_counts.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_oracle_tag_counts.erase(itr);
		}
		_sum_oracle_tags_count -= 1;
		assert(_sum_oracle_tags_count >= 0);
	}
	void InfiniteHMM::decrement_tag_bigram_count(int context_tag, int tag){
		auto itr = _sum_bigram_destination.find(context_tag);
		assert(itr != _sum_bigram_destination.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_sum_bigram_destination.erase(itr);
		}

		auto itr_context = _bigram_tag_table.find(context_tag);
		assert(itr_context != _bigram_tag_table.end());
		unordered_map<int, Table*> &tables = itr_context->second;
		auto itr_table = tables.find(tag);
		assert(itr_table != tables.end());
		Table* table = itr_table->second;
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(empty_table_deleted){
			decrement_oracle_tag_count(tag);
		}
		if(table->is_empty()){
			delete table;
			tables.erase(itr_table);
		}
		if(tables.size() == 0){
			_bigram_tag_table.erase(itr_context);
		}
	}
	void InfiniteHMM::decrement_tag_word_count(int tag, int word_id){
		auto itr_sum = _sum_word_count_of_tag.find(tag);
		assert(itr_sum != _sum_word_count_of_tag.end());
		itr_sum->second -= 1;
		assert(itr_sum->second >= 0);
		if(itr_sum->second == 0){
			_sum_word_count_of_tag.erase(itr_sum);
		}

		auto itr_tag = _tag_word_table.find(tag);
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
			delete table;
			tables.erase(itr_table);
		}
		if(tables.size() == 0){
			_tag_word_table.erase(itr_tag);
		}
	}

}