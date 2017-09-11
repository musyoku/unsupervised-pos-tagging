#pragma once
#include <boost/python.hpp>
#include <unordered_map>
#include <vector>
#include "../bhmm/common.h"
#include "corpus.h"
#include "dictionary.h"

namespace bhmm {
	class Dataset{
	private:
		void _before_python_add_sentence_str(boost::python::list &py_word_str_list, std::vector<std::wstring> &word_str_vec);
		void _add_words_to_dataset(std::vector<std::wstring> &word_str_vec, std::vector<std::vector<Word*>> &dataset, Corpus* corpus, int unknown_count);
		void _mark_low_frequency_words_as_unknown(int threshold, std::vector<std::vector<Word*>> &word_sequence_vec);
	public:
		Dictionary* _dict;
		std::vector<std::vector<Word*>> _word_sequences_train;
		std::vector<std::vector<Word*>> _word_sequences_dev;
		int _max_num_words_in_line;
		int _min_num_words_in_line;
		Dataset(Corpus* corpus, double dev_split, int unknown_count);
		~Dataset();
		int get_num_words();
		int get_count_of_word(id word_id);
		Dictionary &get_dict_obj();
	};
}