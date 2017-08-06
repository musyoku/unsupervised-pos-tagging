#pragma once
#include <boost/python.hpp>
#include <unordered_map>
#include <vector>
#include "../bhmm/common.h"
#include "dictionary.h"

namespace bhmm {
	class Dataset{
	public:
		Dictionary* _dict;
		std::unordered_map<id, int> _word_count;		// 事前に単語数がわからないのでmapを使う
		std::vector<std::vector<Word*>> _word_sequences_train;
		std::vector<std::vector<Word*>> _word_sequences_dev;
		int _max_num_words_in_line;
		int _min_num_words_in_line;
		Dataset(Dictionary* dict);
		~Dataset();
		void add_textfile(std::string filename, double train_split_ratio);
		void add_sentence_str_train(std::wstring sentence_str);
		void add_sentence_str_dev(std::wstring sentence_str);
		void _before_python_add_sentence_str(boost::python::list &py_word_str_list, std::vector<std::wstring> &word_str_vec);
		void python_add_words_train(boost::python::list py_word_str_list);
		void python_add_words_dev(boost::python::list py_word_str_list);
		void _add_words_to_dataset(std::vector<std::wstring> &word_str_vec, std::vector<std::vector<Word*>> &dataset);
		void mark_low_frequency_words_as_unknown(int threshold = 1);
		void _mark_low_frequency_words_as_unknown(int threshold, std::vector<std::vector<Word*>> &word_sequence_vec);
		int get_num_words();
		int get_count_of_word(id word_id);
	};
}