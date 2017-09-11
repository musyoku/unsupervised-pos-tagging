#pragma once
#include <boost/python.hpp>
#include <unordered_map>
#include <vector>
#include "../bhmm/common.h"

namespace bhmm {
	class Corpus{
	private:
		void _add_words_to_corpus(std::vector<std::wstring> &word_str_vec);
		void _before_python_add_sentence_str(boost::python::list &py_word_str_list, std::vector<std::wstring> &word_str_vec);
	public:
		std::unordered_map<std::wstring, int> _word_count;		// 事前に単語数がわからないのでmapを使う
		std::vector<std::vector<std::wstring>> _word_sequences;
		int _max_num_words_in_line;
		int _min_num_words_in_line;
		Corpus();
		~Corpus();
		void add_textfile(std::string filename);
		void add_sentence_str(std::wstring sentence_str);
		void python_add_words(boost::python::list py_word_str_list);
		int get_num_words();
		int get_count_of_word(std::wstring word_str);
	};
}