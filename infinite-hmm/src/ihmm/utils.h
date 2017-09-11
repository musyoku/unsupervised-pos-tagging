#pragma once
#include <boost/python.hpp>
#include <unordered_map>
#include <vector>

namespace bhmm {
	namespace utils{
		template<class T>
		boost::python::list list_from_vector(std::vector<T> &vec);
		template<class T>
		std::vector<T> vector_from_list(boost::python::list &list);
		void split_word_by(const std::wstring &str, wchar_t delim, std::vector<std::wstring> &word_str_vec);
	}
}