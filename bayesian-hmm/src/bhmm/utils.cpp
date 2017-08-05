#include <iostream>
#include "utils.h"

namespace bhmm {
	namespace utils{
		template<class T>
		boost::python::list list_from_vector(std::vector<T> &vec){  
			 boost::python::list list;
			 typename std::vector<T>::const_iterator it;
			 for(it = vec.begin(); it != vec.end(); ++it){
				  list.append(*it);
			 }
			 return list;
		}
		template boost::python::list list_from_vector(std::vector<boost::python::tuple> &vec);
		template boost::python::list list_from_vector(std::vector<boost::python::list> &vec);
		void split_word_by(const std::wstring &str, wchar_t delim, std::vector<std::wstring> &word_str_vec){
			word_str_vec.clear();
		    std::wstring word_str;
		    for(wchar_t ch: str){
		        if (ch == delim){
		            if (!word_str.empty()){
		                word_str_vec.push_back(word_str);
		            }
		            word_str.clear();
		        }
		        else{
		            word_str += ch;
		        }
		    }
		    if (!word_str.empty()){
		        word_str_vec.push_back(word_str);
		    }
		}
	}
}