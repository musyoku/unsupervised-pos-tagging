#include <iostream>
#include "utils.h"

namespace ithmm {
	namespace utils{
		void split_word_by(const std::wstring &str, wchar_t delim, std::vector<std::wstring> &elems){
			elems.clear();
			std::wstring item;
			for(wchar_t ch: str){
				if (ch == delim){
					if (!item.empty()){
						elems.push_back(item);
					}
					item.clear();
				}
				else{
					item += ch;
				}
			}
			if (!item.empty()){
				elems.push_back(item);
			}
		}
	}
}