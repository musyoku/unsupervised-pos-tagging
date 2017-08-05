#pragma once
#include <unordered_map>
#include <vector>

namespace utils{
	void split_word_by(const std::wstring &str, wchar_t delim, std::vector<std::wstring> &elems);
}