#pragma once
#include <unordered_map>
#include <string>
#include "../bhmm/common.h"

#define ID_UNK 0

namespace bhmm {
	class Dictionary{
	public:
		std::unordered_map<id, std::wstring> _id_to_str;
		std::unordered_map<std::wstring, id> _str_to_id;
		id _autoincrement;
		Dictionary();
		id add_word_string(std::wstring word);
		id string_to_word_id(std::wstring word);
		std::wstring word_id_to_string(id word_id);
		int get_vocabrary_size();
		bool is_unk(std::wstring word);
		bool load(std::string filename);
		bool save(std::string filename);
	};
}