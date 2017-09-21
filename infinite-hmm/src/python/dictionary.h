#pragma once
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "../ihmm/common.h"

#define ID_UNK 0

namespace ihmm {
	class Dictionary{
	public:
		std::unordered_map<int, std::wstring> _id_to_str;
		std::unordered_map<std::wstring, int> _str_to_id;
		int _autoincrement;
		Dictionary();
		int add_word_string(std::wstring word);
		int string_to_word_id(std::wstring word);
		std::wstring word_id_to_string(int word_id);
		void remove_ids(std::unordered_set<int> word_ids);
		int get_vocabrary_size();
		bool is_string_unk(std::wstring word);
		bool is_id_unk(int word_id);
		bool load(std::string filename);
		bool save(std::string filename);
	};
}