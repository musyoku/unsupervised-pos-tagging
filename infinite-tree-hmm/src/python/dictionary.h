#pragma once
#include <boost/python.hpp>
#include <unordered_map>
#include <vector>
#include "../ithmm/common.h"

#define ID_BOS 0
#define ID_EOS 1
#define ID_UNK 2

class Dictionary{
public:
	std::unordered_map<id, std::wstring> _id_to_str;
	std::unordered_map<std::wstring, id> _str_to_id;
	id _autoincrement;
	Dictionary();
	id get_eos_id();
	id add_word_string(std::wstring word);
	id string_to_word_id(std::wstring word);
	bool load(std::string filename);
	bool save(std::string filename);
};