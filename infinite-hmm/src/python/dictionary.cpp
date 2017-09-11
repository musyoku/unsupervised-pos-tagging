#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <fstream>
#include <cassert>
#include "dictionary.h"

namespace ihmm {
	Dictionary::Dictionary(){
		_id_to_str[ID_UNK] = L"<unk>";
		_str_to_id[L"<unk>"] = ID_UNK;
		_autoincrement = ID_UNK + 1;
	}
	id Dictionary::add_word_string(std::wstring word){
		auto itr = _str_to_id.find(word);
		if(itr == _str_to_id.end()){
			_id_to_str[_autoincrement] = word;
			_str_to_id[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	id Dictionary::string_to_word_id(std::wstring word){
		auto itr = _str_to_id.find(word);
		if(itr == _str_to_id.end()){
			return ID_UNK;
		}
		return itr->second;
	}
	std::wstring Dictionary::word_id_to_string(id word_id){
		auto itr = _id_to_str.find(word_id);
		assert(itr != _id_to_str.end());
		return itr->second;
	}
	// <unk>に置き換える
	void Dictionary::remove_ids(std::unordered_set<id> word_ids){
		for(id word_id: word_ids){
			std::wstring word = word_id_to_string(word_id);
			_id_to_str[word_id] = L"<unk>";
			_str_to_id[word] = ID_UNK;
		}
	}
	int Dictionary::get_vocabrary_size(){
		return _str_to_id.size();
	}
	bool Dictionary::is_unk(std::wstring word){
		id word_id = string_to_word_id(word);
		return word_id == ID_UNK;
	}
	bool Dictionary::load(std::string filename){
		std::string dictionary_filename = filename;
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _id_to_str;
			iarchive >> _str_to_id;
			iarchive >> _autoincrement;
			ifs.close();
			return true;
		}
		return false;
	}
	bool Dictionary::save(std::string filename){
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _id_to_str;
		oarchive << _str_to_id;
		oarchive << _autoincrement;
		ofs.close();
		return true;
	}
}