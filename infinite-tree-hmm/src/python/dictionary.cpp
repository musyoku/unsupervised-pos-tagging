#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <fstream>
#include "dictionary.h"

namespace ithmm {
	Dictionary::Dictionary(){
		_id_to_str[ID_UNK] = L"<unk>";
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
	int Dictionary::get_vocab_size(){
		return _id_to_str.size();
	}
	bool Dictionary::load(std::string filename){
		std::ifstream ifs(filename);
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