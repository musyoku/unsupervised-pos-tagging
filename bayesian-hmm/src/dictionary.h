#pragma once

#define ID_BOS 0
#define ID_EOS 1
#define ID_UNK 2

class Dictionary{
public:
	std::unordered_map<id, std::wstring> _id_to_str;
	std::unordered_map<std::wstring, id> _str_to_id;
	id _autoincrement;
	Dictionary(){
		_id_to_str[ID_BOS] = L"<bos>";
		_id_to_str[ID_EOS] = L"<eos>";
		_id_to_str[ID_UNK] = L"<unk>";
		_autoincrement = ID_UNK + 1;
	}
	id get_eos_id(){
		return ID_EOS;
	}
	id add_word_string(std::wstring word){
		auto itr = _str_to_id.find(word);
		if(itr == _str_to_id.end()){
			_id_to_str[_autoincrement] = word;
			_str_to_id[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	id string_to_word_id(std::wstring word){
		auto itr = _str_to_id.find(word);
		if(itr == _str_to_id.end()){
			return ID_UNK;
		}
		return itr->second;
	}
	std::wstring word_id_to_string(id word_id){
		auto itr = _id_to_str.find(word_id);
		assert(itr != _id_to_str.end());
		return itr->second;
	}
	int get_vocabrary_size(){
		return _str_to_id.size();
	}
	bool load(std::string filename){
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
	bool save(std::string filename){
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _id_to_str;
		oarchive << _str_to_id;
		oarchive << _autoincrement;
		ofs.close();
		return true;
	}
};