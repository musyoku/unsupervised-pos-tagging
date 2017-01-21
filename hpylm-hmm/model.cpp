#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <string>
#include <set>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "core/hpylm.h"
#include "core/util.h"
using namespace std;
using namespace boost;

struct value_comparator {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.second > b.second;
	}   
};

typedef struct Word {
	int word_id;
	int tag_id;
} Word;

class PyHpylmHMM{
private:
	HPYLM* _word_hpylm;
	HPYLM** _pos_hpylm_for_tag;
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	vector<vector<Word*>> _dataset;
	vector<int> _rand_indices;
	int _autoincrement;
	int _bos_id;
	int _eos_id;
	int _num_tags;
public:
	PyHpylmHMM(int num_tags){
		// 日本語周り
		// ただのテンプレ
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);

		// HPYLMは全て3-gram
		_word_hpylm = new HPYLM(3);
		_pos_hpylm_for_tag = new HPYLM*[num_tags];
		for(int tag = 0;tag < _num_tags;tag++){
			_pos_hpylm_for_tag[tag] = new HPYLM(3);
		}
		_bos_id = 0;
		_dictionary[_bos_id] = L"<bos>";
		_eos_id = 1;
		_dictionary[_eos_id] = L"<eos>";
		_autoincrement = _eos_id + 1;
	}
	int string_to_word_id(wstring &word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			_dictionary[_autoincrement] = word;
			_dictionary_inv[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	void load_textfile(string filename){
		c_printf("[*]%s\n", (boost::format("%sを読み込んでいます ...") % filename.c_str()).str().c_str());
		wifstream ifs(filename.c_str());
		wstring line_str;
		if (ifs.fail()){
			c_printf("[R]%s [*]%s", "エラー", (boost::format("%sを開けません.") % filename.c_str()).str().c_str());
			exit(1);
		}
		while (getline(ifs, line_str) && !line_str.empty()){
			vector<wstring> word_strs = split_word_by(line_str, ' ');	// スペースで分割
			if(word_strs.size() > 0){
				vector<Word*> words;
				// <bos>
				// 3-gramなので2つ
				for(int n = 0;n < 2;n++){
					Word* bos = new Word();
					bos->word_id = _bos_id;
					bos->tag_id = 0;
					words.push_back(bos);
				}
				for(auto &word_str: word_strs){
					if(word_str.size() == 0){
						continue;
					}
					Word* word = new Word();
					word->word_id = string_to_word_id(word_str);
					word->tag_id = 0;
					words.push_back(word);
				}
				Word* eos = new Word();
				eos->word_id = _eos_id;
				eos->tag_id = 0;
				words.push_back(eos);
				// 訓練データに追加
				_dataset.push_back(words);
			}
		}
		c_printf("[*]%s\n", (boost::format("%sを読み込みました.") % filename.c_str()).str().c_str());
	}
	void load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/hmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _dictionary;
			ifs.close();
		}
		// モデルパラメータの読み込み
		_word_hpylm->load(dirname + "/word.hpylm");
		for(int tag = 0;tag < _num_tags;tag++){
			_pos_hpylm_for_tag[tag]->load((boost::format("%s/pos.%d.hpylm") % dirname.c_str() % tag).str());
		}
	}
	void save(string dirname){
		// 辞書を保存
		string dictionary_filename = dirname + "/hmm.dict";
		std::ofstream ofs(dictionary_filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		ofs.close();
		// モデルパラメータの保存
		_word_hpylm->save(dirname + "/word.hpylm");
		for(int tag = 0;tag < _num_tags;tag++){
			_pos_hpylm_for_tag[tag]->save((boost::format("%s/pos.%d.hpylm") % dirname.c_str() % tag).str());
		}
	}
	void perform_gibbs_sampling(){
		if(_rand_indices.size() != _dataset.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < _dataset.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		shuffle(_rand_indices.begin(), _rand_indices.end(), Sampler::mt);	// データをシャッフル
		for(int n = 0;n < _dataset.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			vector<Word*> &line = _dataset[data_index];
		}
	}
	int get_num_tags(){
		return _num_tags;
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyHpylmHMM>("hpylm_hmm", python::init<int>())
	.def("string_to_word_id", &PyHpylmHMM::string_to_word_id)
	.def("perform_gibbs_sampling", &PyHpylmHMM::perform_gibbs_sampling)
	.def("load", &PyHpylmHMM::load)
	.def("save", &PyHpylmHMM::save)
	.def("get_num_tags", &PyHpylmHMM::get_num_tags)
	.def("load_textfile", &PyHpylmHMM::load_textfile);
}