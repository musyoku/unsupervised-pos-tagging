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
#include "core/ihmm.h"
#include "core/util.h"
using namespace std;
using namespace boost;

struct value_comparator {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.second > b.second;
	}   
};

class PyBayesianHMM{
private:
	InfiniteHMM* _hmm;
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	vector<vector<Word*>> _dataset;
	vector<int> _rand_indices;
	int _autoincrement;
	int _bos_id;
	int _eos_id;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
public:
	PyBayesianHMM(int initial_num_tags){
		// 日本語周り
		// ただのテンプレ
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);

		_hmm = new InfiniteHMM(initial_num_tags);
		_bos_id = 0;
		_dictionary[_bos_id] = L"<bos>";
		_eos_id = 1;
		_dictionary[_eos_id] = L"<eos>";
		_autoincrement = _eos_id + 1;

		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
	}
	int string_to_word_id(wstring word){
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
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			add_line(line_str);
		}
		c_printf("[*]%s\n", (boost::format("%sを読み込みました.") % filename.c_str()).str().c_str());
	}
	void add_line(wstring line_str){
		vector<wstring> word_strs = split_word_by(line_str, L' ');	// スペースで分割
		int num_words = word_strs.size();
		if(num_words > _max_num_words_in_line){
			_max_num_words_in_line = num_words;
		}
		if(num_words < _max_num_words_in_line || _min_num_words_in_line == -1){
			_min_num_words_in_line = num_words;
		}
		if(word_strs.size() > 0){
			vector<Word*> words;
			// <bos>
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
			// <eos>も2つ追加しておくとt_{i+1}, t_{i+2}が常に存在するのでギブスサンプリング時に場合分けしなくてもいいかもしれない
			for(int n = 0;n < 2;n++){
				Word* eos = new Word();
				eos->word_id = _eos_id;
				eos->tag_id = 0;
				words.push_back(eos);
			}
			// 訓練データに追加
			_dataset.push_back(words);
		}
	}
	void initialize(){
		_hmm->initialize(_dataset);
	}
	bool load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/hmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _dictionary;
			iarchive >> _dictionary_inv;
			iarchive >> _autoincrement;
			ifs.close();
		}
		return _hmm->load(dirname);
	}
	bool save(string dirname){
		// 辞書を保存
		std::ofstream ofs(dirname + "/hmm.dict");
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		oarchive << _dictionary_inv;
		oarchive << _autoincrement;
		ofs.close();
		return _hmm->save(dirname);
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
			_hmm->perform_gibbs_sampling_with_line(line);
		}
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyBayesianHMM>("bayesian_hmm", python::init<int>())
	.def("string_to_word_id", &PyBayesianHMM::string_to_word_id)
	.def("perform_gibbs_sampling", &PyBayesianHMM::perform_gibbs_sampling)
	.def("initialize", &PyBayesianHMM::initialize)
	.def("load", &PyBayesianHMM::load)
	.def("save", &PyBayesianHMM::save)
	.def("add_line", &PyBayesianHMM::add_line)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}