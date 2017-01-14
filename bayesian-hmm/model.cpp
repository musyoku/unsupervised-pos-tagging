#include <boost/python.hpp>
#include <boost/format.hpp>
#include <string>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "core/bhmm.h"
#include "util.h"
using namespace boost;

class PyBayesianHMM{
private:
	BayesianHMM* _hmm;
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	vector<vector<int>> _dataset;
	int _max_epoch;
	int _autoincrement;
public:
	PyBayesianHMM(int vector_length){
		// 日本語周り
		// ただのテンプレ
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);
		
		_hmm = BayesianHMM();
		_autoincrement = 0;
	}
	int string_to_token_id(wstring &word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			_dictionary[_autoincrement] = word;
			_dictionary_inv[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr.second;
	}
	void load_textfile(string filename){
		wifstream ifs(filename.c_str());
		wstring str;
		if (ifs.fail()){
			c_printf("[R]%s [*]%s", "エラー", (boost::format("%sを開けません.") % filename.c_str()).str().c_str());
			exit(1);
		}
		while (getline(ifs, str) && !str.empty()){
			vector<wstring> words = split_word_by(str, ' ');	// スペースで分割
			if(words.size() > 0){
				vector<int> token_ids;
				for(auto &word: words){
					if(word.size() == 0){
						continue;
					}
					int token_id = string_to_token_id(word);
					token_ids.push_back(token_id);
				}
				_dataset.push_back(token_ids);
			}
		}
		c_printf("[*]%s", (boost::format("%sを読み込みました. (%d行)") % filename.c_str() % _word2vec->_dataset.size()).str().c_str());
	}
	void train(){
		c_printf("[*]%s\n", (boost::format("学習を開始します. (スレッド数: %d)") % num_threads).str().c_str());
		for(int epoch = 0;epoch < _max_epoch;epoch++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
		}
	}
	void set_num_phrases(int number){
		_hmm->set_num_phrases(lr);
	}
	void set_num_words(int number){
		_hmm->set_num_words(number);
	}
	void set_max_epoch(int epoch){
		_max_epoch = epoch;
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyBayesianHMM>("word2vec", python::init<int>())
	.def("string_to_token_id", &PyBayesianHMM::string_to_token_id)
	.def("train", &PyBayesianHMM::train)
	.def("set_num_phrases", &PyBayesianHMM::set_num_phrases)
	.def("set_num_words", &PyBayesianHMM::set_num_words)
	.def("set_max_epoch", &PyBayesianHMM::set_max_epoch)
	.def("load_textfile", &PyBayesianHMM::load_textfile)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}