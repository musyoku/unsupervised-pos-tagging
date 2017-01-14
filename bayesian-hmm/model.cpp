#include <boost/python.hpp>
#include <boost/format.hpp>
#include <string>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "core/bhmm.h"
#include "core/util.h"
using namespace std;
using namespace boost;

class PyBayesianHMM{
private:
	BayesianHMM* _hmm;
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	vector<vector<Word*>> _dataset;
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

		_hmm = new BayesianHMM();
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
		return itr->second;
	}
	void load_textfile(string filename){
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
				for(auto &word_str: word_strs){
					if(word_str.size() == 0){
						continue;
					}
					Word* word = new Word();
					word->token_id = string_to_token_id(word_str);
					word->phrase_id = -1;
					words.push_back(word);
				}
				_dataset.push_back(words);
			}
		}
		c_printf("[*]%s", (boost::format("%sを読み込みました. (%d行)") % filename.c_str() % _dataset.size()).str().c_str());
	}
	void perform_gibbs_sampling(){
		for(int epoch = 0;epoch < _max_epoch;epoch++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
		}
	}
	void set_num_phrases(int number){
		_hmm->set_num_phrases(number);
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
	.def("perform_gibbs_sampling", &PyBayesianHMM::perform_gibbs_sampling)
	.def("set_num_phrases", &PyBayesianHMM::set_num_phrases)
	.def("set_num_words", &PyBayesianHMM::set_num_words)
	.def("set_max_epoch", &PyBayesianHMM::set_max_epoch)
	.def("load_textfile", &PyBayesianHMM::load_textfile)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}