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
	int _bos_id;
	int _eos_id;
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
		_bos_id = 0;
		_eos_id = 1;
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
				for(int n = 0;n < 2;n++){
					Word* bos = new Word();
					bos->word_id = _bos_id;
					bos->tag_id = -1;
					words.push_back(bos);
				}
				for(auto &word_str: word_strs){
					if(word_str.size() == 0){
						continue;
					}
					Word* word = new Word();
					word->word_id = string_to_word_id(word_str);
					word->tag_id = -1;
					words.push_back(word);
				}
				// <eos>も2つ追加しておくとt_{i+1}, t_{i+2}が常に存在するのでギブスサンプリング時に場合分けしなくてもいいかもしれない
				for(int n = 0;n < 2;n++){
					Word* eos = new Word();
					eos->word_id = _eos_id;
					eos->tag_id = -1;
					words.push_back(eos);
				}
				// 訓練データに追加
				_dataset.push_back(words);
			}
		}
		c_printf("[*]%s", (boost::format("%sを読み込みました. (%d行)") % filename.c_str() % _dataset.size()).str().c_str());
	}
	void initialize(){
		_hmm->init_ngram_counts(_dataset);
	}
	void perform_gibbs_sampling(){
		vector<int> rand_indices;
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			rand_indices.push_back(data_index);
		}
		for(int epoch = 0;epoch < _max_epoch;epoch++){
			shuffle(rand_indices.begin(), rand_indices.end(), Sampler::mt);	// データをシャッフル
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
		}
	}
	void set_num_tags(int number){
		_hmm->set_num_tags(number);
	}
	void set_max_epoch(int epoch){
		_max_epoch = epoch;
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyBayesianHMM>("word2vec", python::init<int>())
	.def("string_to_word_id", &PyBayesianHMM::string_to_word_id)
	.def("perform_gibbs_sampling", &PyBayesianHMM::perform_gibbs_sampling)
	.def("initialize", &PyBayesianHMM::initialize)
	.def("set_num_tags", &PyBayesianHMM::set_num_tags)
	.def("set_max_epoch", &PyBayesianHMM::set_max_epoch)
	.def("load_textfile", &PyBayesianHMM::load_textfile)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}