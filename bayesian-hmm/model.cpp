#include <boost/python.hpp>
#include <boost/format.hpp>
#include <string>
#include <set>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "core/bhmm.h"
#include "core/util.h"
using namespace std;
using namespace boost;

struct value_comparator {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		if (a.second > b.second) return true;
		return false;
	}   
};

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
	PyBayesianHMM(){
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
		c_printf("[*]%s\n", (boost::format("%sを読み込みました.") % filename.c_str()).str().c_str());
	}
	void initialize(){
		_hmm->init_ngram_counts(_dataset);
	}
	void perform_gibbs_sampling(){
		vector<int> rand_indices;
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			rand_indices.push_back(data_index);
		}
		_hmm->set_temperature(1);
		for(int epoch = 1;epoch <= _max_epoch;epoch++){
			shuffle(rand_indices.begin(), rand_indices.end(), Sampler::mt);	// データをシャッフル
			for(int n = 0;n < _dataset.size();n++){
				if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
					return;
				}
				int data_index = rand_indices[n];
				vector<Word*> &line = _dataset[data_index];
				_hmm->perform_gibbs_sampling_with_line(line);
			}
			show_progress(epoch, _max_epoch);
			if(epoch % 100 == 0){
				show_random_line(10);
				_hmm->dump_word_types();
				show_typical_words_for_each_tag(20);
			}
			// _hmm->anneal_temperature(0.99989);
		}
	}
	void show_random_line(int num_to_show, bool show_most_co_occurring_tag = true){
		for(int n = 0;n < num_to_show;n++){
			int data_index = Sampler::uniform_int(0, _dataset.size() - 1);
			vector<Word*> &line = _dataset[data_index];
			for(int pos = 2;pos < line.size() - 2;pos++){
				Word* word = line[pos];
				int tag_id = word->tag_id;
				if(show_most_co_occurring_tag){
					tag_id = _hmm->get_most_co_occurring_tag(word->word_id);
				}
				wcout << _dictionary[word->word_id] << L"/" << tag_id << L" ";
			}
			wcout << endl;
		}
	}
	void show_typical_words_for_each_tag(int number_to_show_for_each_tag){
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			map<int, int> &word_counts = _hmm->_tag_word_counts[tag];
			int n = 0;
			wcout << L"tag " << tag << L":" << endl;
			wcout << L"	";
			multiset<pair<int, int>, value_comparator> ranking;
			for(auto elem: word_counts){
				ranking.insert(std::make_pair(elem.first, elem.second));
			}
			for(auto elem: ranking){
				wstring word = _dictionary[elem.first];
				wcout << word << L"/" << elem.second << L", ";
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
				if(elem.second < 10){
					break;
				}
			}
			cout << endl;
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
	python::class_<PyBayesianHMM>("bayesian_hmm")
	.def("string_to_word_id", &PyBayesianHMM::string_to_word_id)
	.def("perform_gibbs_sampling", &PyBayesianHMM::perform_gibbs_sampling)
	.def("initialize", &PyBayesianHMM::initialize)
	.def("set_num_tags", &PyBayesianHMM::set_num_tags)
	.def("set_max_epoch", &PyBayesianHMM::set_max_epoch)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}