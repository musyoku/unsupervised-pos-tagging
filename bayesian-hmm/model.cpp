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
	vector<int> _rand_indices;
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
		_hmm->initialize(_dataset);
	}
	bool load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/hmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _dictionary;
			ifs.close();
		}
		// 獲得された品詞と単語のリストを読み込み
		string tags_filename = dirname + "/hmm.tags";
		return _hmm->load(tags_filename);
	}
	bool save(string dirname){
		// 辞書を保存
		string dictionary_filename = dirname + "/hmm.dict";
		std::ofstream ofs(dictionary_filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		ofs.close();
		// 獲得された品詞と単語のリストを保存
		string tags_filename = dirname + "/hmm.tags";
		return _hmm->save(tags_filename);
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
	void sample_new_beta(){
		_hmm->sample_new_beta(_dataset);
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
	python::list get_all_words_for_each_tag(int threshold = 0){
		vector<python::list> result;
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			vector<python::tuple> words;
			unordered_map<int, int> &word_counts = _hmm->_tag_word_counts[tag];
			multiset<pair<int, int>, value_comparator> ranking;
			for(auto elem: word_counts){
				ranking.insert(std::make_pair(elem.first, elem.second));
			}
			for(auto elem: ranking){
				if(elem.second <= threshold){
					continue;
				}
				wstring word = _dictionary[elem.first];
				words.push_back(python::make_tuple(word, elem.second));
			}
			result.push_back(list_from_vector(words));
		}
		return list_from_vector(result);
	}
	void show_typical_words_for_each_tag(int number_to_show_for_each_tag){
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			unordered_map<int, int> &word_counts = _hmm->_tag_word_counts[tag];
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
	void set_alpha(double alpha){
		_hmm->_alpha = alpha;
	}
	void set_num_tags(int number){
		_hmm->_num_tags = number;
	}
	int get_num_tags(){
		return _hmm->_num_tags;
	}
	void set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void set_Wt(python::list Wt){
		int length = python::len(Wt);
		for(int tag = 0;tag < length;tag++){
			_hmm->set_Wt_for_tag(tag, python::extract<int>(Wt[tag]));
		}
	}
	void set_minimum_temperature(double temperature){
		_hmm->_minimum_temperature = temperature;
	}
	void anneal_temperature(double temperature){
		_hmm->anneal_temperature(temperature);
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyBayesianHMM>("bayesian_hmm")
	.def("string_to_word_id", &PyBayesianHMM::string_to_word_id)
	.def("perform_gibbs_sampling", &PyBayesianHMM::perform_gibbs_sampling)
	.def("initialize", &PyBayesianHMM::initialize)
	.def("load", &PyBayesianHMM::load)
	.def("save", &PyBayesianHMM::save)
	.def("get_num_tags", &PyBayesianHMM::get_num_tags)
	.def("get_all_words_for_each_tag", &PyBayesianHMM::get_all_words_for_each_tag)
	.def("set_temperature", &PyBayesianHMM::set_temperature)
	.def("set_num_tags", &PyBayesianHMM::set_num_tags)
	.def("set_minimum_temperature", &PyBayesianHMM::set_minimum_temperature)
	.def("set_Wt", &PyBayesianHMM::set_Wt)
	.def("sample_new_beta", &PyBayesianHMM::sample_new_beta)
	.def("anneal_temperature", &PyBayesianHMM::anneal_temperature)
	.def("show_typical_words_for_each_tag", &PyBayesianHMM::show_typical_words_for_each_tag)
	.def("show_random_line", &PyBayesianHMM::show_random_line)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}