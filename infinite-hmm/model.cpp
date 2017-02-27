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

class PyInfiniteHMM{
private:
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	unordered_map<int, int> _word_count;
	vector<vector<Word*>> _dataset;
	vector<int> _rand_indices;
	int _autoincrement;
	int _bos_id;
	int _eos_id;
	int _unk_id;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
	double _minimum_temperature;
public:
	InfiniteHMM* _hmm;
	PyInfiniteHMM(int initial_num_tags){
		// 日本語周り
		// ただのテンプレ
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);

		_hmm = new InfiniteHMM(initial_num_tags + 1);
		_bos_id = 0;
		_dictionary[_bos_id] = L"<bos>";
		_eos_id = 1;
		_dictionary[_eos_id] = L"<eos>";
		_unk_id = 2;
		_dictionary[_unk_id] = L"<unk>";
		_autoincrement = _unk_id + 1;

		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;

		_minimum_temperature = 0.08;
	}
	int add_string(wstring word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			_dictionary[_autoincrement] = word;
			_dictionary_inv[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	int string_to_word_id(wstring word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			return _unk_id;
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
			// Word* bos = new Word();
			// bos->word_id = -1;
			// bos->tag_id = BOP;
			// words.push_back(bos);
			// _word_count[_bos_id] += 1;

			for(auto &word_str: word_strs){
				if(word_str.size() == 0){
					continue;
				}
				Word* word = new Word();
				word->word_id = add_string(word_str);
				words.push_back(word);
				_word_count[word->word_id] += 1;
			}

			Word* eos = new Word();
			eos->word_id = _eos_id;
			words.push_back(eos);
			_word_count[_eos_id] += 1;

			// Word* eop = new Word();
			// eop->word_id = -1;
			// eop->tag_id = EOP;
			// words.push_back(eop);

			// 訓練データに追加
			_dataset.push_back(words);
		}
	}
	int get_count_for_word(int word_id){
		auto itr = _word_count.find(word_id);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
	int get_num_tags(){
		return _hmm->get_num_tags();
	}
	void mark_low_frequency_words_as_unknown(int threshold = 1){
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			vector<Word*> &line = _dataset[data_index];
			for(auto word = line.begin(), end = line.end();word != end;word++){
				int word_id = (*word)->word_id;
				int count = get_count_for_word(word_id);
				if(count <= threshold){
					(*word)->word_id = _unk_id;
				}
			}
		}
	}
	void initialize(){
		_hmm->initialize(_dataset);
	}
	bool load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/ihmm.dict";
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
		std::ofstream ofs(dirname + "/ihmm.dict");
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		oarchive << _dictionary_inv;
		oarchive << _autoincrement;
		ofs.close();
		return _hmm->save(dirname);
	}
	int argmax_Ptag_context_word(int context_tag_id, int word_id){
		return _hmm->argmax_Ptag_context_word(context_tag_id, word_id);
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
	void perform_beam_sampling(){
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
			_hmm->perform_beam_sampling_with_line(line);
		}
	}
	void set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void anneal_temperature(double multiplier){
		if(_hmm->_temperature > _minimum_temperature){
			_hmm->_temperature *= multiplier;
		}
	}
	void show_temperature(){
		c_printf("[*]%s: %lf\n", "temperature", _hmm->_temperature);
	}
	void show_log_Pdata(){
		double log_p = 0;
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			vector<Word*> &line = _dataset[data_index];
			log_p += _hmm->compute_log_Pdata(line);
		}
		c_printf("[*]%s: %lf\n", "log_Pdata", log_p);
	}
	void show_typical_words_for_each_tag(int number_to_show_for_each_tag){
		auto pair = std::make_pair(0, 0);
		for(int tag = 0;tag < _hmm->_tag_unigram_count.size();tag++){
			if(_hmm->_tag_unigram_count[tag] == 0){
				continue;
			}
			multiset<std::pair<int, int>, value_comparator> ranking;
			unordered_map<int, Table*> &tables = _hmm->_tag_word_table[tag];
			int n = 0;
			c_printf("[*]%s\n", (boost::format("tag %d:") % tag).str().c_str());
			wcout << L"\t";
			for(const auto &elem: tables){
				pair.first = elem.first;
				pair.second = elem.second->_num_customers;
				ranking.insert(pair);
			}
			for(const auto &elem: ranking){
				wstring &word = _dictionary[elem.first];
				wcout << word << L"/" << elem.second << L", ";
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			wcout << endl;
			ranking.clear();
		}
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyInfiniteHMM>("ihmm", python::init<int>())
	.def("string_to_word_id", &PyInfiniteHMM::string_to_word_id)
	.def("add_string", &PyInfiniteHMM::add_string)
	.def("perform_gibbs_sampling", &PyInfiniteHMM::perform_gibbs_sampling)
	.def("perform_beam_sampling", &PyInfiniteHMM::perform_beam_sampling)
	.def("initialize", &PyInfiniteHMM::initialize)
	.def("set_temperature", &PyInfiniteHMM::set_temperature)
	.def("anneal_temperature", &PyInfiniteHMM::anneal_temperature)
	.def("load", &PyInfiniteHMM::load)
	.def("save", &PyInfiniteHMM::save)
	.def("add_line", &PyInfiniteHMM::add_line)
	.def("mark_low_frequency_words_as_unknown", &PyInfiniteHMM::mark_low_frequency_words_as_unknown)
	.def("show_typical_words_for_each_tag", &PyInfiniteHMM::show_typical_words_for_each_tag)
	.def("show_log_Pdata", &PyInfiniteHMM::show_log_Pdata)
	.def("show_temperature", &PyInfiniteHMM::show_temperature)
	.def("argmax_Ptag_context_word", &PyInfiniteHMM::argmax_Ptag_context_word)
	.def("get_num_tags", &PyInfiniteHMM::get_num_tags)
	.def("load_textfile", &PyInfiniteHMM::load_textfile);
}