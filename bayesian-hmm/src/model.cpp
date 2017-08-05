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
#include <iostream>
#include <cassert>
#include "bhmm/hmm.h"
#include "bhmm/utils.h"
#include "dictionary.h"

struct value_comparator {
	bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
		return a.second > b.second;
	}   
};

class Trainer{
private:
	HMM* _hmm;
	Dictionary* _dict;
	std::unordered_map<int, int> _word_count;
	std::vector<std::vector<Word*>> _dataset;
	std::vector<int> _rand_indices;
	int _autoincrement;
	int _bos_id;
	int _eos_id;
	int _unk_id;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
public:
	Trainer(Dictionary* dict){
		// 日本語周り
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_hmm = new HMM();
		_dict = dict;
		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
	}
	void load_textfile(std::string filename){
		c_printf("[*]%s\n", (boost::format("%sを読み込んでいます ...") % filename.c_str()).str().c_str());
		std::wifstream ifs(filename.c_str());
		std::wstring line_str;
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
	void add_line(std::wstring line_str){
		std::vector<std::wstring> word_strs;
		utils::split_word_by(line_str, L' ', word_strs);	// スペースで分割
		int num_words = word_strs.size();
		if(num_words > _max_num_words_in_line){
			_max_num_words_in_line = num_words;
		}
		if(num_words < _max_num_words_in_line || _min_num_words_in_line == -1){
			_min_num_words_in_line = num_words;
		}
		if(word_strs.size() > 0){
			std::vector<Word*> words;
			// <bos>
			for(int n = 0;n < 2;n++){
				Word* bos = new Word();
				bos->word_id = _bos_id;
				bos->tag_id = 0;
				words.push_back(bos);
				_word_count[_bos_id] += 1;
			}
			for(auto &word_str: word_strs){
				if(word_str.size() == 0){
					continue;
				}
				Word* word = new Word();
				word->word_id = _dict->add_word_string(word_str);
				word->tag_id = 0;
				words.push_back(word);
				_word_count[word->word_id] += 1;
			}
			// <eos>も2つ追加しておくとt_{i+1}, t_{i+2}が常に存在するのでギブスサンプリング時に場合分けしなくてもいいかもしれない
			for(int n = 0;n < 2;n++){
				Word* eos = new Word();
				eos->word_id = _eos_id;
				eos->tag_id = 0;
				words.push_back(eos);
				_word_count[_eos_id] += 1;
			}
			// 訓練データに追加
			_dataset.push_back(words);
		}
	}
	void initialize(){
		_hmm->initialize(_dataset);
	}
	void mark_low_frequency_words_as_unknown(int threshold = 1){
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			std::vector<Word*> &line = _dataset[data_index];
			for(auto word = line.begin(), end = line.end();word != end;word++){
				int word_id = (*word)->word_id;
				int count = get_count_for_word(word_id);
				if(count <= threshold){
					(*word)->word_id = _unk_id;
				}
			}
		}
	}
	int get_count_for_word(int word_id){
		auto itr = _word_count.find(word_id);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
	bool load(std::string dirname){
		// 辞書を読み込み
		std::string dictionary_filename = dirname + "/hmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _autoincrement;
			ifs.close();
		}
		return _hmm->load(dirname);
	}
	bool save(std::string dirname){
		// 辞書を保存
		std::ofstream ofs(dirname + "/hmm.dict");
		boost::archive::binary_oarchive oarchive(ofs);
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
		shuffle(_rand_indices.begin(), _rand_indices.end(), sampler::mt);	// データをシャッフル
		for(int n = 0;n < _dataset.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			std::vector<Word*> &line = _dataset[data_index];
			_hmm->perform_gibbs_sampling_with_line(line);
		}
	}
	int sample_tag_from_Pt_w(int ti_2, int ti_1, int wi){
		return _hmm->sample_tag_from_Pt_w(ti_2, ti_1, wi);
	}
	int argmax_tag_from_Pt_w(int ti_2, int ti_1, int wi){
		return _hmm->argmax_tag_from_Pt_w(ti_2, ti_1, wi);
	}
	void sample_new_alpha(){
		_hmm->sample_new_alpha(_dataset);
	}
	void show_alpha(){
		std::cout << (boost::format("alpha <- %e") % _hmm->_alpha).str() << std::endl;
	}
	void sample_new_beta(){
		_hmm->sample_new_beta(_dataset);
	}
	void show_beta(){
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			std::cout << (boost::format("beta[%d] <- %e") % tag % _hmm->_beta[tag]).str() << std::endl;
		}
	}
	void show_random_line(int num_to_show, bool show_most_co_occurring_tag = true){
		for(int n = 0;n < num_to_show;n++){
			int data_index = sampler::uniform_int(0, _dataset.size() - 1);
			std::vector<Word*> &line = _dataset[data_index];
			for(int pos = 2;pos < line.size() - 2;pos++){
				Word* word = line[pos];
				int tag_id = word->tag_id;
				if(show_most_co_occurring_tag){
					tag_id = _hmm->get_most_co_occurring_tag(word->word_id);
				}
				std::wcout << _dict->word_id_to_string(word->word_id) << L"/" << tag_id << L" ";
			}
			std::wcout << std::endl;
		}
	}
	boost::python::list get_all_words_for_each_tag(int threshold = 0){
		std::vector<boost::python::list> result;
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			std::vector<boost::python::tuple> words;
			std::unordered_map<int, int> &word_counts = _hmm->_tag_word_counts[tag];
			std::multiset<std::pair<int, int>, value_comparator> ranking;
			for(auto elem: word_counts){
				ranking.insert(std::make_pair(elem.first, elem.second));
			}
			for(auto elem: ranking){
				if(elem.second <= threshold){
					continue;
				}
				std::wstring word = _dict->word_id_to_string(elem.first);
				words.push_back(boost::python::make_tuple(word, elem.second));
			}
			result.push_back(utils::list_from_vector(words));
		}
		return utils::list_from_vector(result);
	}
	void show_typical_words_for_each_tag(int number_to_show_for_each_tag){
		for(int tag = 0;tag < _hmm->_num_tags;tag++){
			std::unordered_map<int, int> &word_counts = _hmm->_tag_word_counts[tag];
			int n = 0;
			c_printf("[*]%s\n", (boost::format("tag %d:") % tag).str().c_str());
			std::wcout << L"\t";
			std::multiset<std::pair<int, int>, value_comparator> ranking;
			for(auto elem: word_counts){
				ranking.insert(std::make_pair(elem.first, elem.second));
			}
			for(auto elem: ranking){
				std::wstring word = _dict->word_id_to_string(elem.first);
				std::wcout << word << L"/" << elem.second << L", ";
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			std::wcout << std::endl;
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
	double get_temperature(){
		return _hmm->_temperature;
	}
	void set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void set_Wt(boost::python::list Wt){
		int length = boost::python::len(Wt);
		for(int tag = 0;tag < length;tag++){
			_hmm->set_Wt_for_tag(tag, boost::python::extract<int>(Wt[tag]));
		}
	}
	void set_minimum_temperature(double temperature){
		_hmm->_minimum_temperature = temperature;
	}
	void anneal_temperature(double temperature){
		_hmm->anneal_temperature(temperature);
	}
	int get_max_num_words_in_line(){
		return _max_num_words_in_line;
	}
	int get_min_num_words_in_line(){
		return _min_num_words_in_line;
	}
};

BOOST_PYTHON_MODULE(bhmm){
	boost::python::class_<Trainer>("trainer", boost::python::init<Dictionary*>())
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling)
	.def("initialize", &Trainer::initialize)
	.def("mark_low_frequency_words_as_unknown", &Trainer::mark_low_frequency_words_as_unknown)
	.def("load", &Trainer::load)
	.def("save", &Trainer::save)
	.def("get_num_tags", &Trainer::get_num_tags)
	.def("get_all_words_for_each_tag", &Trainer::get_all_words_for_each_tag)
	.def("get_temperature", &Trainer::get_temperature)
	.def("get_max_num_words_in_line", &Trainer::get_max_num_words_in_line)
	.def("get_min_num_words_in_line", &Trainer::get_min_num_words_in_line)
	.def("set_temperature", &Trainer::set_temperature)
	.def("set_num_tags", &Trainer::set_num_tags)
	.def("set_minimum_temperature", &Trainer::set_minimum_temperature)
	.def("set_Wt", &Trainer::set_Wt)
	.def("set_alpha", &Trainer::set_alpha)
	.def("add_line", &Trainer::add_line)
	.def("sample_new_alpha", &Trainer::sample_new_alpha)
	.def("sample_new_beta", &Trainer::sample_new_beta)
	.def("sample_tag_from_Pt_w", &Trainer::sample_tag_from_Pt_w)
	.def("argmax_tag_from_Pt_w", &Trainer::argmax_tag_from_Pt_w)
	.def("anneal_temperature", &Trainer::anneal_temperature)
	.def("show_typical_words_for_each_tag", &Trainer::show_typical_words_for_each_tag)
	.def("show_random_line", &Trainer::show_random_line)
	.def("show_alpha", &Trainer::show_alpha)
	.def("show_beta", &Trainer::show_beta)
	.def("load_textfile", &Trainer::load_textfile);
}