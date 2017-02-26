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
#include "core/hpylm.hpp"
#include "core/ithmm.h"
#include "core/util.h"
using namespace std;
using namespace boost;

class PyInfiniteTreeHMM{
private:
	unordered_map<id, wstring> _dictionary;
	unordered_map<wstring, id> _dictionary_inv;
	unordered_map<id, int> _word_count;
	vector<vector<Word*>> _dataset;
	vector<int> _rand_indices;
	id _autoincrement;
	id _bos_id;
	id _eos_id;
	id _unk_id;
	int _max_num_words_in_line;
	int _min_num_words_in_line;
public:
	iTHMM* _ithmm;
	PyInfiniteTreeHMM(){
		// 日本語周り
		// ただのテンプレ
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		ios_base::sync_with_stdio(false);
		locale default_loc("ja_JP.UTF-8");
		locale::global(default_loc);
		locale ctype_default(locale::classic(), default_loc, locale::ctype); //※
		wcout.imbue(ctype_default);
		wcin.imbue(ctype_default);

		_ithmm = new iTHMM();
		_bos_id = 0;
		_dictionary[_bos_id] = L"<bos>";
		_eos_id = 1;
		_dictionary[_eos_id] = L"<eos>";
		_unk_id = 2;
		_dictionary[_unk_id] = L"<unk>";
		_autoincrement = _unk_id + 1;

		_max_num_words_in_line = -1;
		_min_num_words_in_line = -1;
	}
	id add_string(wstring word){
		auto itr = _dictionary_inv.find(word);
		if(itr == _dictionary_inv.end()){
			_dictionary[_autoincrement] = word;
			_dictionary_inv[word] = _autoincrement;
			_autoincrement++;
			return _autoincrement - 1;
		}
		return itr->second;
	}
	id string_to_word_id(wstring word){
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

			for(auto word_str: word_strs){
				if(word_str.size() == 0){
					continue;
				}
				Word* word = new Word();
				word->id = add_string(word_str);
				word->state = NULL;
				words.push_back(word);
				_word_count[word->id] += 1;
			}

			Word* eos = new Word();
			eos->id = _eos_id;
			eos->state = NULL;
			words.push_back(eos);
			_word_count[_eos_id] += 1;

			// 訓練データに追加
			_dataset.push_back(words);
		}
	}
	int get_num_words(){
		return _word_count.size();
	}
	int get_count_for_word(id word_id){
		auto itr = _word_count.find(word_id);
		if(itr == _word_count.end()){
			return 0;
		}
		return itr->second;
	}
	void mark_low_frequency_words_as_unknown(int threshold = 1){
		for(int data_index = 0;data_index < _dataset.size();data_index++){
			vector<Word*> &line = _dataset[data_index];
			for(auto word = line.begin(), end = line.end();word != end;word++){
				id word_id = (*word)->id;
				int count = get_count_for_word(word_id);
				if(count <= threshold){
					(*word)->id = _unk_id;
				}
			}
		}
	}
	void compile(){
		_ithmm->set_word_g0(1.0 / _word_count.size());
		_ithmm->initialize_data(_dataset);
	}
	void remove_all_data(){
		_ithmm->remove_all_data(_dataset);
		_ithmm->delete_invalid_children_on_structure_tssb(_ithmm->_structure_tssb);
	}
	bool load(string dirname){
		// 辞書を読み込み
		string dictionary_filename = dirname + "/ithmm.dict";
		std::ifstream ifs(dictionary_filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> _dictionary;
			iarchive >> _dictionary_inv;
			iarchive >> _autoincrement;
			ifs.close();
		}
		return _ithmm->load(dirname);
	}
	bool save(string dirname){
		// 辞書を保存
		std::ofstream ofs(dirname + "/ithmm.dict");
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << _dictionary;
		oarchive << _dictionary_inv;
		oarchive << _autoincrement;
		ofs.close();
		return _ithmm->save(dirname);
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
			_ithmm->perform_gibbs_sampling_line(line);
		}
	}
	void update_hyperparameters(){
		_ithmm->sample_hpylm_hyperparameters();
	}
	void show_typical_words_for_each_tag(int number_to_show_for_each_tag, bool show_probability = true){
		auto pair = std::make_pair(0, 0);
		vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int count = 0;
			string indices = node->_dump_indices();
			wstring tab = L"";
			for(int i = 0;i < node->_depth_v;i++){
				tab += L"	";
			}
			wcout << tab;
			c_printf("[*]%s\n", (boost::format("[%s]") % indices.c_str()).str().c_str());
			wcout << tab;
			for(const auto &elem: ranking){
				wstring &word = _dictionary[elem.first];
				wcout << word << L" ";
				if(show_probability){
					wcout << elem.second << L" ";
				} 
				count++;
				if(count > number_to_show_for_each_tag){
					break;
				}
			}
			wcout << endl;
			ranking.clear();
		}
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyInfiniteTreeHMM>("ithmm")
	.def("string_to_word_id", &PyInfiniteTreeHMM::string_to_word_id)
	.def("add_string", &PyInfiniteTreeHMM::add_string)
	.def("perform_gibbs_sampling", &PyInfiniteTreeHMM::perform_gibbs_sampling)
	.def("compile", &PyInfiniteTreeHMM::compile)
	.def("load", &PyInfiniteTreeHMM::load)
	.def("save", &PyInfiniteTreeHMM::save)
	.def("add_line", &PyInfiniteTreeHMM::add_line)
	.def("mark_low_frequency_words_as_unknown", &PyInfiniteTreeHMM::mark_low_frequency_words_as_unknown)
	.def("load_textfile", &PyInfiniteTreeHMM::load_textfile)
	.def("update_hyperparameters", &PyInfiniteTreeHMM::update_hyperparameters)
	.def("get_num_words", &PyInfiniteTreeHMM::get_num_words)
	.def("get_count_for_word", &PyInfiniteTreeHMM::get_count_for_word)
	.def("remove_all_data", &PyInfiniteTreeHMM::remove_all_data);
}