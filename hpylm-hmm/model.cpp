#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
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
#include "core/lattice.h"
#include "core/util.h"
using namespace std;
using namespace boost;

struct value_comparator {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.second > b.second;
	}   
};

class PyHpylmHMM{
private:
	HPYLM** _word_hpylm_for_tag;
	HPYLM* _pos_hpylm;
	Lattice* _lattice;
	unordered_map<int, wstring> _dictionary;
	unordered_map<wstring, int> _dictionary_inv;
	vector<vector<Word*>> _dataset;
	vector<int> _rand_indices;
	set<int> _types_of_words;
	int _autoincrement;
	int _num_tags;
	int _max_num_words_in_sentence;	// 1文あたりの最大単語数
	bool _is_ready;
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
		_pos_hpylm = new HPYLM(3);
		_word_hpylm_for_tag = new HPYLM*[num_tags];
		for(int tag = 0;tag < num_tags;tag++){
			_word_hpylm_for_tag[tag] = new HPYLM(3);
		}
		_lattice = NULL;
		_num_tags = num_tags;
		_dictionary[BEGIN_OF_SENTENSE] = L"<bos>";
		_dictionary[END_OF_SENTENSE] = L"<eos>";
		_autoincrement = END_OF_SENTENSE + 1;
		_max_num_words_in_sentence = 0;
		_is_ready = false;
	}
	~PyHpylmHMM(){
		delete _pos_hpylm;
		for(int tag = 0;tag < _num_tags;tag++){
			delete _word_hpylm_for_tag[tag];
		}
		delete[] _word_hpylm_for_tag;
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
		int prev_dataset_size = _dataset.size();
		while (getline(ifs, line_str) && !line_str.empty()){
			vector<wstring> word_strs = split_word_by(line_str, ' ');	// スペースで分割
			if(word_strs.size() > 0){
				vector<Word*> words;
				// <bos>
				// 3-gramなので2つ
				for(int n = 0;n < 2;n++){
					Word* bos = new Word();
					bos->word_id = BEGIN_OF_SENTENSE;
					bos->tag_id = BEGIN_OF_POS;
					words.push_back(bos);
				}
				for(auto &word_str: word_strs){
					if(word_str.size() == 0){
						continue;
					}
					Word* word = new Word();
					word->word_id = string_to_word_id(word_str);
					words.push_back(word);
					_types_of_words.insert(word->word_id);
				}
				Word* eos = new Word();
				eos->word_id = END_OF_SENTENSE;
				eos->tag_id = END_OF_POS;
				words.push_back(eos);
				// 単語数を記録しておく
				if(words.size() > _max_num_words_in_sentence){
					_max_num_words_in_sentence = words.size();
				}
				// 訓練データに追加
				_dataset.push_back(words);
			}
		}
		c_printf("[*]%s\n", (boost::format("%sを読み込みました. (%d行)") % filename.c_str() % (_dataset.size() - prev_dataset_size)).str().c_str());
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
		_pos_hpylm->load(dirname + "/pos.hpylm");
		for(int tag = 0;tag < _num_tags;tag++){
			_word_hpylm_for_tag[tag]->load((boost::format("%s/word.%d.hpylm") % dirname.c_str() % tag).str());
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
		_pos_hpylm->save(dirname + "/pos.hpylm");
		for(int tag = 0;tag < _num_tags;tag++){
			_word_hpylm_for_tag[tag]->save((boost::format("%s/word.%d.hpylm") % dirname.c_str() % tag).str());
		}
	}
	void prepare_for_training(){
		if(_lattice == NULL){
			_lattice = new Lattice(_max_num_words_in_sentence, _num_tags, _pos_hpylm, _word_hpylm_for_tag);
		}
		if(_rand_indices.size() != _dataset.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < _dataset.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		// 基底分布G0を設定
		assert(_pos_hpylm != NULL);
		_pos_hpylm->set_g0(1.0 / _num_tags);
		for(int tag = 0;tag < _num_tags;tag++){
			HPYLM* word_hpylm = _word_hpylm_for_tag[tag];
			assert(word_hpylm != NULL);
			word_hpylm->set_g0(1.0 / get_num_types_of_word());
		}
		_is_ready = true;
	}
	void perform_gibbs_sampling(){
		assert(_is_ready);
		assert(_rand_indices.size() == _dataset.size());
		shuffle(_rand_indices.begin(), _rand_indices.end(), Sampler::mt);	// データをシャッフル
		for(int n = 0;n < _dataset.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			vector<Word*> &sentence = _dataset[data_index];
			_lattice->perform_blocked_gibbs_sampling(sentence, false);	// argmax=false
			show_progress(n, _dataset.size());
		}
	}
	int get_num_types_of_word(){
		return _types_of_words.size();
	}
	int get_num_tags(){
		return _num_tags;
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyHpylmHMM>("hpylm_hmm", python::init<int>())
	.def("string_to_word_id", &PyHpylmHMM::string_to_word_id)
	.def("prepare_for_training", &PyHpylmHMM::prepare_for_training)
	.def("perform_gibbs_sampling", &PyHpylmHMM::perform_gibbs_sampling)
	.def("load", &PyHpylmHMM::load)
	.def("save", &PyHpylmHMM::save)
	.def("get_num_tags", &PyHpylmHMM::get_num_tags)
	.def("load_textfile", &PyHpylmHMM::load_textfile);
}