#include <boost/python.hpp>
#include <boost/format.hpp>
#include <string>
#include <functional>
#include <fstream>
#include <cassert>
#include "util.h"
using namespace boost;

class PyBayesianHMM{
private:
	Word2Vec* _word2vec;
	std::hash<wstring> _hash_func;
	unordered_map<id, wstring> _dictionary;
	wstring _refference_word;
	bool _show_simliar_words_on_each_epoch;
	int _vector_length;
	int _max_epoch;
	int _num_updates_per_epoch;
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
		_word2vec = new Word2Vec();

		_vector_length = vector_length;		// 単語ベクトルの長さ
		_max_epoch = 1000;
		_num_updates_per_epoch = 500;
		_show_simliar_words_on_each_epoch = false; 
	}
	id string_to_token_id(wstring str){
		return (id)_hash_func(str);
	}
	void load_textfile(string filename){
		wifstream ifs(filename.c_str());
		wstring str;
		if (ifs.fail()){
			c_printf("[R]%s [*]%s", "エラー", (boost::format("%sを開けません.") % filename.c_str()).str().c_str());
			exit(1);
		}
		id bos_id = string_to_token_id(L"<bos>");
		id eos_id = string_to_token_id(L"<eos>");
		id unk_id = string_to_token_id(L"<unk>");
		while (getline(ifs, str) && !str.empty()){
			vector<wstring> words = split_word_by(str, ' ');	// スペースで分割
			if(words.size() > 0){
				vector<id> token_ids = {bos_id};
				for(auto word: words){
					if(word.size() == 0){
						continue;
					}
					id token_id = string_to_token_id(word);
					token_ids.push_back(token_id);
					_dictionary[token_id] = word;
				}
				token_ids.push_back(eos_id);		// bosとeosで挟む
				_word2vec->add_sentense(token_ids);	// 文を追加
			}
		}
		c_printf("[*]%s", (boost::format("%sを読み込みました. (%d行)") % filename.c_str() % _word2vec->_dataset.size()).str().c_str());
		_word2vec->generate_vectors_from_dataset();	// 読み込んだデータから単語ベクトルを生成
		_word2vec->generate_bigram_table();			// バイグラムテーブルの生成
	}
	void train(int num_threads){
		c_printf("[*]%s\n", (boost::format("学習を開始します. (スレッド数: %d)") % num_threads).str().c_str());
		for(int epoch = 0;epoch < _max_epoch;epoch++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			_word2vec->train(num_threads, _num_updates_per_epoch);
			c_printf("[n]%s\n", (boost::format("Epoch %d/%d - loss %lf") % epoch % _max_epoch % (_word2vec->_sum_loss / num_threads)).str().c_str());

			if(_show_simliar_words_on_each_epoch){
				set<pair<wstring, vectype>, compare> ranking;
				id token_id = string_to_token_id(_refference_word);
				vectype* vec = _word2vec->get_vector_for_token_id(token_id);
				assert(vec != NULL);
				for(auto &elem: _dictionary){
					id target_id = elem.first;
					vectype* _vec = _word2vec->get_vector_for_token_id(target_id);
					assert(_vec != NULL);
					vectype similarity = _word2vec->compute_cosine_similarity(vec, _vec);
					if(similarity > 0.7){
						pair<wstring, vectype> rank(elem.second, similarity);
						ranking.insert(rank);
					}
				}
			    for(auto itr = ranking.begin(); itr != ranking.end(); ++itr) {
			    	wcout << itr->first << L": " << itr->second << endl;
			    }
			}
		}
	}
	void set_learning_rate(vectype lr){
		_word2vec->set_learning_rate(lr);
	}
	void set_negative_marginal_window_size(int size){
		_word2vec->set_negative_marginal_window_size(size);
	}
	void set_num_positive_samples(int size){
		_word2vec->set_num_positive_samples(size);
	}
	void set_num_negative_samples(int size){
		_word2vec->set_num_negative_samples(size);
	}
	void set_num_updates_per_epoch(int updates){
		_num_updates_per_epoch = updates;
	}
	void set_max_epoch(int epoch){
		_max_epoch = epoch;
	}
	void show_similar_words_to(wstring word){
		_show_simliar_words_on_each_epoch = true;
		_refference_word = word;
	}
};

BOOST_PYTHON_MODULE(model){
	python::class_<PyBayesianHMM>("word2vec", python::init<int>())
	.def("string_to_token_id", &PyBayesianHMM::string_to_token_id)
	.def("train", &PyBayesianHMM::train)
	.def("set_negative_marginal_window_size", &PyBayesianHMM::set_negative_marginal_window_size)
	.def("set_learning_rate", &PyBayesianHMM::set_learning_rate)
	.def("set_num_positive_samples", &PyBayesianHMM::set_num_positive_samples)
	.def("set_num_negative_samples", &PyBayesianHMM::set_num_negative_samples)
	.def("set_num_updates_per_epoch", &PyBayesianHMM::set_num_updates_per_epoch)
	.def("set_max_epoch", &PyBayesianHMM::set_max_epoch)
	.def("show_similar_words_to", &PyBayesianHMM::show_similar_words_to)
	.def("load_textfile", &PyBayesianHMM::load_textfile);
}