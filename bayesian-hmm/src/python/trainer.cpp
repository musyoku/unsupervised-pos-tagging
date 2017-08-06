#include <cassert>
#include <iostream>
#include "../bhmm/sampler.h"
#include "../bhmm/utils.h"
#include "trainer.h"

namespace bhmm {
	struct value_comparator {
		bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
			return a.second > b.second;
		}   
	};
	Trainer::Trainer(Dataset* dataset, Model* model, Dictionary* dict, boost::python::list py_Wt){
		// 日本語周り
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype);
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_model = model;
		std::vector<int> Wt = utils::vector_from_list<int>(py_Wt);
		_model->_hmm->initialize_with_training_corpus(dataset->_word_sequences_train, Wt);
		_dict = dict;
		_dataset = dataset;
	}
	void Trainer::perform_gibbs_sampling(){
		std::vector<std::vector<Word*>> &dataset = _dataset->_word_sequences_train;
		if(_rand_indices.size() != dataset.size()){
			_rand_indices.clear();
			for(int data_index = 0;data_index < dataset.size();data_index++){
				_rand_indices.push_back(data_index);
			}
		}
		shuffle(_rand_indices.begin(), _rand_indices.end(), sampler::mt);	// データをシャッフル
		for(int n = 0;n < dataset.size();n++){
			if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
				return;
			}
			int data_index = _rand_indices[n];
			std::vector<Word*> &word_vec = dataset[data_index];
			_model->_hmm->perform_gibbs_sampling_with_words(word_vec);
		}
	}
	void Trainer::update_hyperparameters(){
		_model->_hmm->sample_new_alpha(_dataset->_word_sequences_train);
		_model->_hmm->sample_new_beta(_dataset->_word_sequences_train);
	}
	boost::python::list Trainer::get_all_words_for_each_tag(int threshold){
		std::vector<boost::python::list> result;
		for(int tag = 0;tag < _model->_hmm->_num_tags;tag++){
			std::vector<boost::python::tuple> words;
			std::unordered_map<int, int> &word_counts = _model->_hmm->_tag_word_counts[tag];
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
	void Trainer::show_typical_words_for_each_tag(int number_to_show_for_each_tag){
		for(int tag = 0;tag < _model->_hmm->_num_tags;tag++){
			std::unordered_map<int, int> &word_counts = _model->_hmm->_tag_word_counts[tag];
			int n = 0;
			std::cout << "tag " << tag << std::endl;
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
}