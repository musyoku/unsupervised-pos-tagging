#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/split_free.hpp>
#include <iostream>
#include <fstream>
#include "hmm.h"
#include "sampler.h"
#include "utils.h"

namespace bhmm {
	HMM::HMM(){
		_trigram_counts = NULL;
		_bigram_counts = NULL;
		_unigram_counts = NULL;
		_sampling_table = NULL;
		_Wt = NULL;
		_num_tags = -1;
		_num_words = -1;
		_alpha = 1;
		_beta = NULL;
		_temperature = 1;
		_minimum_temperature = 1;
		_allocated = false;
	}
	HMM::~HMM(){
		if(_allocated == false){
			return;
		}
		for(int tri_tag = 0;tri_tag <= _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag <= _num_tags;bi_tag++){
				delete[] _trigram_counts[tri_tag][bi_tag];
			}
			delete[] _trigram_counts[tri_tag];
		}
		delete[] _trigram_counts;
		for(int bi_tag = 0;bi_tag <= _num_tags;bi_tag++){
			delete[] _bigram_counts[bi_tag];
		}
		delete[] _bigram_counts;
		delete[] _unigram_counts;
		delete[] _sampling_table;
		delete[] _beta;
		delete[] _Wt;
	}
	void HMM::anneal_temperature(double multiplier){
		if(_temperature > _minimum_temperature){
			_temperature *= multiplier;
		}
	}
	void HMM::initialize_with_training_corpus(std::vector<std::vector<Word*>> &dataset, std::vector<int> &Wt){
		int length = Wt.size();
		assert(length > 0);
		_num_tags = length;
		alloc_count_tables(_num_tags);
		init_ngram_counts_with_corpus(dataset);
		for(int tag = 1;tag <= length;tag++){
			set_Wt_for_tag(tag, Wt[tag - 1]);
		}
	}
	void HMM::alloc_count_tables(int num_tags){
		assert(_allocated == false);
		assert(num_tags > 0);
		// 各タグの可能な単語数
		_Wt = new int[num_tags + 1];
		for(int tag = 0;tag <= num_tags;tag++){
			_Wt[tag] = 0;
		}
		// Betaの初期化
		// 初期値は1
		_beta = new double[num_tags + 1];
		for(int tag = 0;tag <= num_tags;tag++){
			_beta[tag] = 1;
		}
		// 3-gram		
		_trigram_counts = new int**[num_tags + 1];
		for(int tri_tag = 0;tri_tag <= num_tags;tri_tag++){
			_trigram_counts[tri_tag] = new int*[num_tags + 1];
			for(int bi_tag = 0;bi_tag <= num_tags;bi_tag++){
				_trigram_counts[tri_tag][bi_tag] = new int[num_tags + 1];
				for(int uni_tag = 0;uni_tag <= num_tags;uni_tag++){
					_trigram_counts[tri_tag][bi_tag][uni_tag] = 0;
				}
			}
		}
		// 2-gram
		_bigram_counts = new int*[num_tags + 1];
		for(int bi_tag = 0;bi_tag <= num_tags;bi_tag++){
			_bigram_counts[bi_tag] = new int[num_tags + 1];
			for(int uni_tag = 0;uni_tag <= num_tags;uni_tag++){
				_bigram_counts[bi_tag][uni_tag] = 0;
			}
		}
		// 1-gram
		_unigram_counts = new int[num_tags + 1];
		for(int tag = 0;tag <= num_tags;tag++){
			_unigram_counts[tag] = 0;
		}
		_allocated = true;
	}
	void HMM::init_ngram_counts_with_corpus(std::vector<std::vector<Word*>> &dataset){
		// 最初は品詞をランダムに割り当てる
		assert(_num_tags != -1);
		std::set<int> word_set;
		// std::unordered_map<int, int> tag_for_word;
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &word_vec = dataset[data_index];
			for(int pos = 2;pos < word_vec.size();pos++){	// 3-gramなので3番目から.
				Word* word = word_vec[pos];
				// auto itr = tag_for_word.find(word->_id);
				// if(itr == tag_for_word.end()){
					int state = sampler::uniform_int(1, _num_tags);
					assert(1 <= state && state <= _num_tags);
					word->_state = state;
					// tag_for_word[word->_id] = state;
				// }else{
				// 	word->_state = itr->second;
				// }
				increment_tag_trigram_count_by_words(word_vec[pos - 2], word_vec[pos - 1], word_vec[pos]);
				word_set.insert(word->_id);
				// 同じタグの単語集合をカウント
				increment_tag_word_count(word->_state, word->_id);
			}
		}
		_num_words = word_set.size();
	}
	void HMM::set_Wt_for_tag(int tag_id, int number){
		assert(_Wt != NULL);
		assert(1 <= tag_id && tag_id <= _num_tags);
		_Wt[tag_id] = number;
	}
	void HMM::set_num_tags(int n){
		assert(n > 0);
		_num_tags = n;
	}
	void HMM::increment_tag_trigram_count_by_words(Word* tri_word, Word* bi_word, Word* uni_word){
		assert(uni_word != NULL);
		assert(bi_word != NULL);
		assert(tri_word != NULL);
		_unigram_counts[uni_word->_state] += 1;
		_bigram_counts[bi_word->_state][uni_word->_state] += 1;
		_trigram_counts[tri_word->_state][bi_word->_state][uni_word->_state] += 1;
	}
	void HMM::increment_tag_word_count(int tag_id, int word_id){
		assert(1 <= tag_id && tag_id <= _num_tags);
		std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			word_counts[word_id] = 1;
			return;
		}
		itr->second += 1;
	}
	void HMM::decrement_tag_word_count(int tag_id, int word_id){
		assert(1 <= tag_id && tag_id <= _num_tags);
		std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			std::cout << "品詞-単語ペアのカウントが正しく実装されていません." << std::endl;
			exit(1);
		}
		itr->second -= 1;
		if(itr->second <= 0){
			word_counts.erase(itr);
		}
	}
	int HMM::get_count_of_tag_word(int tag_id, int word_id){
		assert(1 <= tag_id && tag_id <= _num_tags);
		std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		auto itr = word_counts.find(word_id);
		if(itr == word_counts.end()){
			return 0;
		}
		return itr->second;
	}
	int HMM::get_word_types_of_tag(int tag_id){
		assert(1 <= tag_id && tag_id <= _num_tags);
		std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
		return word_counts.size();
	}
	double HMM::compute_log_p_t_given_alpha(std::vector<Word*> &word_vec, double alpha){
		double log_Pt_alpha = 0;
		for(int i = 2;i < word_vec.size() - 2;i++){	// <bos>と<eos>の内側だけ考える
			int ti_2 = word_vec[i - 2]->_state;
			int ti_1 = word_vec[i - 1]->_state;
			int ti = word_vec[i]->_state;
			double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
			double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
			double Pt_i_alpha = (n_ti_2_ti_1_ti + alpha) / (n_ti_2_ti_1 + _num_tags * alpha);
			log_Pt_alpha += log(Pt_i_alpha);
		}
		return log_Pt_alpha;
	}
	// 正規化定数で割る前の値
	double HMM::compute_p_wi_given_ti_beta(int ti, int wi, double beta){
		double n_ti_wi = get_count_of_tag_word(ti, wi);
		double n_ti = _unigram_counts[ti];
		double W_ti = _Wt[ti];
		return (n_ti_wi + beta) / (n_ti + W_ti * beta);
	}
	double HMM::compute_p_ti_given_t_alpha(int ti, int ti_1, int ti_2, double alpha){
		double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
		double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
		return (n_ti_2_ti_1_ti + alpha) / (n_ti_2_ti_1 + _num_tags * alpha);
	}
	// in:  t_{i-2},t_{i-1},ti,t_{i+1},t_{i+2},w_i
	void HMM::add_tag_trigram_to_model(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
		assert(1 <= ti && ti <= _num_tags);
		// 1-gram
		_unigram_counts[ti] += 1;
		// 2-gram
		_bigram_counts[ti_1][ti] += 1;
		_bigram_counts[ti][ti1] += 1;
		// 3-gram
		_trigram_counts[ti_2][ti_1][ti] += 1;
		_trigram_counts[ti_1][ti][ti1] += 1;
		_trigram_counts[ti][ti1][ti2] += 1;
		// 品詞-単語ペア
		increment_tag_word_count(ti, wi);
	}
	// in:  t_{i-2},t_{i-1},ti,t_{i+1},t_{i+2},w_i
	void HMM::remove_tag_trigram_from_model(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
		assert(1 <= ti && ti <= _num_tags);
		// 1-gram
		_unigram_counts[ti] -= 1;
		assert(_unigram_counts[ti] >= 0);
		// 2-gram
		_bigram_counts[ti_1][ti] -= 1;
		assert(_bigram_counts[ti_1][ti] >= 0);
		_bigram_counts[ti][ti1] -= 1;
		assert(_bigram_counts[ti][ti1] >= 0);
		// 3-gram
		_trigram_counts[ti_2][ti_1][ti] -= 1;
		assert(_trigram_counts[ti_2][ti_1][ti] >= 0);
		_trigram_counts[ti_1][ti][ti1] -= 1;
		assert(_trigram_counts[ti_1][ti][ti1] >= 0);
		_trigram_counts[ti][ti1][ti2] -= 1;
		assert(_trigram_counts[ti][ti1][ti2] >= 0);
		// 品詞-単語ペア
		decrement_tag_word_count(ti, wi);
	}
	void HMM::perform_gibbs_sampling_with_sequence(std::vector<Word*> &word_vec){
		if(_sampling_table == NULL){
			_sampling_table = new double[_num_tags + 1];
		}
		for(int i = 2;i < word_vec.size() - 2;i++){	// <s>と</s>の内側だけ考える
			int ti_2 = word_vec[i - 2]->_state;
			int ti_1 = word_vec[i - 1]->_state;
			int ti = word_vec[i]->_state;
			int wi = word_vec[i]->_id;
			int ti1 = word_vec[i + 1]->_state;
			int ti2 = word_vec[i + 2]->_state;
			// t_iをモデルパラメータから除去
			remove_tag_trigram_from_model(ti_2, ti_1, ti, ti1, ti2, wi);
			// t_iを再サンプリング
			double sum_prob = 0;
			int new_ti = 0;
			for(int tag = 1;tag <= _num_tags;tag++){
				_sampling_table[tag] = 1;
				double n_ti_wi = get_count_of_tag_word(tag, wi);
				double W_ti = _Wt[tag];
				// n.
				double n_ti_2_ti_1_ti 	= _trigram_counts[ti_2][ti_1][tag];
				double n_ti_2_ti_1 		= _bigram_counts[ti_2][ti_1];
				double n_ti_1_ti_ti1 	= _trigram_counts[ti_1][tag][ti1];
				double n_ti_1_ti 		= _bigram_counts[ti_1][tag];
				double n_ti 			= _unigram_counts[tag];
				double n_ti_ti1_ti2 	= _trigram_counts[tag][ti1][ti2];
				double n_ti_ti1 		= _bigram_counts[tag][ti1];
				// I(.)
				double I_ti_2_ti_1_ti_ti1 			= (ti_2 == ti_1 == tag == ti1) ? 1 : 0;
				double I_ti_2_ti_1_ti 				= (ti_2 == ti_1 == tag) ? 1 : 0;
				double I_ti_2_ti_ti2_and_ti_1_ti1 	= (ti_2 == tag == ti2 && ti_1 == ti1) ? 1 : 0;
				double I_ti_1_ti_ti1_ti2 			= (ti_1 == tag == ti1 == ti2) ? 1 : 0;
				double I_ti_2_ti_and_ti_1_ti1 		= (ti_2 == tag && ti_1 == ti1) ? 1 : 0;
				double I_ti_1_ti_ti1 				= (ti_1 == tag == ti1) ? 1 : 0;
				// 確率を計算
				_sampling_table[tag] *= (n_ti_wi + _beta[tag]) / (n_ti + W_ti * _beta[tag]);
				_sampling_table[tag] *= (n_ti_2_ti_1_ti + _alpha) / (n_ti_2_ti_1 + _num_tags * _alpha);
				_sampling_table[tag] *= (n_ti_1_ti_ti1 + I_ti_2_ti_1_ti_ti1 + _alpha) / (n_ti_1_ti + I_ti_2_ti_1_ti + _num_tags * _alpha);
				_sampling_table[tag] *= (n_ti_ti1_ti2 + I_ti_2_ti_ti2_and_ti_1_ti1 + I_ti_1_ti_ti1_ti2 + _alpha) / (n_ti_ti1 + I_ti_2_ti_and_ti_1_ti1 + I_ti_1_ti_ti1 + _num_tags * _alpha);
				// アニーリング
				_sampling_table[tag] = pow(_sampling_table[tag], 1.0 / _temperature);
				sum_prob += _sampling_table[tag];
			}
			assert(sum_prob > 0);
			double normalizer = 1.0 / sum_prob;
			double bernoulli = sampler::uniform(0, 1);
			double stack = 0;
			for(int tag = 1;tag <= _num_tags;tag++){
				stack += _sampling_table[tag] * normalizer;
				if(stack >= bernoulli){
					new_ti = tag;
					break;
				}
			}
			assert(1 <= new_ti && new_ti <= _num_tags);
			// 新しいt_iをモデルパラメータに追加
			add_tag_trigram_to_model(ti_2, ti_1, new_ti, ti1, ti2, wi);
			word_vec[i]->_state = new_ti;
		}
	}
	Word* HMM::get_random_word_with_tag_from_corpus(int tag, std::vector<std::vector<Word*>> &dataset){
		int random_index = sampler::uniform_int(0, dataset.size() - 1);
		std::vector<Word*> &word_vec = dataset[random_index];
		for(int t = 0;t < word_vec.size();t++){
			int ti = word_vec[t]->_state;
			if(ti == tag){
				return word_vec[t];
			}
		}
		return NULL;
	}
	// 新しいAlphaをサンプリング
	void HMM::sample_new_alpha(std::vector<std::vector<Word*>> &dataset){
		double new_alpha = sampler::normal(_alpha, 0.1 * _alpha);
		int random_index = sampler::uniform_int(0, dataset.size() - 1);
		std::vector<Word*> &word_vec = dataset[random_index];
		// メトロポリス・ヘイスティングス法
		// http://ebsa.ism.ac.jp/ebooks/sites/default/files/ebook/1881/pdf/vol3_ch10.pdf
		// 提案分布は正規分布
		double log_Pt_alpha = compute_log_p_t_given_alpha(word_vec, _alpha);
		double log_Pt_new_alpha = compute_log_p_t_given_alpha(word_vec, new_alpha);
		// q(alpha|new_alpha) / q(new_alpha|alpha)の計算
		double sigma_alpha = 0.1 * _alpha;
		double sigma_new_alpha = 0.1 * new_alpha;
		double var_alpha = sigma_alpha * sigma_alpha;
		double var_new_alpha = sigma_new_alpha * sigma_new_alpha;
		double correcting_term = (_alpha / new_alpha) * exp(
			  0.5 * (new_alpha - _alpha) * (new_alpha - _alpha) / var_alpha
			+ 0.5 * (_alpha - new_alpha) * (_alpha - new_alpha) / var_new_alpha
		);
		if(log_Pt_new_alpha == 0){
			return;
		}
		// 採択率
		double adoption_rate = std::min(1.0, exp(log_Pt_new_alpha - log_Pt_alpha) * correcting_term);
		double bernoulli = sampler::uniform(0, 1);
		if(bernoulli < adoption_rate){
			_alpha = new_alpha;
		}
	}
	// 新しいBetaをサンプリング
	void HMM::sample_new_beta(std::vector<std::vector<Word*>> &dataset){
		for(int tag = 0;tag < _num_tags;tag++){
			double beta = _beta[tag];
			double new_beta = sampler::normal(beta, 0.1 * beta);
			Word* random_word = NULL;
			int limit = 100;
			while(random_word == NULL){
				if (PyErr_CheckSignals() != 0) {		// ctrl+cが押されたかチェック
					return;
				}
				random_word = get_random_word_with_tag_from_corpus(tag, dataset);
				limit--;
				if(limit < 0){
					break;
				}
			}
			if(random_word == NULL){
				continue;
			}
			// メトロポリス・ヘイスティングス法
			// http://ebsa.ism.ac.jp/ebooks/sites/default/files/ebook/1881/pdf/vol3_ch10.pdf
			// 提案分布は正規分布
			double Pti_wi_beta = compute_p_wi_given_ti_beta(random_word->_state, random_word->_id, beta);
			double Pti_wi_new_beta = compute_p_wi_given_ti_beta(random_word->_state, random_word->_id, new_beta);
			// q(beta|new_beta) / q(new_beta|beta)の計算
			double sigma_beta = 0.1 * beta;
			double sigma_new_beta = 0.1 * new_beta;
			double var_beta = sigma_beta * sigma_beta;
			double var_new_beta = sigma_new_beta * sigma_new_beta;
			double correcting_term = (beta / new_beta) * exp(
				  0.5 * (new_beta - beta) * (new_beta - beta) / var_beta
				+ 0.5 * (beta - new_beta) * (beta - new_beta) / var_new_beta
			);
			// 採択率
			double adoption_rate = std::min(1.0, Pti_wi_new_beta * correcting_term / Pti_wi_beta);
			double bernoulli = sampler::uniform(0, 1);
			if(bernoulli < adoption_rate){
				_beta[tag] = new_beta;
			}
		}
	}
	int HMM::get_most_co_occurring_tag(int word_id){
		int max_count = 0;
		int most_co_occurring_tag_id = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			int count = get_count_of_tag_word(tag, word_id);
			if(count > max_count){
				max_count = count;
				most_co_occurring_tag_id = tag;
			}
		}
		return most_co_occurring_tag_id;
	}
	void HMM::dump_trigram_counts(){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
					std::cout << (boost::format("3-gram [%d][%d][%d] = %d") % tri_tag % bi_tag % uni_tag % _trigram_counts[tri_tag][bi_tag][uni_tag]).str() << std::endl;
				}
			}
		}
	}
	void HMM::dump_bigram_counts(){
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
				std::cout << (boost::format("2-gram [%d][%d] = %d") % bi_tag % uni_tag % _bigram_counts[bi_tag][uni_tag]).str() << std::endl;
			}
		}
	}
	void HMM::dump_unigram_counts(){
		for(int uni_tag = 0;uni_tag < _num_tags;uni_tag++){
			std::cout << (boost::format("1-gram [%d] = %d") % uni_tag % _unigram_counts[uni_tag]).str() << std::endl;
		}
	}
	void HMM::dump_word_types(){
		for(int tag = 0;tag < _num_tags;tag++){
			std::cout << tag << ": " << get_word_types_of_tag(tag) << std::endl;
		}
	}
	template <class Archive>
	void HMM::serialize(Archive &ar, unsigned int version){
		boost::serialization::split_free(ar, *this, version);
	}
	bool HMM::save(std::string filename){
		bool success = false;
		std::ofstream ofs(filename);
		if(ofs.good()){
			boost::archive::binary_oarchive oarchive(ofs);
			oarchive << *this;
			success = true;
		}
		ofs.close();
		return success;
	}
	bool HMM::load(std::string filename){
		bool success = false;
		std::ifstream ifs(filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> *this;
			success = true;
		}
		ifs.close();
		return success;
	}
}

namespace boost { 
	namespace serialization {
		template<class Archive>
		void save(Archive &ar, const bhmm::HMM &hmm, unsigned int version) {
			ar & hmm._num_tags;
			ar & hmm._num_words;
			ar & hmm._alpha;
			ar & hmm._temperature;
			ar & hmm._minimum_temperature;
			ar & hmm._tag_word_counts;
			ar & hmm._allocated;
			if(hmm._allocated){
				assert(hmm._num_tags > 0);
				int num_tags = hmm._num_tags;
				// 各タグの可能な単語数
				for(int tag = 0;tag < num_tags;tag++){
					ar & hmm._Wt[tag];
				}
				// Betaの初期化
				// 初期値は1
				for(int tag = 0;tag < num_tags;tag++){
					ar & hmm._beta[tag];
				}
				// 3-gram		
				for(int tri_tag = 0;tri_tag < num_tags;tri_tag++){
					for(int bi_tag = 0;bi_tag < num_tags;bi_tag++){
						for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
							ar & hmm._trigram_counts[tri_tag][bi_tag][uni_tag];
						}
					}
				}
				// 2-gram
				for(int bi_tag = 0;bi_tag < num_tags;bi_tag++){
					for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
						ar & hmm._bigram_counts[bi_tag][uni_tag];
					}
				}
				// 1-gram
				for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
					ar & hmm._unigram_counts[uni_tag];
				}
			}
		}
		template<class Archive>
		void load(Archive &ar, bhmm::HMM &hmm, unsigned int version) {
			ar & hmm._num_tags;
			ar & hmm._num_words;
			ar & hmm._alpha;
			ar & hmm._temperature;
			ar & hmm._minimum_temperature;
			ar & hmm._tag_word_counts;
			ar & hmm._allocated;
			if(hmm._allocated){
				assert(hmm._num_tags > 0);
				int num_tags = hmm._num_tags;
				// 各タグの可能な単語数
				hmm._Wt = new int[num_tags];
				for(int tag = 0;tag < num_tags;tag++){
					ar & hmm._Wt[tag];
				}
				// Betaの初期化
				// 初期値は1
				hmm._beta = new double[num_tags];
				for(int tag = 0;tag < num_tags;tag++){
					ar & hmm._beta[tag];
				}
				// 3-gram		
				hmm._trigram_counts = new int**[num_tags];
				for(int tri_tag = 0;tri_tag < num_tags;tri_tag++){
					hmm._trigram_counts[tri_tag] = new int*[num_tags];
					for(int bi_tag = 0;bi_tag < num_tags;bi_tag++){
						hmm._trigram_counts[tri_tag][bi_tag] = new int[num_tags];
						for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
							ar & hmm._trigram_counts[tri_tag][bi_tag][uni_tag];
						}
					}
				}
				// 2-gram
				hmm._bigram_counts = new int*[num_tags];
				for(int bi_tag = 0;bi_tag < num_tags;bi_tag++){
					hmm._bigram_counts[bi_tag] = new int[num_tags];
					for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
						ar & hmm._bigram_counts[bi_tag][uni_tag];
					}
				}
				// 1-gram
				hmm._unigram_counts = new int[num_tags];
				for(int uni_tag = 0;uni_tag < num_tags;uni_tag++){
					ar & hmm._unigram_counts[uni_tag];
				}
			}
		}
	}
}