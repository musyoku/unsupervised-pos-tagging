#include <iostream>
#include "../bhmm/utils.h"
#include "model.h"

namespace bhmm {
	Model::Model(int num_tags, Dataset* dataset, boost::python::list py_Wt){
		_set_locale();
		_hmm = new HMM(num_tags, dataset->get_num_words());
		std::vector<int> Wt = utils::vector_from_list<int>(py_Wt);
		_hmm->initialize_with_training_corpus(dataset->_word_sequences_train, Wt);
	}
	Model::Model(int num_tags, Dataset* dataset, std::vector<int> &Wt){
		_set_locale();
		_hmm = new HMM(num_tags, dataset->get_num_words());
		_hmm->initialize_with_training_corpus(dataset->_word_sequences_train, Wt);
	}
	Model::Model(std::string filename){
		_set_locale();
		_hmm = new HMM();
		assert(load(filename) == true);
	}
	Model::~Model(){
		delete _hmm;
	}
	// 日本語周り
	void Model::_set_locale(){
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);
	}
	bool Model::load(std::string filename){
		return _hmm->load(filename);
	}
	bool Model::save(std::string filename){
		return _hmm->save(filename);
	}
	void Model::set_initial_alpha(double alpha){
		_hmm->set_alpha(alpha);
	}
	void Model::set_initial_beta(double beta){
		_hmm->set_beta(beta);
	}
	int Model::get_num_tags(){
		return _hmm->_num_tags;
	}
	double Model::get_temperature(){
		return _hmm->_temperature;
	}
	void Model::set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void Model::set_minimum_temperature(double temperature){
		_hmm->_minimum_temperature = temperature;
	}
	void Model::anneal_temperature(double temperature){
		_hmm->anneal_temperature(temperature);
	}
	// 文の確率
	// 前向きアルゴリズムの拡張
	double Model::compute_p_sentence(std::vector<Word*> &sentence, double*** forward_table){
		assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		int tag_bos = 0;	// <s>
		for(int ti = 1;ti <= _hmm->_num_tags;ti++){
			int ti_2 = tag_bos;	// <s>
			int ti_1 = tag_bos;	// <s>
			id wi = sentence[2]->_id;
			double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
			double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
			assert(p_s_given_prev > 0);
			assert(p_w_given_s > 0);
			forward_table[2][tag_bos][ti] = p_w_given_s * p_s_given_prev;
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				forward_table[2][ti_1][ti] = 0;
			}
		}
		for(int i = 3;i < sentence.size() - 2;i++){
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				for(int ti = 1;ti <= _hmm->_num_tags;ti++){

					id wi = sentence[i]->_id;
					double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
					assert(p_w_given_s > 0);
					forward_table[i][ti_1][ti] = 0;
					if(i == 3){
						forward_table[i][ti_1][ti] += forward_table[i - 1][tag_bos][ti_1] * _hmm->compute_p_ti_given_t(ti, ti_1, tag_bos);
					}else{
						for(int ti_2 = 1;ti_2 <= _hmm->_num_tags;ti_2++){
							forward_table[i][ti_1][ti] += forward_table[i - 1][ti_2][ti_1] * _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
						}
					}
					forward_table[i][ti_1][ti] *= p_w_given_s;

				}
			}
		}
		int i = sentence.size() - 3;
		double p_x = 0;
		for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
			for(int ti = 1;ti <= _hmm->_num_tags;ti++){
				p_x += forward_table[i][ti_1][ti];
			}
		}
		return p_x;
	}
	void Model::_alloc_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table){
		forward_table = new double**[sentence_length];
		decode_table = new double**[sentence_length];
		for(int i = 0;i < sentence_length;i++){
			forward_table[i] = new double*[_hmm->_num_tags + 1];
			decode_table[i] = new double*[_hmm->_num_tags + 1];
			for(int k = 0;k <= _hmm->_num_tags;k++){
				forward_table[i][k] = new double[_hmm->_num_tags + 1];
				decode_table[i][k] = new double[_hmm->_num_tags + 1];
			}
		}
	}
	void Model::_free_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table){
		for(int i = 0;i < sentence_length;i++){
			for(int k = 0;k <= _hmm->_num_tags;k++){
				delete[] forward_table[i][k];
				delete[] decode_table[i][k];
			}
			delete[] forward_table[i];
			delete[] decode_table[i];
		}
		delete[] forward_table;
		delete[] decode_table;
	}
	boost::python::list Model::python_viterbi_decode(boost::python::list py_word_ids){
		// デコード用のテーブルを確保
		int num_words = boost::python::len(py_word_ids);
		double*** forward_table = NULL;
		double*** decode_table = NULL;
		_alloc_viterbi_tables(num_words + 4, forward_table, decode_table);
		// Python側から渡された単語IDリストを変換
		std::vector<Word*> sentence;
		// <s>を2つセット
		for(int i = 0;i < 2;i++){
			Word* bos = new Word();
			bos->_state = 0;
			sentence.push_back(bos);
		}
		for(int i = 0;i < num_words;i++){
			Word* word = new Word();
			word->_id = boost::python::extract<id>(py_word_ids[i]);
			word->_state = 0;
			sentence.push_back(word);
		}
		// </s>を2つセット
		for(int i = 0;i < 2;i++){
			Word* eos = new Word();
			eos->_state = 0;
			sentence.push_back(eos);
		}
		// ビタビアルゴリズム
		std::vector<int> sampled_state_sequence;
		viterbi_decode(sentence, sampled_state_sequence, forward_table, decode_table);
		// 結果を返す
		boost::python::list result;
		for(int i = 0;i < sampled_state_sequence.size();i++){
			result.append(sampled_state_sequence[i]);
		}
		_free_viterbi_tables(num_words + 4, forward_table, decode_table);
		for(int i = 0;i < sentence.size();i++){
			delete sentence[i];
		}
		return result;
	}
	// 状態系列の復号
	// ビタビアルゴリズムの拡張
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence){
		double*** forward_table = NULL;
		double*** decode_table = NULL;
		_alloc_viterbi_tables(sentence.size(), forward_table, decode_table);
		viterbi_decode(sentence, sampled_state_sequence, forward_table, decode_table);
		_free_viterbi_tables(sentence.size(), forward_table, decode_table);
	}
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double*** forward_table, double*** decode_table){
		assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		int tag_bos = 0;	// <s>
		for(int ti = 1;ti <= _hmm->_num_tags;ti++){
			int ti_2 = tag_bos;	// <s>
			int ti_1 = tag_bos;	// <s>
			id wi = sentence[2]->_id;
			double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
			double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
			assert(p_s_given_prev > 0);
			assert(p_w_given_s > 0);
			forward_table[2][tag_bos][ti] = log(p_w_given_s) + log(p_s_given_prev);
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				forward_table[2][ti_1][ti] = -10000000;
			}
		}
		for(int i = 3;i < sentence.size() - 2;i++){
			for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
				for(int ti = 1;ti <= _hmm->_num_tags;ti++){

					id wi = sentence[i]->_id;
					double p_w_given_s = _hmm->compute_p_wi_given_ti(wi, ti);
					assert(p_w_given_s > 0);
					if(i == 3){
						double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, tag_bos);
						forward_table[i][ti_1][ti] = forward_table[i - 1][tag_bos][ti_1] + log(p_s_given_prev) + log(p_w_given_s);
						decode_table[i][ti_1][ti] = tag_bos;
					}else{
						double max_value = 0;
						for(int ti_2 = 1;ti_2 <= _hmm->_num_tags;ti_2++){
							double p_s_given_prev = _hmm->compute_p_ti_given_t(ti, ti_1, ti_2);
							double value = log(p_s_given_prev) + forward_table[i - 1][ti_2][ti_1];
							if(max_value == 0 || value > max_value){
								max_value = value;
								forward_table[i][ti_1][ti] = value + log(p_w_given_s);
								decode_table[i][ti_1][ti] = ti_2;
							}
						}
					}
				}
			}
		}
		int i = sentence.size() - 3;
		double max_p_x_s = 0;
		double argmax_ti_1 = 0;
		double argmax_ti = 0;
		for(int ti_1 = 1;ti_1 <= _hmm->_num_tags;ti_1++){
			for(int ti = 1;ti <= _hmm->_num_tags;ti++){
				double log_p_x_s = forward_table[i][ti_1][ti];
				if(max_p_x_s == 0 || log_p_x_s > max_p_x_s){
					max_p_x_s = log_p_x_s;
					argmax_ti_1 = ti_1;
					argmax_ti = ti;
				}
			}
		}
		assert(1 <= argmax_ti_1 && argmax_ti_1 <= _hmm->_num_tags);
		assert(1 <= argmax_ti && argmax_ti <= _hmm->_num_tags);
		sampled_state_sequence.clear();
		sampled_state_sequence.push_back(argmax_ti_1);
		sampled_state_sequence.push_back(argmax_ti);
		int ti_1 = argmax_ti_1;
		int ti = argmax_ti;
		for(int i = sentence.size() - 3;i >= 4;i--){
			int ti_2 = decode_table[i][ti_1][ti];
			sampled_state_sequence.push_back(ti_2);
			ti = ti_1;
			ti_1 = ti_2;
		}
		std::reverse(sampled_state_sequence.begin(), sampled_state_sequence.end());
		assert(sampled_state_sequence.size() == sentence.size() - 4);
	}
	struct value_comparator {
		bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
			return a.second > b.second;
		}   
	};
	void Model::show_typical_words_of_each_tag(int number_to_show, Dictionary* dict){
		using std::wcout;
		using std::endl;
		for(int tag = 1;tag <= _hmm->_num_tags;tag++){
			int n = 0;
			wcout << "\x1b[32;1m" << "[" << tag << "]" << "\x1b[0m" << std::endl;
			std::multiset<std::pair<int, int>, value_comparator> ranking;
			for(id word_id = 0;word_id < _hmm->_num_words;word_id++){
				int count = _hmm->_tag_word_counts[tag][word_id];
				if(count > 0){
					ranking.insert(std::make_pair(word_id, count));
				}
			}
			for(auto elem: ranking){
				std::wstring word = dict->word_id_to_string(elem.first);
				wcout << "\x1b[1m" << word << "\x1b[0m" << L"(" << elem.second << L") ";
				n++;
				if(n > number_to_show){
					break;
				}
			}
			wcout << endl;
		}
	}
}