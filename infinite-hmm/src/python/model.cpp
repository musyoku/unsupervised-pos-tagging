#include <iostream>
#include "../ihmm/utils.h"
#include "model.h"

namespace ihmm {
	Model::Model(int num_initial_tags, Dataset* dataset){
		_set_locale();
		_hmm = new InfiniteHMM(num_initial_tags, dataset->get_num_words());
		_hmm->initialize_with_training_dataset(dataset->_word_sequences_train);
	}
	Model::Model(std::string filename){
		_set_locale();
		_hmm = new InfiniteHMM();
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
		_hmm->_alpha = alpha;
	}
	void Model::set_initial_beta(double beta){
		_hmm->_beta = beta;
	}
	void Model::set_initial_gamma(double gamma){
		_hmm->_gamma = gamma;
	}
	void Model::set_initial_gamma_emission(double gamma_emission){
		_hmm->_gamma_emission = gamma_emission;
	}
	void Model::set_initial_beta_emission(double beta_emission){
		_hmm->_beta_emission = beta_emission;
	}
	int Model::get_num_tags(){
		return _hmm->get_num_tags();
	}
	// 文の確率
	// 前向きアルゴリズム
	double Model::compute_p_sentence(std::vector<Word*> &sentence, double** forward_table){
		assert(sentence.size() > 2);	// <s>と</s>
		int tag_bos = 0;	// <s>
		for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
			if(_hmm->is_tag_new(ti)){
				forward_table[1][ti] = 0;
				continue;
			}
			int ti_1 = tag_bos;	// <s>
			int wi = sentence[1]->_id;
			double p_transition = _hmm->compute_p_tag_given_context(ti, ti_1);
			double p_emission = _hmm->compute_p_word_given_tag(wi, ti);
			assert(p_transition > 0);
			assert(p_emission > 0);
			forward_table[1][ti] = p_emission * p_transition;
		}
		for(int i = 2;i < sentence.size() - 1;i++){
			for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
				if(_hmm->is_tag_new(ti)){
					forward_table[i][ti] = 0;
					continue;
				}
				int wi = sentence[i]->_id;
				double p_emission = _hmm->compute_p_word_given_tag(wi, ti);
				assert(p_emission > 0);
				forward_table[i][ti] = 0;
				for(int ti_1 = 1;ti_1 <= _hmm->get_num_tags();ti_1++){
					forward_table[i][ti] += forward_table[i - 1][ti_1] * _hmm->compute_p_tag_given_context(ti, ti_1);
				}
				forward_table[i][ti] *= p_emission;
			}
		}
		int i = sentence.size() - 2;
		double p_x = 0;
		for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
			p_x += forward_table[i][ti];
		}
		return p_x;
	}
	void Model::_alloc_viterbi_tables(int sentence_length, double** &forward_table, double** &decode_table){
		forward_table = new double*[sentence_length];
		decode_table = new double*[sentence_length];
		for(int i = 0;i < sentence_length;i++){
			forward_table[i] = new double[_hmm->get_num_tags() + 1];
			decode_table[i] = new double[_hmm->get_num_tags() + 1];
		}
	}
	void Model::_free_viterbi_tables(int sentence_length, double** &forward_table, double** &decode_table){
		for(int i = 0;i < sentence_length;i++){
			delete[] forward_table[i];
			delete[] decode_table[i];
		}
		delete[] forward_table;
		delete[] decode_table;
	}
	boost::python::list Model::python_viterbi_decode(boost::python::list py_word_ids){
		// デコード用のテーブルを確保
		int num_words = boost::python::len(py_word_ids);
		double** forward_table = NULL;
		double** decode_table = NULL;
		_alloc_viterbi_tables(num_words + 2, forward_table, decode_table);
		// Python側から渡された単語IDリストを変換
		std::vector<Word*> sentence;
		// <s>をセット
		Word* bos = new Word();
		bos->_tag = 0;
		sentence.push_back(bos);
		for(int i = 0;i < num_words;i++){
			Word* word = new Word();
			word->_id = boost::python::extract<int>(py_word_ids[i]);
			word->_tag = 0;
			sentence.push_back(word);
		}
		// </s>をセット
		Word* eos = new Word();
		eos->_tag = 0;
		sentence.push_back(eos);
		// ビタビアルゴリズム
		std::vector<int> sampled_state_sequence;
		viterbi_decode(sentence, sampled_state_sequence, forward_table, decode_table);
		// 結果を返す
		boost::python::list result;
		for(int i = 0;i < sampled_state_sequence.size();i++){
			result.append(sampled_state_sequence[i]);
		}
		_free_viterbi_tables(num_words + 2, forward_table, decode_table);
		for(int i = 0;i < sentence.size();i++){
			delete sentence[i];
		}
		return result;
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence){
		double** forward_table = NULL;
		double** decode_table = NULL;
		_alloc_viterbi_tables(sentence.size(), forward_table, decode_table);
		viterbi_decode(sentence, sampled_state_sequence, forward_table, decode_table);
		_free_viterbi_tables(sentence.size(), forward_table, decode_table);
	}
	void Model::viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double** forward_table, double** decode_table){
		assert(sentence.size() > 4);	// <s>と</s>それぞれ2つづつ
		int tag_bos = 0;	// <s>
		for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
			int ti_1 = tag_bos;	// <s>
			int wi = sentence[1]->_id;
			double p_transition = _hmm->compute_p_tag_given_context(ti, ti_1);
			double p_emission = _hmm->compute_p_word_given_tag(wi, ti);
			assert(p_transition > 0);
			double log_p_emission = -1000000;
			if(p_emission > 0){
				log_p_emission = log(p_emission);
			}
			forward_table[1][ti] = log_p_emission + log(p_transition);
		}
		for(int i = 2;i < sentence.size() - 1;i++){
			for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
				int wi = sentence[i]->_id;
				double p_emission = _hmm->compute_p_word_given_tag(wi, ti);
				double log_p_emission = -1000000;
				if(p_emission > 0){
					log_p_emission = log(p_emission);
				}
				double max_value = 0;
				for(int ti_1 = 1;ti_1 <= _hmm->get_num_tags();ti_1++){
					double p_transition = _hmm->compute_p_tag_given_context(ti, ti_1);
					double value = log(p_transition) + forward_table[i - 1][ti_1];
					if(max_value == 0 || value > max_value){
						max_value = value;
						forward_table[i][ti] = value + log_p_emission;
						decode_table[i][ti] = ti_1;
					}
				}
			}
		}
		int i = sentence.size() - 2;
		double max_p_x_s = 0;
		double argmax_ti = 0;
		for(int ti = 1;ti <= _hmm->get_num_tags();ti++){
			double log_p_x_s = forward_table[i][ti];
			if(max_p_x_s == 0 || log_p_x_s > max_p_x_s){
				max_p_x_s = log_p_x_s;
				argmax_ti = ti;
			}
		}
		assert(1 <= argmax_ti && argmax_ti <= _hmm->get_num_tags());
		sampled_state_sequence.clear();
		sampled_state_sequence.push_back(argmax_ti);
		int ti = argmax_ti;
		for(int i = sentence.size() - 2;i >= 2;i--){
			int ti_1 = decode_table[i][ti];
			sampled_state_sequence.push_back(ti_1);
			ti = ti_1;
		}
		std::reverse(sampled_state_sequence.begin(), sampled_state_sequence.end());
		assert(sampled_state_sequence.size() == sentence.size() - 2);
	}
	struct value_comparator {
		bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) {
			return a.second > b.second;
		}   
	};
	void Model::print_typical_words_assigned_to_each_tag(int number_to_show, Dictionary* dict){
		using std::wcout;
		using std::endl;
		for(int tag = 1;tag <= _hmm->get_num_tags();tag++){
			if(_hmm->is_tag_new(tag)){
				continue;
			}
			int n = 0;
			wcout << "\x1b[32;1m" << "[" << tag << "]" << "\x1b[0m" << std::endl;
			std::multiset<std::pair<int, int>, value_comparator> ranking;
			for(int word_id = 0;word_id < _hmm->get_num_words();word_id++){
				if(dict->is_unk(word_id)){
					continue;
				}
				int count = _hmm->_m_iq_tables[tag][word_id]->get_num_customers();
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
	// void Model::print_alpha_and_beta(){
	// 	using std::cout;
	// 	using std::endl;
	// 	cout << "\x1b[1m" << "alpha" << "\x1b[0m " << _hmm->_alpha << std::endl;
	// 	for(int tag = 1;tag <= _hmm->get_num_tags();tag++){
	// 		cout << "\x1b[1m" << "beta[" << tag << "]" << "\x1b[0m " << _hmm->_beta[tag] << std::endl;
	// 	}
	// }
}