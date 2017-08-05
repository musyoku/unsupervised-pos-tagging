#include <iostream>
#include <fstream>
#include "hmm.h"

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
}
HMM::~HMM(){
	if(_trigram_counts != NULL){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				free(_trigram_counts[tri_tag][bi_tag]);
			}
			free(_trigram_counts[tri_tag]);
		}
		free(_trigram_counts);
	}
	if(_bigram_counts != NULL){
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			free(_bigram_counts);
		}
		free(_bigram_counts);
	}
	if(_unigram_counts != NULL){
		free(_unigram_counts);
	}
	if(_sampling_table != NULL){
		free(_sampling_table);
	}
	if(_beta != NULL){
		free(_beta);
	}
	if(_Wt != NULL){
		free(_Wt);
	}
}

template <class Archive>
void HMM::serialize(Archive& ar, unsigned int version)
{
	static_cast<void>(version);
	ar & _num_tags;
	ar & _num_words;
	ar & _alpha;
	ar & _temperature;
	ar & _minimum_temperature;
	ar & _tag_word_counts;
}
void HMM::anneal_temperature(double multiplier){
	if(_temperature > _minimum_temperature){
		_temperature *= multiplier;
	}
}
void HMM::initialize(std::vector<std::vector<Word*>> &dataset){
	// その他
	alloc_table();
	// nグラムのカウントテーブル
	init_ngram_counts(dataset);
}
void HMM::alloc_table(){
	assert(_num_tags != -1);
	// 各タグの可能な単語数
	_Wt = (int*)calloc(_num_tags, sizeof(int));
	// Betaの初期化
	// 初期値は1
	_beta = (double*)malloc(_num_tags * sizeof(double));
	for(int tag = 0;tag < _num_tags;tag++){
		_beta[tag] = 1;
	}
	// 3-gram		
	_trigram_counts = (int***)malloc(_num_tags * sizeof(int**));
	for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
		_trigram_counts[tri_tag] = (int**)malloc(_num_tags * sizeof(int*));
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			_trigram_counts[tri_tag][bi_tag] = (int*)calloc(_num_tags, sizeof(int));
		}
	}
	// 2-gram
	_bigram_counts = (int**)malloc(_num_tags * sizeof(int*));
	for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
		_bigram_counts[bi_tag] = (int*)calloc(_num_tags, sizeof(int));
	}
	// 1-gram
	_unigram_counts = (int*)calloc(_num_tags, sizeof(int));
}
void HMM::init_ngram_counts(std::vector<std::vector<Word*>> &dataset){
	assert(_num_tags != -1);
	// 最初は品詞をランダムに割り当てる
	std::set<int> word_set;
	std::unordered_map<int, int> tag_for_word;
	for(int data_index = 0;data_index < dataset.size();data_index++){
		std::vector<Word*> &line = dataset[data_index];
		// pos < 2
		// <bos>2つ
		_bigram_counts[line[0]->tag_id][line[1]->tag_id] += 1;
		_unigram_counts[line[0]->tag_id] += 1;
		_unigram_counts[line[1]->tag_id] += 1;
		word_set.insert(line[0]->word_id);
		word_set.insert(line[1]->word_id);
		increment_tag_word_count(line[0]->tag_id, line[0]->word_id);
		increment_tag_word_count(line[1]->tag_id, line[1]->word_id);
		// line.size() - 2 > pos >= 2
		for(int pos = 2;pos < line.size() - 2;pos++){	// 3-gramなので3番目から.
			Word* word = line[pos];
			auto itr = tag_for_word.find(word->word_id);
			if(itr == tag_for_word.end()){
				word->tag_id = sampler::uniform_int(0, _num_tags - 1);
				tag_for_word[word->word_id] = word->tag_id;
			}else{
				word->tag_id = itr->second;
			}
			update_ngram_count(line[pos - 2], line[pos - 1], line[pos]);
			word_set.insert(word->word_id);
			// 同じタグの単語集合をカウント
			increment_tag_word_count(word->tag_id, word->word_id);
		}
		// pos >= line.size() - 2
		// <eos>2つ
		int end_index = line.size() - 1;
		int t_end = line[end_index]->tag_id;
		int t_end_1 = line[end_index - 1]->tag_id;
		int t_end_2 = line[end_index - 2]->tag_id;
		int t_end_3 = line[end_index - 3]->tag_id;
		_unigram_counts[t_end] += 1;
		_unigram_counts[t_end_1] += 1;
		_bigram_counts[t_end_2][t_end_1] += 1;
		_bigram_counts[t_end_1][t_end] += 1;
		_trigram_counts[t_end_3][t_end_2][t_end_1] += 1;
		_trigram_counts[t_end_2][t_end_1][t_end] += 1;
		int w_end = line[end_index]->word_id;
		int w_end_1 = line[end_index - 1]->word_id;
		word_set.insert(w_end);
		word_set.insert(w_end_1);
		increment_tag_word_count(t_end, w_end);
		increment_tag_word_count(t_end_1, w_end_1);
	}
	_num_words = word_set.size();
}
void HMM::set_Wt_for_tag(int tag_id, int number){
	assert(_Wt != NULL);
	assert(tag_id < _num_tags);
	_Wt[tag_id] = number;
}
void HMM::update_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word){
	_trigram_counts[tri_word->tag_id][bi_word->tag_id][uni_word->tag_id] += 1;
	_bigram_counts[bi_word->tag_id][uni_word->tag_id] += 1;
	_unigram_counts[uni_word->tag_id] += 1;
}
void HMM::increment_tag_word_count(int tag_id, int word_id){
	std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
	auto itr = word_counts.find(word_id);
	if(itr == word_counts.end()){
		word_counts[word_id] = 1;
		return;
	}
	itr->second += 1;
}
void HMM::decrement_tag_word_count(int tag_id, int word_id){
	std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
	auto itr = word_counts.find(word_id);
	if(itr == word_counts.end()){
		c_printf("[R]%s [*]%s", "エラー", "品詞-単語ペアのカウントが正しく実装されていません.");
		exit(1);
	}
	itr->second -= 1;
	if(itr->second <= 0){
		word_counts.erase(itr);
	}
}
int HMM::get_count_for_tag_word(int tag_id, int word_id){
	std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
	auto itr = word_counts.find(word_id);
	if(itr == word_counts.end()){
		return 0;
	}
	return itr->second;
}
int HMM::get_word_types_for_tag(int tag_id){
	std::unordered_map<int, int> &word_counts = _tag_word_counts[tag_id];
	return word_counts.size();
}
double HMM::compute_log_Pt_alpha(std::vector<Word*> &line, double alpha){
	double log_Pt_alpha = 0;
	for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
		int ti_2 = line[pos - 2]->tag_id;
		int ti_1 = line[pos - 1]->tag_id;
		int ti = line[pos]->tag_id;
		double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
		double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
		double Pt_i_alpha = (n_ti_2_ti_1_ti + alpha) / (n_ti_2_ti_1 + _num_tags * alpha);
		log_Pt_alpha += log(Pt_i_alpha);
	}
	return log_Pt_alpha;
}
double HMM::compute_log_Pw_t_alpha(std::vector<Word*> &line, double alpha){
	double log_Pt_alpha = 0;
	for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
		int ti_2 = line[pos - 2]->tag_id;
		int ti_1 = line[pos - 1]->tag_id;
		int ti = line[pos]->tag_id;
		double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][ti];
		double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
		double Pt_i_alpha = (n_ti_2_ti_1_ti + alpha) / (n_ti_2_ti_1 + _num_tags * alpha);
		log_Pt_alpha += log(Pt_i_alpha);
	}
	return log_Pt_alpha;
}
// 正規化定数で割る前の値
double HMM::compute_Pti_wi_beta(int ti, int wi, double beta){
	double n_ti_wi = get_count_for_tag_word(ti, wi);
	double n_ti = _unigram_counts[ti];
	double W_ti = _Wt[ti];
	return (n_ti_wi + beta) / (n_ti + W_ti * beta);
}
// in:  t_{i-2},t_{i-1},ti,t_{i+1},t_{i+2},w_i
void HMM::add_tag_to_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
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
void HMM::remove_tag_from_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi){
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
void HMM::perform_gibbs_sampling_with_line(std::vector<Word*> &line){
	if(_sampling_table == NULL){
		_sampling_table = (double*)malloc(_num_tags * sizeof(double));
	}
	for(int pos = 2;pos < line.size() - 2;pos++){	// <bos>と<eos>の内側だけ考える
		int ti_2 = line[pos - 2]->tag_id;
		int ti_1 = line[pos - 1]->tag_id;
		int ti = line[pos]->tag_id;
		int wi = line[pos]->word_id;
		int ti1 = line[pos + 1]->tag_id;
		int ti2 = line[pos + 2]->tag_id;
		// t_iをモデルパラメータから除去
		remove_tag_from_model_parameters(ti_2, ti_1, ti, ti1, ti2, wi);
		// t_iを再サンプリング
		double sum = 0;
		int new_ti = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			_sampling_table[tag] = 1;
			double n_ti_wi = get_count_for_tag_word(tag, wi);
			double n_ti = _unigram_counts[tag];
			double W_ti = _Wt[tag];
			double n_ti_2_ti_1_ti = _trigram_counts[ti_2][ti_1][tag];
			double n_ti_2_ti_1 = _bigram_counts[ti_2][ti_1];
			double n_ti_1_ti_ti1 = _trigram_counts[ti_1][tag][ti1];
			double n_ti_1_ti = _bigram_counts[ti_1][tag];
			double I_ti_2_ti_1_ti_ti1 = (ti_2 == ti_1 == tag == ti1) ? 1 : 0;
			double I_ti_2_ti_1_ti = (ti_2 == ti_1 == tag) ? 1 : 0;
			double n_ti_ti1_ti2 = _trigram_counts[tag][ti1][ti2];
			double n_ti_ti1 = _bigram_counts[tag][ti1];
			double I_ti_2_ti_ti2_and_ti_1_ti1 = (ti_2 == tag == ti2 && ti_1 == ti1) ? 1 : 0;
			double I_ti_1_ti_ti1_ti2 = (ti_1 == tag == ti1 == ti2) ? 1 : 0;
			double I_ti_2_ti_and_ti_1_ti1 = (ti_2 == tag && ti_1 == ti1) ? 1 : 0;
			double I_ti_1_ti_ti1 = (ti_1 == tag == ti1) ? 1 : 0;
			_sampling_table[tag] *= (n_ti_wi + _beta[tag]) / (n_ti + W_ti * _beta[tag]);
			_sampling_table[tag] *= (n_ti_2_ti_1_ti + _alpha) / (n_ti_2_ti_1 + _num_tags * _alpha);
			_sampling_table[tag] *= (n_ti_1_ti_ti1 + I_ti_2_ti_1_ti_ti1 + _alpha) / (n_ti_1_ti + I_ti_2_ti_1_ti + _num_tags * _alpha);
			_sampling_table[tag] *= (n_ti_ti1_ti2 + I_ti_2_ti_ti2_and_ti_1_ti1 + I_ti_1_ti_ti1_ti2 + _alpha) / (n_ti_ti1 + I_ti_2_ti_and_ti_1_ti1 + I_ti_1_ti_ti1 + _num_tags * _alpha);
			_sampling_table[tag] = pow(_sampling_table[tag], 1.0 / _temperature);
			sum += _sampling_table[tag];
		}
		assert(sum > 0);
		double normalizer = 1.0 / sum;
		double bernoulli = sampler::uniform(0, 1);
		sum = 0;
		for(int tag = 0;tag < _num_tags;tag++){
			sum += _sampling_table[tag] * normalizer;
			if(sum >= bernoulli){
				new_ti = tag;
				break;
			}
		}
		// 新しいt_iをモデルパラメータに追加
		add_tag_to_model_parameters(ti_2, ti_1, new_ti, ti1, ti2, wi);
		line[pos]->tag_id = new_ti;
	}
}
// 論文(6)式と(7)式を掛けたものからtiをサンプリング
int HMM::sample_tag_from_Pt_w(int ti_2, int ti_1, int wi){
	double sum_p = 0;
	for(int tag = 0;tag < _num_tags;tag++){
		double Pt_alpha = (_trigram_counts[ti_2][ti_1][tag] + _alpha) / (_bigram_counts[ti_2][ti_1] + _num_tags * _alpha);
		double Pw_t_beta = (get_count_for_tag_word(tag, wi) + _beta[tag]) / (_unigram_counts[tag] + _Wt[tag] * _beta[tag]);
		double Ptw_alpha_beta = Pw_t_beta * Pt_alpha;
		_sampling_table[tag] = Ptw_alpha_beta;
		sum_p += Ptw_alpha_beta;
	}
	double normalizer = 1.0 / sum_p;
	double bernoulli = sampler::uniform(0, 1);
	sum_p = 0;
	for(int tag = 0;tag < _num_tags;tag++){
		sum_p += _sampling_table[tag] * normalizer;
		if(sum_p > bernoulli){
			return tag;
		}
	}
	return _num_tags - 1;
}
// 論文(6)式と(7)式を掛けたものからtiをサンプリング
int HMM::argmax_tag_from_Pt_w(int ti_2, int ti_1, int wi){
	double max_p = 0;
	double max_tag = 0;
	// std::cout << (boost::format("argmax(%d, %d, %d)") % ti_2 % ti_1 % wi).str() << std::endl;
	for(int tag = 0;tag < _num_tags;tag++){
		double Pt_alpha = (_trigram_counts[ti_2][ti_1][tag] + _alpha) / (_bigram_counts[ti_2][ti_1] + _num_tags * _alpha);
		double Pw_t_beta = (get_count_for_tag_word(tag, wi) + _beta[tag]) / (_unigram_counts[tag] + _Wt[tag] * _beta[tag]);
		double Ptw_alpha_beta = Pw_t_beta * Pt_alpha;
		// std::cout << (boost::format("%f = %f * %f") % Ptw_alpha_beta % Pw_t_beta % Pt_alpha).str() << std::endl;
		if(Ptw_alpha_beta > max_p){
			max_p = Ptw_alpha_beta;
			max_tag = tag;
		}
	}
	// std::cout << "return " << max_tag << std::endl;
	return max_tag;
}
Word* HMM::_get_random_word_with_tag(int tag, std::vector<std::vector<Word*>> &dataset){
	int random_index = sampler::uniform_int(0, dataset.size() - 1);
	std::vector<Word*> &line = dataset[random_index];
	for(int pos = 0;pos < line.size();pos++){
		int ti = line[pos]->tag_id;
		if(ti == tag){
			return line[pos];
		}
	}
	return NULL;
}
// 新しいAlphaをサンプリング
void HMM::sample_new_alpha(std::vector<std::vector<Word*>> &dataset){
	double new_alpha = sampler::normal(_alpha, 0.1 * _alpha);
	int random_index = sampler::uniform_int(0, dataset.size() - 1);
	std::vector<Word*> &line = dataset[random_index];
	// メトロポリス・ヘイスティングス法
	// http://ebsa.ism.ac.jp/ebooks/sites/default/files/ebook/1881/pdf/vol3_ch10.pdf
	// 提案分布は正規分布
	double log_Pt_alpha = compute_log_Pt_alpha(line, _alpha);
	double log_Pt_new_alpha = compute_log_Pt_alpha(line, new_alpha);
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
			random_word = _get_random_word_with_tag(tag, dataset);
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
		double Pti_wi_beta = compute_Pti_wi_beta(random_word->tag_id, random_word->word_id, beta);
		double Pti_wi_new_beta = compute_Pti_wi_beta(random_word->tag_id, random_word->word_id, new_beta);
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
		int count = get_count_for_tag_word(tag, word_id);
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
		std::cout << tag << ": " << get_word_types_for_tag(tag) << std::endl;
	}
}
bool HMM::save(std::string dir = "out"){
	std::ofstream ofs(dir + "/hmm.obj");
	boost::archive::binary_oarchive oarchive(ofs);
	oarchive << static_cast<const HMM&>(*this);
	ofs.close();
	std::ofstream ofs_bin;
	// 3-gram
	ofs_bin.open(dir + "/hmm.trigram", std::ios::binary);
	for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			ofs_bin.write((char*)(_trigram_counts[tri_tag][bi_tag]), _num_tags * sizeof(int));
		}
	}
	ofs_bin.close();
	// 2-gram
	ofs_bin.open(dir + "/hmm.bigram", std::ios::binary);
	for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
		ofs_bin.write((char*)(_bigram_counts[bi_tag]), _num_tags * sizeof(int));
	}
	ofs_bin.close();
	// 1-gram
	ofs_bin.open(dir + "/hmm.unigram", std::ios::binary);
	ofs_bin.write((char*)(_unigram_counts), _num_tags * sizeof(int));
	ofs_bin.close();
	// beta
	ofs_bin.open(dir + "/hmm.beta", std::ios::binary);
	ofs_bin.write((char*)(_beta), _num_tags * sizeof(double));
	ofs_bin.close();
	// Wt
	ofs_bin.open(dir + "/hmm.wt", std::ios::binary);
	ofs_bin.write((char*)(_Wt), _num_tags * sizeof(int));
	ofs_bin.close();
	return true;
}
bool HMM::load(std::string dir = "out"){
	bool complete = true;
	std::ifstream ifs(dir + "/hmm.obj");
	if(ifs.good()){
		boost::archive::binary_iarchive iarchive(ifs);
		iarchive >> *this;
	}else{
		complete = false;
	}
	ifs.close();

	// ポインタの読み込み
	if(_trigram_counts == NULL){
		alloc_table();
	}
	std::ifstream ifs_bin;
	// 3-gram
	ifs_bin.open(dir + "/hmm.trigram", std::ios::binary);
	if(ifs_bin.good()){
		for(int tri_tag = 0;tri_tag < _num_tags;tri_tag++){
			for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
				ifs_bin.read((char*)(_trigram_counts[tri_tag][bi_tag]), _num_tags * sizeof(int));
			}
		}
	}else{
		complete = false;
	}
	ifs_bin.close();
	// 2-gram
	ifs_bin.open(dir + "/hmm.bigram", std::ios::binary);
	if(ifs_bin.good()){
		for(int bi_tag = 0;bi_tag < _num_tags;bi_tag++){
			ifs_bin.read((char*)(_bigram_counts[bi_tag]), _num_tags * sizeof(int));
		}
	}else{
		complete = false;
	}
	ifs_bin.close();
	// 1-gram
	ifs_bin.open(dir + "/hmm.unigram", std::ios::binary);
	if(ifs_bin.good()){
		ifs_bin.read((char*)(_unigram_counts), _num_tags * sizeof(int));
	}else{
		complete = false;
	}
	ifs_bin.close();
	// beta
	ifs_bin.open(dir + "/hmm.beta", std::ios::binary);
	if(ifs_bin.good()){
		ifs_bin.read((char*)(_beta), _num_tags * sizeof(double));
	}else{
		complete = false;
	}
	ifs_bin.close();
	// Wt
	ifs_bin.open(dir + "/hmm.wt", std::ios::binary);
	if(ifs_bin.good()){
		ifs_bin.read((char*)(_Wt), _num_tags * sizeof(int));
	}else{
		complete = false;
	}
	ifs_bin.close();
	return complete;
}