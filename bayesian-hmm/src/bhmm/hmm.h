#ifndef _bhmm_
#define _bhmm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <set>
#include "cprintf.h"
#include "sampler.h"
#include "utils.h"

using id = int;

typedef struct Word {
	id word_id;
	int tag_id;
} Word;

class HMM{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version);
public:
	int _num_tags;			// 品詞数
	int _num_words;			// 単語数
	int*** _trigram_counts;	// 品詞3-gramのカウント
	int** _bigram_counts;	// 品詞2-gramのカウント
	int* _unigram_counts;	// 品詞1-gramのカウント
	int* _Wt;
	std::unordered_map<int, std::unordered_map<int, int>> _tag_word_counts;	// 品詞と単語のペアの出現頻度
	double* _sampling_table;	// キャッシュ
	double _alpha;
	double* _beta;
	double _temperature;
	double _minimum_temperature;
	HMM();
	~HMM();
	void anneal_temperature(double multiplier);
	void initialize(std::vector<std::vector<Word*>> &dataset);
	void alloc_table();
	void init_ngram_counts(std::vector<std::vector<Word*>> &dataset);
	void set_Wt_for_tag(int tag_id, int number);
	void update_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word);
	void increment_tag_word_count(int tag_id, int word_id);
	void decrement_tag_word_count(int tag_id, int word_id);
	int get_count_for_tag_word(int tag_id, int word_id);
	int get_word_types_for_tag(int tag_id);
	double compute_log_Pt_alpha(std::vector<Word*> &line, double alpha);
	double compute_log_Pw_t_alpha(std::vector<Word*> &line, double alpha);
	double compute_Pti_wi_beta(int ti, int wi, double beta);
	void add_tag_to_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi);
	void remove_tag_from_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi);
	void perform_gibbs_sampling_with_line(std::vector<Word*> &line);
	int sample_tag_from_Pt_w(int ti_2, int ti_1, int wi);
	int argmax_tag_from_Pt_w(int ti_2, int ti_1, int wi);
	Word* _get_random_word_with_tag(int tag, std::vector<std::vector<Word*>> &dataset);
	void sample_new_alpha(std::vector<std::vector<Word*>> &dataset);
	void sample_new_beta(std::vector<std::vector<Word*>> &dataset);
	int get_most_co_occurring_tag(int word_id);
	void dump_trigram_counts();
	void dump_bigram_counts();
	void dump_unigram_counts();
	void dump_word_types();
	bool save(std::string dir);
	bool load(std::string dir);
};

#endif