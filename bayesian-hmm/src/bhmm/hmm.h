#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <set>
#include "common.h"

namespace bhmm {
	class HMM{
	private:
		void alloc_count_tables(int num_tags);
		void init_ngram_counts_with_corpus(std::vector<std::vector<Word*>> &dataset);
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
		bool _allocated;
		HMM();
		~HMM();
		void anneal_temperature(double multiplier);
		void initialize_with_training_corpus(std::vector<std::vector<Word*>> &dataset, std::vector<int> &Wt);
		void increment_tag_ngram_count(Word* tri_word, Word* bi_word, Word* uni_word);
		void increment_tag_word_count(int tag_id, int word_id);
		void decrement_tag_word_count(int tag_id, int word_id);
		int get_count_for_tag_word(int tag_id, int word_id);
		int get_word_types_for_tag(int tag_id);
		Word* _get_random_word_with_tag(int tag, std::vector<std::vector<Word*>> &dataset);
		int get_most_co_occurring_tag(int word_id);
		void set_Wt_for_tag(int tag_id, int number);
		void set_num_tags(int n);
		double compute_log_Pt_alpha(std::vector<Word*> &word_vec, double alpha);
		double compute_log_Pw_t_alpha(std::vector<Word*> &word_vec, double alpha);
		double compute_Pti_wi_beta(int ti, int wi, double beta);
		double compute_p_wi_given_ti_beta(int wi, int ti, double beta);
		double compute_p_ti_given_t_alpha(int ti, int ti_1, int ti_2, double alpha);
		void add_tags_to_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi);
		void remove_tags_from_model_parameters(int ti_2, int ti_1, int ti, int ti1, int ti2, int wi);
		void perform_gibbs_sampling_with_sequence(std::vector<Word*> &word_vec);
		int sample_tag_from_Pt_w(int ti_2, int ti_1, int wi);
		int argmax_tag_from_Pt_w(int ti_2, int ti_1, int wi);
		void sample_new_alpha(std::vector<std::vector<Word*>> &dataset);
		void sample_new_beta(std::vector<std::vector<Word*>> &dataset);
		void dump_trigram_counts();
		void dump_bigram_counts();
		void dump_unigram_counts();
		void dump_word_types();
		bool save(std::string filename);
		bool load(std::string filename);
		template <class Archive>
		void serialize(Archive& archive, unsigned int version);
	};
}

namespace boost { 
	namespace serialization {
		template<class Archive>
		void save(Archive &ar, const bhmm::HMM &hmm, unsigned int version);
		template<class Archive>
		void load(Archive &ar, bhmm::HMM &hmm, unsigned int version);
	}
}