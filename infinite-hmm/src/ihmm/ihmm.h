#pragma once
#include <boost/serialization/serialization.hpp>
#include <vector>
#include "common.h"
#include "table.h"

// <s>と</s>のIDは0

namespace ihmm {
	class InfiniteHMM {
	public:
		int _initial_num_tags;
		int _num_words;
		int _prev_num_tags;
		std::vector<std::vector<Table*>> _n_ij_tables;	// 品詞bigramの出現頻度
		std::vector<Table**> _m_iq_tables;	// 品詞と単語のペアの出現頻度
		int* _oracle_m_q_counts;
		int _oracle_sum_n_over_j;			// \sum_j{n_j^o}
		int _oracle_sum_m_over_q;			// \sum_j{m_q^o}
		std::vector<int> _sum_n_i_over_j;	// \sum_j{n_ij}の計算用
		std::vector<int> _oracle_n_j_counts;
		std::vector<int> _sum_m_i_over_q;
		double* _gibbs_sampling_table;
		double _alpha;
		double _beta;
		double _gamma;
		double _beta_emission;
		double _gamma_emission;
		InfiniteHMM();
		InfiniteHMM(int initial_num_tags, int num_words);
		~InfiniteHMM();
		void initialize_with_training_dataset(std::vector<std::vector<Word*>> &dataset);
		void _remove_all_training_dataset(std::vector<std::vector<Word*>> &dataset);
		int get_num_tags() const;
		int get_num_valid_tags() const;
		int get_num_words() const;
		int get_sum_n_i_over_j(int tag);
		int get_n_ij(int context_tag, int tag);
		int get_oracle_sum_n_over_j();
		int get_oracle_n_j(int tag);
		int get_sum_m_i_over_q(int tag);
		int get_m_iq(int tag, int word_id);
		int get_oracle_sum_m_over_q();
		int get_oracle_m_q(int word_id);
		int _get_new_tag();
		bool is_tag_new(int tag);
		double compute_p_tag_given_context(int tag, int context_tag);
		double compute_p_word_given_tag(int word_id, int tag);
		int _add_new_tag();
		void _delete_tag(int tag);
		void _increment_tag_bigram_count(int context_tag, int tag);
		void _decrement_tag_bigram_count(int context_tag, int tag);
		void _increment_tag_word_count(int tag, int word_id);
		void _decrement_tag_word_count(int tag, int word_id);
		void _increment_oracle_tag_count(int tag);
		void _decrement_oracle_tag_count(int tag);
		void _increment_oracle_word_count(int word_id);
		void _decrement_oracle_word_count(int word_id);
		int _perform_gibbs_sampling_on_markov_blanket(int ti_1, int ti1, int wi);
		void perform_gibbs_sampling_with_sequence(std::vector<Word*> &word_vec);
		bool save(std::string filename);
		bool load(std::string filename);
		template <class Archive>
		void serialize(Archive &ar, unsigned int version);
	};
}

namespace boost { 
	namespace serialization {
		template<class Archive>
		void save(Archive &ar, const ihmm::InfiniteHMM &hmm, unsigned int version);
		template<class Archive>
		void load(Archive &ar, ihmm::InfiniteHMM &hmm, unsigned int version);
	}
}