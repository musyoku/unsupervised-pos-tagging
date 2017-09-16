#pragma once
#include <boost/python.hpp>
#include <string>
#include "../ihmm/ihmm.h"
#include "dataset.h"
#include "dictionary.h"

namespace ihmm {
	class Model{
	private:
		void _alloc_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table);
		void _free_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table);
		void _set_locale();
	public:
		InfiniteHMM* _hmm;
		Model(int num_initial_tags, Dataset* dataset);
		Model(std::string filename);
		~Model();
		bool load(std::string filename);
		bool save(std::string filename);
		void set_initial_alpha(double alpha);
		void set_initial_beta(double beta);
		void set_initial_gamma(double gamma);
		void set_initial_gamma_emission(double gamma_emission);
		void set_initial_beta_emission(double beta_emission);
		int get_num_tags();
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence);
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double** forward_table, double** decode_table);
		boost::python::list python_viterbi_decode(boost::python::list py_word_ids);
		double compute_p_sentence(std::vector<Word*> &sentence, double** forward_table);
		void print_typical_words_assigned_to_each_tag(int number_to_show, Dictionary* dict);
		void print_alpha_and_beta();
	};
}