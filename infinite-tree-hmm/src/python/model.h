#pragma once
#include <boost/python.hpp>
#include <string>
#include "../ithmm/ithmm.h"
#include "dataset.h"
#include "dictionary.h"

namespace ithmm {
	class Model{
	private:
		void _alloc_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table);
		void _free_viterbi_tables(int sentence_length, double*** &forward_table, double*** &decode_table);
		void _set_locale();
	public:
		iTHMM* _ithmm;
		Model(Dataset* dataset);
		Model(std::string filename);
		~Model();
		bool load(std::string filename);
		bool save(std::string filename);
		void set_initial_alpha(double alpha);
		void set_initial_beta(double beta);
		int get_num_tags();
		double get_temperature();
		void set_temperature(double temperature);
		void set_minimum_temperature(double temperature);
		void anneal_temperature(double decay);
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence);
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double*** forward_table, double*** decode_table);
		boost::python::list python_viterbi_decode(boost::python::list py_word_ids);
		double compute_p_sentence(std::vector<Word*> &sentence, double*** forward_table);
		void print_typical_words_assigned_to_each_tag(int number_to_show, Dictionary* dict);
		void print_alpha_and_beta();
	};
}