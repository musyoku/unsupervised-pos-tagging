#pragma once
#include <boost/python.hpp>
#include <string>
#include "../bhmm/hmm.h"

namespace bhmm {
	class Model{
	private:
		void _alloc_viterbi_tables(int sentence_length, double*** forward_table, double*** decode_table);
		void _free_viterbi_tables(int sentence_length, double*** forward_table, double*** decode_table);
	public:
		HMM* _hmm;
		Model(int num_tags);
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
		void anneal_temperature(double temperature);
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence);
		void viterbi_decode(std::vector<Word*> &sentence, std::vector<int> &sampled_state_sequence, double*** forward_table, double*** decode_table);
		boost::python::list python_viterbi_decode(boost::python::list py_word_ids);
		double compute_p_sentence(std::vector<Word*> &sentence, double*** forward_table);
		
	};
}