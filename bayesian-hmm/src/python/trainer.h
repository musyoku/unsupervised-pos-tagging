#pragma once
#include <boost/python.hpp>
#include "model.h"
#include "dataset.h"
#include "dictionary.h"

namespace bhmm {
	class Trainer{
	private:
		void _before_viterbi_decode();
		void _after_viterbi_decode();
		void _before_compute_log_p_dataset();
		void _after_compute_log_p_dataset();
		double _compute_log_p_dataset(std::vector<std::vector<Word*>> &dataset);
		double _compute_log2_p_dataset(std::vector<std::vector<Word*>> &dataset);
		double _compute_perplexity(std::vector<std::vector<Word*>> &dataset);
		Model* _model;
		Dictionary* _dict;
		Dataset* _dataset;
		std::vector<int> _rand_indices;
		double*** _forward_table;		// 前向き確率計算用
		double*** _decode_table;			// viterbiデコーディング用
	public:
		Trainer(Dataset* dataset, Model* model, boost::python::list py_Wt);
		Trainer(Dataset* dataset, Model* model, std::vector<int> &Wt);
		void perform_gibbs_sampling();
		void update_hyperparameters();
		boost::python::list python_get_all_words_of_each_tag(int threshold = 0);
		void show_typical_words_of_each_tag(int number_to_show);
		double compute_log_p_dataset_train();
		double compute_log_p_dataset_dev();
	};
}