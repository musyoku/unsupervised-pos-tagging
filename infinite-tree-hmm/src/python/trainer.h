#pragma once
#include <boost/python.hpp>
#include <cassert>
#include "dataset.h"
#include "model.h"
#include "dictionary.h"

namespace ithmm {
	class Trainer{
	public:
		Dataset* _dataset;
		Model* _model;
		Dictionary* _dict;
		std::vector<int> _rand_indices;
		double** _forward_table;		// 前向き確率計算用
		double** _decode_table;			// viterbiデコーディング用
		Trainer(Dataset* dataset, Model* model, Dictionary* dict);
		void remove_all_data();
		void perform_gibbs_sampling();
		void _before_viterbi_decode(std::vector<Node*> &nodes);
		void _after_viterbi_decode();
		void _before_compute_log_Pdataset(std::vector<Node*> &nodes);
		void _after_compute_log_Pdataset();
		double compute_log_Pdataset_train();
		double compute_log_Pdataset_dev();
		double _compute_log_Pdataset(std::vector<std::vector<Word*>> &dataset);
		double compute_log2_Pdataset_train();
		double compute_log2_Pdataset_dev();
		double _compute_log2_Pdataset(std::vector<std::vector<Word*>> &dataset);
		double compute_perplexity_train();
		double compute_perplexity_dev();
		double _compute_perplexity(std::vector<std::vector<Word*>> &dataset);
		void update_hyperparameters();
		void show_assigned_words_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag, bool show_probability = true);
	};
}