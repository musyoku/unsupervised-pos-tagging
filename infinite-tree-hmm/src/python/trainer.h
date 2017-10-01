#pragma once
#include <boost/python.hpp>
#include "model.h"
#include "dataset.h"
#include "dictionary.h"

namespace ithmm {
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
		double** _forward_table;		// 前向き確率計算用
		double** _decode_table;			// viterbiデコーディング用
	public:
		Trainer(Dataset* dataset, Model* model);
		void perform_gibbs_sampling();
		double compute_log_p_dataset_train();
		double compute_log_p_dataset_dev();
		void set_model(Model* model);
	};
}