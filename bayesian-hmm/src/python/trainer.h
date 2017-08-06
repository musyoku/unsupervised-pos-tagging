#pragma once
#include <boost/python.hpp>
#include "model.h"
#include "dataset.h"
#include "dictionary.h"

namespace bhmm {
	class Trainer{
	private:
		Model* _model;
		Dictionary* _dict;
		Dataset* _dataset;
		std::vector<int> _rand_indices;
	public:
		Trainer(Dataset* dataset, Model* model, Dictionary* dict);
		void perform_gibbs_sampling();
		void sample_new_alpha();
		void sample_new_beta();
		boost::python::list get_all_words_for_each_tag(int threshold = 0);
		void show_typical_words_for_each_tag(int number_to_show_for_each_tag);
	};
}