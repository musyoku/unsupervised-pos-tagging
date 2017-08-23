#include "python/model.h"
#include "python/dataset.h"
#include "python/dictionary.h"
#include "python/trainer.h"

using namespace ithmm;

BOOST_PYTHON_MODULE(ithmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
	.def("get_vocab_size", &Dictionary::get_vocab_size)
	.def("save", &Dictionary::save)
	.def("load", &Dictionary::load)
	.def("add_word", &Dictionary::add_word_string);

	boost::python::class_<Dataset>("dataset")
	.def("get_num_words", &Dataset::get_num_words)
	.def("get_count_of_word", &Dataset::get_count_of_word)
	.def("get_dict", &Dataset::get_dict, boost::python::return_internal_reference<>())
	.def("add_words_train", &Dataset::python_add_words_train)
	.def("add_words_dev", &Dataset::python_add_words_dev)
	.def("add_textfile", &Dataset::add_textfile)
	.def("mark_low_frequency_words_as_unknown", &Dataset::mark_low_frequency_words_as_unknown);
	
	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*, Dictionary*>())
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("compute_perplexity_dev", &Trainer::compute_perplexity_dev)
	.def("compute_perplexity_train", &Trainer::compute_perplexity_train)
	.def("compute_log2_p_dataset_dev", &Trainer::compute_log2_p_dataset_dev)
	.def("compute_log2_p_dataset_train", &Trainer::compute_log2_p_dataset_train)
	.def("compute_log_p_dataset_dev", &Trainer::compute_log_p_dataset_dev)
	.def("compute_log_p_dataset_train", &Trainer::compute_log_p_dataset_train)
	.def("show_assigned_words_for_each_tag", &Trainer::show_assigned_words_for_each_tag)
	.def("remove_all_data", &Trainer::remove_all_data);

	boost::python::class_<Model>("model")
	.def("viterbi_decode", &Model::python_viterbi_decode)
	.def("update_hyperparameters", &Model::update_hyperparameters)
	.def("save", &Model::save)
	.def("load", &Model::load)
	.def("get_alpha", &Model::get_alpha)
	.def("get_gamma", &Model::get_gamma)
	.def("get_lambda_alpha", &Model::get_lambda_alpha)
	.def("get_lambda_gamma", &Model::get_lambda_gamma)
	.def("get_strength", &Model::get_strength)
	.def("get_tau0", &Model::get_tau0)
	.def("get_tau1", &Model::get_tau1)
	.def("get_metropolis_hastings_acceptance_rate", &Model::get_metropolis_hastings_acceptance_rate)
	.def("get_all_states", &Model::python_get_all_states)
	.def("set_alpha", &Model::set_alpha)
	.def("set_gamma", &Model::set_gamma)
	.def("set_lambda_alpha", &Model::set_lambda_alpha)
	.def("set_lambda_gamma", &Model::set_lambda_gamma)
	.def("set_strength", &Model::set_strength)
	.def("set_tau0", &Model::set_tau0)
	.def("set_tau1", &Model::set_tau1)
	.def("set_depth_limit", &Model::set_depth_limit)
	.def("set_metropolis_hastings_enabled", &Model::set_metropolis_hastings_enabled)
	.def("show_hpylm_for_each_tag", &Model::show_hpylm_for_each_tag)
	.def("show_sticks", &Model::show_sticks)
	.def("show_assigned_words_for_each_tag", &Model::show_assigned_words_for_each_tag)
	.def("show_assigned_words_and_probability_for_each_tag", &Model::show_assigned_words_and_probability_for_each_tag);
}