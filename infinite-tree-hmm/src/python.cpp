#include "model.h"
#include "dataset.h"
#include "dictionary.h"
#include "trainer.h"

BOOST_PYTHON_MODULE(ithmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
	.def("save", &Dictionary::save)
	.def("load", &Dictionary::load)
	.def("add_string", &Dictionary::add_string);

	boost::python::class_<Dataset>("dataset", boost::python::init<Dictionary*>())
	.def("get_num_words", &Dataset::get_num_words)
	.def("get_count_for_word", &Dataset::get_count_for_word)
	.def("add_train_data", &Dataset::add_train_data)
	.def("add_test_data", &Dataset::add_test_data)
	.def("add_textfile", &Dataset::add_textfile)
	.def("mark_low_frequency_words_as_unknown", &Dataset::mark_low_frequency_words_as_unknown);
	
	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*, Dictionary*>())
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("viterbi_decode_train", &Trainer::viterbi_decode_train)
	.def("viterbi_decode_test", &Trainer::viterbi_decode_test)
	.def("compute_perplexity_test", &Trainer::compute_perplexity_test)
	.def("compute_perplexity_train", &Trainer::compute_perplexity_train)
	.def("compute_log2_Pdataset_test", &Trainer::compute_log2_Pdataset_test)
	.def("compute_log2_Pdataset_train", &Trainer::compute_log2_Pdataset_train)
	.def("compute_log_Pdataset_test", &Trainer::compute_log_Pdataset_test)
	.def("compute_log_Pdataset_train", &Trainer::compute_log_Pdataset_train)
	.def("show_assigned_words_for_each_tag", &Trainer::show_assigned_words_for_each_tag)
	.def("remove_all_data", &Trainer::remove_all_data);

	boost::python::class_<Model>("model")
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
	.def("get_all_tags", &Model::get_all_tags)
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