#include "python/model.h"
#include "python/dataset.h"
#include "python/dictionary.h"
#include "python/trainer.h"

using namespace bhmm;

BOOST_PYTHON_MODULE(bhmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
	.def("is_unk", &Dictionary::is_unk)
	.def("save", &Dictionary::save)
	.def("load", &Dictionary::load);

	boost::python::class_<Dataset>("dataset")
	.def("get_num_words", &Dataset::get_num_words)
	.def("get_count_of_word", &Dataset::get_count_of_word)
	.def("get_dict", &Dataset::get_dict_obj, boost::python::return_internal_reference<>())
	.def("add_words_train", &Dataset::python_add_words_train)
	.def("add_words_dev", &Dataset::python_add_words_dev)
	.def("add_textfile", &Dataset::add_textfile)
	.def("mark_low_frequency_words_as_unknown", &Dataset::mark_low_frequency_words_as_unknown);
	
	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*>())
	.def("compute_log_p_dataset_train", &Trainer::compute_log_p_dataset_train)
	.def("compute_log_p_dataset_dev", &Trainer::compute_log_p_dataset_dev)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("anneal_temperature", &Trainer::anneal_temperature)
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling);

	boost::python::class_<Model>("model", boost::python::init<int, Dataset*, boost::python::list>())
	.def(boost::python::init<std::string>())
	.def("get_num_tags", &Model::get_num_tags)
	.def("get_temperature", &Model::get_temperature)
	.def("set_temperature", &Model::set_temperature)
	.def("set_minimum_temperature", &Model::set_minimum_temperature)
	.def("set_initial_alpha", &Model::set_initial_alpha)
	.def("set_initial_beta", &Model::set_initial_beta)
	.def("anneal_temperature", &Model::anneal_temperature)
	.def("viterbi_decode", &Model::python_viterbi_decode)
	.def("print_typical_words_assigned_to_each_tag", &Model::print_typical_words_assigned_to_each_tag)
	.def("print_alpha_and_beta", &Model::print_alpha_and_beta)
	.def("save", &Model::save)
	.def("load", &Model::load);
}