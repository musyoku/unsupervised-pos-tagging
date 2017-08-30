#include "python/model.h"
#include "python/dataset.h"
#include "python/dictionary.h"
#include "python/trainer.h"

using namespace bhmm;

BOOST_PYTHON_MODULE(bhmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
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
	
	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*, boost::python::list>())
	.def("compute_log_p_dataset_dev", &Trainer::compute_log_p_dataset_dev)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("show_typical_words_of_each_tag", &Trainer::show_typical_words_of_each_tag)
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling);

	boost::python::class_<Model>("model", boost::python::init<int>())
	.def("get_temperature", &Model::get_temperature)
	.def("set_temperature", &Model::set_temperature)
	.def("set_minimum_temperature", &Model::set_minimum_temperature)
	.def("anneal_temperature", &Model::anneal_temperature)
	.def("save", &Model::save)
	.def("load", &Model::load);
}