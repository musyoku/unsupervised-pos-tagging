#include "python/model.h"
#include "python/dataset.h"
#include "python/dictionary.h"
#include "python/trainer.h"

using namespace bhmm;

BOOST_PYTHON_MODULE(bhmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
	.def("get_eos_id", &Dictionary::get_eos_id)
	.def("save", &Dictionary::save)
	.def("load", &Dictionary::load);

	boost::python::class_<Dataset>("dataset", boost::python::init<Dictionary*>())
	.def("get_num_words", &Dataset::get_num_words)
	.def("get_count_of_word", &Dataset::get_count_of_word)
	.def("add_words_train", &Dataset::python_add_words_train)
	.def("add_words_dev", &Dataset::python_add_words_dev)
	.def("add_textfile", &Dataset::add_textfile)
	.def("mark_low_frequency_words_as_unknown", &Dataset::mark_low_frequency_words_as_unknown);
	
	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*, Dictionary*, boost::python::list>())
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling);

	boost::python::class_<Model>("model")
	.def("set_temperature", &Model::set_temperature)
	.def("set_minimum_temperature", &Model::set_minimum_temperature)
	.def("save", &Model::save)
	.def("load", &Model::load);
}