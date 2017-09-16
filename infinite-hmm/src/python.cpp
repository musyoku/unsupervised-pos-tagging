#include "python/model.h"
#include "python/corpus.h"
#include "python/dataset.h"
#include "python/dictionary.h"
#include "python/trainer.h"

using namespace ihmm;

BOOST_PYTHON_MODULE(bhmm){
	boost::python::class_<Dictionary>("dictionary")
	.def("string_to_word_id", &Dictionary::string_to_word_id)
	.def("is_unk", &Dictionary::is_unk)
	.def("save", &Dictionary::save)
	.def("load", &Dictionary::load);

	boost::python::class_<Corpus>("corpus")
	.def("add_words", &Corpus::python_add_words);

	boost::python::class_<Dataset>("dataset", boost::python::init<Corpus*, double, int>())
	.def("get_num_words", &Dataset::get_num_words)
	.def("get_dict", &Dataset::get_dict_obj, boost::python::return_internal_reference<>());

	boost::python::class_<Trainer>("trainer", boost::python::init<Dataset*, Model*>())
	.def("compute_log_p_dataset_train", &Trainer::compute_log_p_dataset_train)
	.def("compute_log_p_dataset_dev", &Trainer::compute_log_p_dataset_dev)
	.def("update_hyperparameters", &Trainer::update_hyperparameters)
	.def("anneal_temperature", &Trainer::anneal_temperature)
	.def("perform_gibbs_sampling", &Trainer::perform_gibbs_sampling);

	boost::python::class_<Model>("model", boost::python::init<int, Dataset*>())
	.def(boost::python::init<std::string>())
	.def("get_num_tags", &Model::get_num_tags)
	.def("set_initial_alpha", &Model::set_initial_alpha)
	.def("set_initial_beta", &Model::set_initial_beta)
	.def("viterbi_decode", &Model::python_viterbi_decode)
	// .def("print_typical_words_assigned_to_each_tag", &Model::print_typical_words_assigned_to_each_tag)
	// .def("print_alpha_and_beta", &Model::print_alpha_and_beta)
	.def("save", &Model::save)
	.def("load", &Model::load);
}