#include "model.h"

namespace bhmm {
	Model::Model(){
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_hmm = new HMM();
	}
	Model::~Model(){
		delete _hmm;
	}
	bool Model::load(std::string filename){
		return _hmm->load(filename);
	}
	bool Model::save(std::string filename){
		return _hmm->save(filename);
	}
	void Model::set_alpha(double alpha){
		_hmm->_alpha = alpha;
	}
	void Model::set_num_tags(int number){
		_hmm->_num_tags = number;
	}
	int Model::get_num_tags(){
		return _hmm->_num_tags;
	}
	double Model::get_temperature(){
		return _hmm->_temperature;
	}
	void Model::set_temperature(double temperature){
		_hmm->_temperature = temperature;
	}
	void Model::set_Wt(boost::python::list Wt){
		int length = boost::python::len(Wt);
		for(int tag = 0;tag < length;tag++){
			_hmm->set_Wt_for_tag(tag, boost::python::extract<int>(Wt[tag]));
		}
	}
	void Model::set_minimum_temperature(double temperature){
		_hmm->_minimum_temperature = temperature;
	}
	void Model::anneal_temperature(double temperature){
		_hmm->anneal_temperature(temperature);
	}
}