#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <string>
#include <set>
#include <unordered_map>
#include <functional>
#include <fstream>
#include <cassert>
#include "core/ihmm.h"
#include "core/util.h"
#include "model.cpp"
using namespace std;
using namespace boost;

int main(){
	PyBayesianHMM* hmm = new PyBayesianHMM(20);
	hmm->load_textfile("../alice.txt");
	hmm->initialize();
	for(int i = 0;i < 1000;i++){
		hmm->perform_gibbs_sampling();
		hmm->_hmm->dump_oracle_tags();
		// hmm->_hmm->dump_oracle_words();
		hmm->_hmm->check_oracle_tag_count();
		hmm->_hmm->check_oracle_word_count();
	}
	return 0;
}