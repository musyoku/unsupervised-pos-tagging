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
	PyBayesianHMM* hmm = new PyBayesianHMM(5);
	hmm->load_textfile("../test.txt");
	hmm->initialize();
	for(int i = 0;i < 1000;i++){
		hmm->perform_gibbs_sampling();
	}
	return 0;
}