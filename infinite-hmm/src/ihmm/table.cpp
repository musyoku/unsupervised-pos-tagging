#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <numeric>
#include "table.h"
#include "sampler.h"

namespace ihmm {
	template <class Archive>
	void Table::serialize(Archive & ar, unsigned int version)
	{
		ar & _arrangement;
		ar & _num_customers;
		ar & _identifier;
	}
	template void Table::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive & ar, unsigned int version);
	template void Table::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive & ar, unsigned int version);
	Table::Table(){
		_num_customers = 0;
		_identifier = 0;
	}
	Table::Table(int identifier){
		_num_customers = 0;
		_identifier = identifier;
	}
	bool Table::is_empty(){
		return _arrangement.size() == 0;
	}
	void Table::add_customer(double concentration_parameter, bool &new_table_generated){
		_num_customers += 1;
		if(_arrangement.size() == 0){
			_arrangement.push_back(1);
			new_table_generated = true;
			return;
		}
		new_table_generated = false;
		double sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0) + concentration_parameter;
		double normalizer = 1 / sum;
		double bernoulli = sampler::uniform(0, 1);
		double stack = 0;
		for(int i = 0;i < _arrangement.size();i++){
			stack += _arrangement[i] * normalizer;
			if(bernoulli <= stack){
				_arrangement[i] += 1;
				return;
			}
		}
		_arrangement.push_back(1);
		new_table_generated = true;
	}
	void Table::remove_customer(bool &empty_table_deleted){
		assert(_arrangement.size() > 0);
		empty_table_deleted = false;
		_num_customers -= 1;
		int sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
		int bernoulli = sampler::uniform_int(0, sum);
		double stack = 0;
		int target_index = _arrangement.size() - 1;
		for(int i = 0;i < _arrangement.size();i++){
			stack += _arrangement[i];
			if(bernoulli <= stack){
				target_index = i;
				break;
			}
		}
		_arrangement[target_index] -= 1;
		if(_arrangement[target_index] == 0){
			_arrangement.erase(_arrangement.begin() + target_index);
			empty_table_deleted = true;
		}
	}
	int Table::get_num_customers(){
		return _num_customers;
	}
}