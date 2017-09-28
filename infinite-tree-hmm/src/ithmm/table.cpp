#include "table.h"
#include "sampler.h"
#include <iostream>

namespace ithmm {
	Table::Table(){
		_num_customers = 0;
		_token_id = 0;
		_last_added_index = 0;
	}
	Table::Table(int token_id){
		_num_customers = 0;
		_token_id = token_id;
		_last_added_index = 0;
	}
	template <class Archive>
	void Table::serialize(Archive & ar, unsigned int version)
	{
		static_cast<void>(version);
		ar & _arrangement;
		ar & _num_customers;
		ar & _token_id;
		ar & _last_added_index;
	}
	template void Table::serialize(boost::archive::binary_iarchive &ar, unsigned int version);
	template void Table::serialize(boost::archive::binary_oarchive &ar, unsigned int version);
	bool Table::is_empty(){
		return _arrangement.size() == 0;
	}
	void Table::add_customer(double concentration_parameter, double g0, int num_total_customers, bool &new_table_generated){
		_num_customers += 1;
		if(_arrangement.size() == 0){
			_arrangement.push_back(1);
			new_table_generated = true;
			_last_added_index = _arrangement.size() - 1;
			return;
		}
		new_table_generated = false;
		double sum_arrangement = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
		// g0を掛けるかどうかがよく分からない
		// double total_counts = num_total_customers + concentration_parameter * g0;
		double total_counts = num_total_customers + concentration_parameter;
		double bernoulli = sampler::uniform(0, 1);
		// まず親から生成されたかどうかを判定
		if(bernoulli >= (num_total_customers / total_counts)){
			_arrangement.push_back(1);
			new_table_generated = true;
			_last_added_index = _arrangement.size() - 1;
			return;
		}
		// 自身のテーブルのどこに座るかを決定
		bernoulli = sampler::uniform(0, 1);
		double normalizer = 1.0 / sum_arrangement;
		double stack = 0;
		for(int i = 0;i < _arrangement.size();i++){
			stack += _arrangement[i] * normalizer;
			if(bernoulli <= stack){
				_arrangement[i] += 1;
				_last_added_index = i;
				return;
			}
		}
		assert(false);
	}
	void Table::remove_customer(bool &empty_table_deleted){
		assert(_arrangement.size() > 0);
		empty_table_deleted = false;
		_num_customers -= 1;
		int sum_arrangement = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
		double normalizer = 1.0 / sum_arrangement;
		int bernoulli = sampler::uniform_int(0, 1);
		double stack = 0;
		int target_index = _arrangement.size() - 1;
		for(int i = 0;i < _arrangement.size();i++){
			stack += _arrangement[i] * normalizer;
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
	void Table::remove_last_customer(bool &empty_table_deleted){
		empty_table_deleted = false;
		assert(_arrangement.size() > 0);
		assert(_last_added_index < _arrangement.size());
		_arrangement[_last_added_index] -= 1;
		if(_arrangement[_last_added_index] == 0){
			_arrangement.erase(_arrangement.begin() + _last_added_index);
			empty_table_deleted = true;
		}
	}
}
