#include "table.hpp"
#include "sampler.h"
#include <iostream>
using namespace std;

Table::Table(){
	_num_customers = 0;
	_token_id = 0;
}
Table::Table(int token_id){
	_num_customers = 0;
	_token_id = token_id;
}
bool Table::is_empty(){
	return _arrangement.size() == 0;
}
void Table::add_customer(double concentration_parameter, double g0, int num_total_customers, bool &new_table_generated){
	_num_customers += 1;
	if(_arrangement.size() == 0){
		_arrangement.push_back(1);
		new_table_generated = true;
		return;
	}
	new_table_generated = false;
	double sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
	double scale = num_total_customers / sum;
	sum = num_total_customers + concentration_parameter * g0;
	double normalizer = 1.0 / sum;
	double bernoulli = Sampler::uniform(0, 1);
	sum = 0;
	for(int i = 0;i < _arrangement.size();i++){
		sum += _arrangement[i] * scale * normalizer;
		if(bernoulli <= sum){
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
	int bernoulli = Sampler::uniform_int(0, sum);
	sum = 0;
	int target_index = _arrangement.size() - 1;
	for(int i = 0;i < _arrangement.size();i++){
		sum += _arrangement[i];
		if(bernoulli <= sum){
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
