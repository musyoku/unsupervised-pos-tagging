#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <cassert>
#include <iostream>
#include "sampler.h"
#include "node.h"
#include "hpylm.h"

namespace ithmm {
	using std::vector;

	HPYLM::HPYLM(){
		_num_tables = 0;
		_num_customers = 0;
		_state_node = NULL;
		_depth = -1;
		_parent = NULL;
	}
	HPYLM::HPYLM(Node* node){
		assert(node != NULL);
		_num_tables = 0;
		_num_customers = 0;
		_state_node = node;
		_depth = node->_depth_v;
		_parent = NULL;
		if(node->_parent != NULL){
			assert(node->_parent->_hpylm != NULL);
			_parent = node->_parent->_hpylm;
		}
	}
	template <class Archive>
	void HPYLM::serialize(Archive &ar, unsigned int version)
	{
		ar & _arrangement;
		ar & _num_tables;
		ar & _num_customers;
		ar & _depth;
		ar & _parent;
		ar & _state_node;
	}
	template void HPYLM::serialize(boost::archive::binary_iarchive &ar, unsigned int version);
	template void HPYLM::serialize(boost::archive::binary_oarchive &ar, unsigned int version);
	// 客をテーブルに追加
	bool HPYLM::_add_customer_to_table(int token_id, int table_k, double g0, vector<double> &d_m, vector<double> &theta_m){
		auto itr = _arrangement.find(token_id);
		assert(itr != _arrangement.end());
		vector<int> &num_customers_at_table = itr->second;
		assert(table_k < num_customers_at_table.size());
		num_customers_at_table[table_k]++;
		_num_customers++;
		return true;
	}
	bool HPYLM::_add_customer_to_new_table(int token_id, double g0, vector<double> &d_m, vector<double> &theta_m){
		auto itr = _arrangement.find(token_id);
		if(itr == _arrangement.end()){
			vector<int> tables = {1};
			_arrangement[token_id] = tables;
		}else{
			vector<int> &tables = itr->second;
			tables.push_back(1);
		}
		_num_tables++;
		_num_customers++;
		if(_parent != NULL){
			bool success = _parent->add_customer(token_id, g0, d_m, theta_m);
			assert(success == true);
		}
		return true;
	}
	bool HPYLM::_remove_customer_from_table(int token_id, int table_k, vector<int> &num_customers_at_table){
		assert(table_k < num_customers_at_table.size());
		num_customers_at_table[table_k]--;
		_num_customers--;
		assert(_num_customers >= 0);
		assert(num_customers_at_table[table_k] >= 0);
		if(num_customers_at_table[table_k] == 0){
			if(_parent != NULL){
				bool success = _parent->remove_customer(token_id);
				assert(success == true);
			}
			num_customers_at_table.erase(num_customers_at_table.begin() + table_k);
			_num_tables--;
			assert(_num_tables >= 0);
			if(num_customers_at_table.size() == 0){
				_arrangement.erase(token_id);
			}
		}
		return true;
	}
	int HPYLM::get_num_tables_serving_word(int token_id){
		if(_arrangement.find(token_id) == _arrangement.end()){
			return 0;
		}
		return _arrangement[token_id].size();
	}
	int HPYLM::get_num_customers_eating_word(int token_id){
		auto itr = _arrangement.find(token_id);
		if(itr == _arrangement.end()){
			return 0;
		}
		vector<int> &num_customers_at_table = itr->second;
		int sum = 0;
		for(int i = 0;i < num_customers_at_table.size();i++){
			sum += num_customers_at_table[i];
		}
		return sum;
	}
	bool HPYLM::add_customer(int token_id, double g0, vector<double> &d_m, vector<double> &theta_m){
		assert(_depth < d_m.size());
		double d_u = d_m[_depth];
		double theta_u = theta_m[_depth];
		double parent_p_w = g0;
		if(_parent){
			parent_p_w = _parent->compute_p_w(token_id, g0, d_m, theta_m);
		}
		auto itr = _arrangement.find(token_id);
		if(itr == _arrangement.end()){
			_add_customer_to_new_table(token_id, g0, d_m, theta_m);
			return true;
		}
		vector<int> &num_customers_at_table = itr->second;
		double sum = 0.0;
		for(int k = 0;k < num_customers_at_table.size();k++){
			sum += std::max(0.0, num_customers_at_table[k] - d_u);
		}
		double t_u = _num_tables;
		sum += (theta_u + d_u * t_u) * parent_p_w;
		double normalizer = 1.0 / sum;
		double r = sampler::uniform(0, 1);
		double stack = 0;
		for(int k = 0;k < num_customers_at_table.size();k++){
			stack += std::max(0.0, num_customers_at_table[k] - d_u) * normalizer;
			if(r <= stack){
				_add_customer_to_table(token_id, k, g0, d_m, theta_m);
				return true;
			}
		}
		_add_customer_to_new_table(token_id, g0, d_m, theta_m);
		return true;
	}
	bool HPYLM::remove_customer(int token_id){
		auto itr = _arrangement.find(token_id);
		assert(itr != _arrangement.end());
		vector<int> &num_customers_at_table = itr->second;
		double sum = std::accumulate(num_customers_at_table.begin(), num_customers_at_table.end(), 0);
		double normalizer = 1.0 / sum;
		double r = sampler::uniform(0, 1);
		sum = 0;
		for(int k = 0;k < num_customers_at_table.size();k++){
			sum += num_customers_at_table[k] * normalizer;
			if(r <= sum){
				_remove_customer_from_table(token_id, k, num_customers_at_table);
				return true;
			}
		}
		_remove_customer_from_table(token_id, num_customers_at_table.size() - 1, num_customers_at_table);
		return true;
	}
	double HPYLM::compute_p_w(int token_id, double g0, vector<double> &d_m, vector<double> &theta_m){
		assert(_depth < d_m.size());
		assert(_depth < theta_m.size());
		double d_u = d_m[_depth];
		double theta_u = theta_m[_depth];
		double t_u = _num_tables;
		double c_u = _num_customers;
		auto itr = _arrangement.find(token_id);
		if(itr == _arrangement.end()){
			double coeff = (theta_u + d_u * t_u) / (theta_u + c_u);
			if(_parent != NULL){
				return _parent->compute_p_w(token_id, g0, d_m, theta_m) * coeff;
			}
			return g0 * coeff;
		}
		double parent_p_w = g0;
		if(_parent != NULL){
			parent_p_w = _parent->compute_p_w(token_id, g0, d_m, theta_m);
		}
		vector<int> &num_customers_at_table = itr->second;
		double c_uw = std::accumulate(num_customers_at_table.begin(), num_customers_at_table.end(), 0);
		double t_uw = num_customers_at_table.size();
		double first_term = std::max(0.0, c_uw - d_u * t_uw) / (theta_u + c_u);
		double coeff = (theta_u + d_u * t_u) / (theta_u + c_u);
		return first_term + coeff * parent_p_w;
	}
	int HPYLM::get_num_tables(){
		int num = 0;
		for(const auto &elem: _arrangement){
			num += elem.second.size();
		}
		assert(num == _num_tables);
		return num;
	}
	int HPYLM::get_num_customers(){
		int num = 0;
		for(const auto &elem: _arrangement){
			num += std::accumulate(elem.second.begin(), elem.second.end(), 0);
		}
		assert(num == _num_customers);
		return num;
	}
	double HPYLM::auxiliary_log_x_u(double theta_u){
		if(_num_customers >= 2){
			double x_u = sampler::beta(theta_u + 1, _num_customers - 1);
			return log(x_u + 1e-8);
		}
		return 0;
	}
	double HPYLM::auxiliary_y_ui(double d_u, double theta_u){
		if(_num_tables >= 2){
			double sum_y_ui = 0;
			for(int i = 1;i <= _num_tables - 1;i++){
				double denominator = theta_u + d_u * i;
				assert(denominator > 0);
				sum_y_ui += sampler::bernoulli(theta_u / denominator);;
			}
			return sum_y_ui;
		}
		return 0;
	}
	double HPYLM::auxiliary_1_y_ui(double d_u, double theta_u){
		if(_num_tables >= 2){
			double sum_1_y_ui = 0;
			for(int i = 1;i <= _num_tables - 1;i++){
				double denominator = theta_u + d_u * i;
				assert(denominator > 0);
				sum_1_y_ui += 1.0 - sampler::bernoulli(theta_u / denominator);
			}
			return sum_1_y_ui;
		}
		return 0;
	}
	double HPYLM::auxiliary_1_z_uwkj(double d_u){
		double sum_z_uwkj = 0;
		// c_u..
		for(auto &elem: _arrangement){
			// c_uw.
			vector<int> &num_customers_at_table = elem.second;
			for(int k = 0;k < num_customers_at_table.size();k++){
				// c_uwk
				int c_uwk = num_customers_at_table[k];
				if(c_uwk >= 2){
					for(int j = 1;j <= c_uwk - 1;j++){
						assert(j - d_u > 0);
						sum_z_uwkj += 1 - sampler::bernoulli((j - 1) / (j - d_u));
					}
				}
			}
		}
		return sum_z_uwkj;
	}
	void HPYLM::dump(){
		assert(_state_node != NULL);
		std::string indices_str = "";
		for(int i = 0;i < _state_node->_depth_v;i++){
			indices_str += std::to_string(_state_node->_horizontal_indices_from_root[i]);
			indices_str += ",";
		}
		std::cout << (boost::format("HPYLM: %d [tables:%d,customers:%d][%s]") % _state_node->_identifier % _num_tables % _num_customers % indices_str.c_str()).str() << std::endl;
	}
}