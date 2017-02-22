#include "hpylm.hpp"
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_map>
#include <map>
#include <cmath>
#include <cstdlib>
#include "cprintf.h"
#include "sampler.h"
#include "const.h"
#include "vocab.h"

// 客をテーブルに追加
bool HPYLM::add_customer_to_table(id &token_id, int table_k, double parent_Pw, vector<double> &d_m, vector<double> &theta_m, int &added_to_table_k){
	if(_arrangement.find(token_id) == _arrangement.end()){
		return add_customer_to_new_table(token_id, parent_Pw, d_m, theta_m, added_to_table_k);
	}
	vector<int> &num_customers_at_table = _arrangement[token_id];
	if(table_k < num_customers_at_table.size()){
		num_customers_at_table[table_k]++;
		_num_customers++;
		return true;
	}
	c_printf("[r]%s [*]%s\n", "エラー:", "客を追加できません. table_k < _arrangement[token_id].size()");
	exit(1);
	return false;
}
bool HPYLM::add_customer_to_new_table(id &token_id, double parent_Pw, vector<double> &d_m, vector<double> &theta_m, int &added_to_table_k){
	if(_arrangement.find(token_id) == _arrangement.end()){
		vector<int> tables = {1};
		_arrangement[token_id] = tables;
	}else{
		_arrangement[token_id].push_back(1);
	}
	_num_tables++;
	_num_customers++;
	if(_parent != NULL){
		bool success = _parent->add_customer(token_id, parent_Pw, d_m, theta_m, false, added_to_table_k);
		if(success == false){
			c_printf("[r]%s [*]%s\n", "エラー:", "客を追加できません. success == false");
			exit(1);
		}
	}
	return true;
}
bool HPYLM::remove_customer_from_table(id &token_id, int table_k, int &removed_from_table_k){
	if(_arrangement.find(token_id) == _arrangement.end()){
		c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. _arrangement.find(token_id) == _arrangement.end()");
		exit(1);
	}
	if(table_k >= _arrangement[token_id].size()){
		c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. table_k >= _arrangement[token_id].size()");
		exit(1);
	}
	vector<int> &num_customers_at_table = _arrangement[token_id];
	num_customers_at_table[table_k]--;
	_num_customers--;
	if(num_customers_at_table[table_k] < 0){
		c_printf("[r]%s [*]%s\n", "エラー:", "客の管理に不具合があります. num_customers_at_table[table_k] < 0");
		exit(1);
	}
	if(num_customers_at_table[table_k] == 0){
		if(_parent != NULL){
			bool success = _parent->remove_customer(token_id, false, removed_from_table_k);
			if(success == false){
				c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. success == false");
				exit(1);
			}
		}
		num_customers_at_table.erase(num_customers_at_table.begin() + table_k);
		_num_tables--;
		if(num_customers_at_table.size() == 0){
			_arrangement.erase(token_id);
		}
	}
	return true;
}
HPYLM::HPYLM(){
	_num_tables = 0;
	_num_customers = 0;
}
bool HPYLM::parent_exists(){
	return !(_parent == NULL);
}
bool HPYLM::child_exists(id &token_id){
	return !(_children.find(token_id) == _children.end());
}
bool HPYLM::need_to_remove_from_parent(){
	if(_parent == NULL){
		return false;
	}
	if(_children.size() == 0 and _arrangement.size() == 0){
		return true;
	}
	return false;
}
int HPYLM::get_num_tables_serving_word(id &token_id){
	if(_arrangement.find(token_id) == _arrangement.end()){
		return 0;
	}
	return _arrangement[token_id].size();
}
int HPYLM::get_num_customers_eating_word(id &token_id){
	if(_arrangement.find(token_id) == _arrangement.end()){
		return 0;
	}
	vector<int> &num_customers_at_table = _arrangement[token_id];
	int sum = 0;
	for(int i = 0;i < num_customers_at_table.size();i++){
		sum += num_customers_at_table[i];
	}
	return sum;
}
Node* HPYLM::find_child_node(id &token_id, bool generate_if_not_exist = false){
	auto itr = _children.find(token_id);
	if (itr != _children.end()) {
		return itr->second;
	}
	if(generate_if_not_exist == false){
		return NULL;
	}
	Node* child = new Node(token_id);
	child->_parent = this;
	child->_depth = _depth + 1;
	_children[token_id] = child;
	return child;
}
bool HPYLM::add_customer(id &token_id, double g0, vector<double> &d_m, vector<double> &theta_m, bool update_n, int &added_to_table_k){
	init_hyperparameters_at_depth_if_needed(_depth, d_m, theta_m);
	double d_u = d_m[_depth];
	double theta_u = theta_m[_depth];
	double parent_Pw = g0;
	if(_parent){
		parent_Pw = _parent->Pw(token_id, g0, d_m, theta_m);
	}
	if(_arrangement.find(token_id) == _arrangement.end()){
		add_customer_to_new_table(token_id, parent_Pw, d_m, theta_m, added_to_table_k);
		if(update_n == true){
			increment_stop_count();
		}
		if(_depth == 0){	// if root node
			added_to_table_k = 0;
		}
		return true;
	}
	vector<int> &num_customers_at_table = _arrangement[token_id];
	double sum_props = 0.0;
	for(int k = 0;k < num_customers_at_table.size();k++){
		sum_props += std::max(0.0, num_customers_at_table[k] - d_u);
	}
	double t_u = _num_tables;
	sum_props += (theta_u + d_u * t_u) * parent_Pw;
	double normalizer = 1.0 / sum_props;
	double r = Sampler::uniform(0, 1);
	double sum_normalized_probs = 0.0;
	for(int k = 0;k < num_customers_at_table.size();k++){
		sum_normalized_probs += std::max(0.0, num_customers_at_table[k] - d_u) * normalizer;
		if(r <= sum_normalized_probs){
			add_customer_to_table(token_id, k, parent_Pw, d_m, theta_m, added_to_table_k);
			if(update_n){
				increment_stop_count();
			}
			if(_depth == 0){
				added_to_table_k = k;
			}
			return true;
		}
	}
	add_customer_to_new_table(token_id, parent_Pw, d_m, theta_m, added_to_table_k);
	if(update_n){
		increment_stop_count();
	}
	if(_depth == 0){
		added_to_table_k = num_customers_at_table.size() - 1;
	}
	return true;
}
bool HPYLM::remove_customer(id &token_id, bool update_n, int &removed_from_table_k){
	if(_arrangement.find(token_id) == _arrangement.end()){
		c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. _arrangement.find(token_id) == _arrangement.end()");
		exit(1);
	}
	vector<int> &num_customers_at_table = _arrangement[token_id];
	double sum_props = std::accumulate(num_customers_at_table.begin(), num_customers_at_table.end(), 0);		
	double normalizer = 1.0 / sum_props;
	double r = Sampler::uniform(0, 1);
	double sum_normalized_probs = 0.0;
	for(int k = 0;k < num_customers_at_table.size();k++){
		sum_normalized_probs += num_customers_at_table[k] * normalizer;
		if(r <= sum_normalized_probs){
			remove_customer_from_table(token_id, k, removed_from_table_k);
			if(update_n == true){
				decrement_stop_count();
			}
			if(_depth == 0){
				removed_from_table_k = k;
			}
			return true;
		}
	}
	remove_customer_from_table(token_id, num_customers_at_table.size() - 1, removed_from_table_k);
	if(update_n == true){
		decrement_stop_count();
	}
	if(_depth == 0){
		removed_from_table_k = num_customers_at_table.size() - 1;
	}
	return true;
}
double HPYLM::Pw(id &token_id, double g0, vector<double> &d_m, vector<double> &theta_m){
	init_hyperparameters_at_depth_if_needed(_depth, d_m, theta_m);
	double d_u = d_m[_depth];
	double theta_u = theta_m[_depth];
	double t_u = _num_tables;
	double c_u = _num_customers;
	auto itr = _arrangement.find(token_id);
	if(itr == _arrangement.end()){
		double coeff = (theta_u + d_u * t_u) / (theta_u + c_u);
		if(_parent != NULL){
			return _parent->Pw(token_id, g0, d_m, theta_m) * coeff;
		}
		return g0 * coeff;
	}
	double parent_Pw = g0;
	if(_parent != NULL){
		parent_Pw = _parent->Pw(token_id, g0, d_m, theta_m);
	}
	vector<int> &num_customers_at_table = itr->second;
	double c_uw = std::accumulate(num_customers_at_table.begin(), num_customers_at_table.end(), 0);
	double t_uw = num_customers_at_table.size();
	double first_coeff = std::max(0.0, c_uw - d_u * t_uw) / (theta_u + c_u);
	double second_coeff = (theta_u + d_u * t_u) / (theta_u + c_u);
	return first_coeff + second_coeff * parent_Pw;
}
double HPYLM::stop_probability(double beta_stop, double beta_pass){
	double p = (_stop_count + beta_stop) / (_stop_count + _pass_count + beta_stop + beta_pass);
	if(_parent != NULL){
		p *= _parent->pass_probability(beta_stop, beta_pass);
	}
	return p;
}
double HPYLM::pass_probability(double beta_stop, double beta_pass){
	double p = (_pass_count + beta_pass) / (_stop_count + _pass_count + beta_stop + beta_pass);
	if(_parent != NULL){
		p *= _parent->pass_probability(beta_stop, beta_pass);
	}
	return p;
}
void HPYLM::increment_stop_count(){
	_stop_count++;
	if(_parent != NULL){
		_parent->increment_pass_count();
	}
}
void HPYLM::decrement_stop_count(){
	_stop_count--;
	if(_stop_count < 0){
		c_printf("[r]%s [*]%s\n", "エラー:", "停止回数の管理に不具合があります. _stop_count < 0");
		exit(1);
	}
	if(_parent != NULL){
		_parent->decrement_passC_count();
	}
}
void HPYLM::increment_pass_count(){
	_pass_count++;
	if(_parent != NULL){
		_parent->increment_pass_count();
	}
}
void HPYLM::decrement_passC_count(){
	_pass_count--;
	if(_pass_count < 0){
		c_printf("[r]%s [*]%s\n", "エラー:", "通過回数の管理に不具合があります. _pass_count < 0");
		exit(1);
	}
	if(_parent != NULL){
		_parent->decrement_passC_count();
	}
}
bool HPYLM::remove_from_parent(){
	if(_parent == NULL){
		return false;
	}
	_parent->delete_child_node(_token_id);
	return true;
}
void HPYLM::delete_child_node(id &token_id){
	Node* child = find_child_node(token_id);
	if(child){
		_children.erase(token_id);
		// __children->delete_key(token_id);
		delete child;
	}
	if(_children.size() == 0 && _arrangement.size() == 0){
		remove_from_parent();
	}
}
int HPYLM::get_num_nodes(){
	int num = _children.size();
	for(const auto &elem: _children){
		num += elem.second->get_num_nodes();
	}
	return num;
}
int HPYLM::get_num_tables(){
	int num = 0;
	for(const auto &elem: _arrangement){
		num += elem.second.size();
	}
	if(num != _num_tables){
		c_printf("[r]%s [*]%s\n", "エラー:", "テーブルの管理に不具合があります. num != _num_tables");
		exit(1);
	}
	for(const auto &elem: _children){
		num += elem.second->get_num_tables();
	}
	return num;
}
int HPYLM::get_num_customers(){
	int num = 0;
	for(const auto &elem: _arrangement){
		num += std::accumulate(elem.second.begin(), elem.second.end(), 0);
	}
	if(num != _num_customers){
		c_printf("[r]%s [*]%s\n", "エラー:", "客の管理に不具合があります. num != _num_customers");
		exit(1);
	}
	for(const auto &elem: _children){
		num += elem.second->get_num_customers();
	}
	return num;
}
// dとθの推定用
// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
double HPYLM::auxiliary_log_x_u(double theta_u){
	if(_num_customers >= 2){
		double x_u = Sampler::beta(theta_u + 1, _num_customers - 1);
		return log(x_u + 1e-8);
	}
	return 0;
}
double HPYLM::auxiliary_y_ui(double d_u, double theta_u){
	if(_num_tables >= 2){
		double sum_y_ui = 0;
		for(int i = 1;i <= _num_tables - 1;i++){
			double denominator = theta_u + d_u * i;
			if(denominator == 0){
				c_printf("[r]%s [*]%s\n", "エラー:", "0除算です. denominator == 0");
				exit(1);
			}
			sum_y_ui += Sampler::bernoulli(theta_u / denominator);;
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
			if(denominator == 0){
				c_printf("[r]%s [*]%s\n", "エラー:", "0除算です. denominator == 0");
				exit(1);
			}
			sum_1_y_ui += 1.0 - Sampler::bernoulli(theta_u / denominator);
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
					if(j - d_u == 0){
						c_printf("[r]%s [*]%s\n", "エラー:", "0除算です. j - d_u == 0");
						exit(1);
					}
					sum_z_uwkj += 1 - Sampler::bernoulli((j - 1) / (j - d_u));
				}
			}
		}
	}
	return sum_z_uwkj;
}