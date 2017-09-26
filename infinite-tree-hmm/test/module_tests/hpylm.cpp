#include  <iostream>
#include  <cassert>
#include  <string>
#include "../../src/ithmm/hpylm.h"
#include "../../src/ithmm/node.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;
using std::vector;

void test_add_customer(){
	Node* node = new Node();
	HPYLM* hpylm = new HPYLM(node);
	vector<double> d_m;
	d_m.push_back(0.2);
	vector<double> theta_m;
	theta_m.push_back(1);
	double g0 = 0.001;
	for(id token_id = 0;token_id < 100;token_id++){
		hpylm->add_customer(token_id, g0, d_m, theta_m);
	}
	double p_w_1 = hpylm->compute_p_w(0, g0, d_m, theta_m);
	double p_w_2 = hpylm->compute_p_w(100, g0, d_m, theta_m);
	assert(p_w_1 > p_w_2);
}

void test_remove_customer(){
	Node* node = new Node();
	HPYLM* hpylm = new HPYLM(node);
	vector<double> d_m;
	d_m.push_back(0.2);
	vector<double> theta_m;
	theta_m.push_back(1);
	double g0 = 0.001;
	for(int repeat = 0;repeat < 10;repeat++){
		for(int n = 0;n < 100;n++){
			for(id token_id = 0;token_id < 100;token_id++){
				hpylm->add_customer(token_id, g0, d_m, theta_m);
			}
		}
		assert(hpylm->_num_customers == 10000);
		for(int n = 0;n < 100;n++){
			for(id token_id = 0;token_id < 100;token_id++){
				hpylm->remove_customer(token_id);
			}
		}
		assert(hpylm->_num_customers == 0);
		for(id token_id = 0;token_id < 100;token_id++){
			double p_w = hpylm->compute_p_w(token_id, g0, d_m, theta_m);
			assert(p_w == g0);
		}
	}
}

void test_parent_child(){
	Node* parent_node = new Node();
	HPYLM* parent_hpylm = new HPYLM(parent_node);
	Node* child_node = new Node();
	HPYLM* child_hpylm = new HPYLM(child_node);
	child_hpylm->_parent = parent_hpylm;
	vector<double> d_m;
	d_m.push_back(0.2);
	vector<double> theta_m;
	theta_m.push_back(1);
	double g0 = 0.001;
	for(int repeat = 0;repeat < 10;repeat++){
		for(int n = 0;n < 100;n++){
			for(id token_id = 0;token_id < 100;token_id++){
				child_hpylm->add_customer(token_id, g0, d_m, theta_m);
			}
		}
		assert(child_hpylm->_num_customers == 10000);
		assert(parent_hpylm->_num_customers == child_hpylm->_num_tables);

		double p_w_1 = child_hpylm->compute_p_w(0, g0, d_m, theta_m);
		double p_w_2 = child_hpylm->compute_p_w(100, g0, d_m, theta_m);
		assert(p_w_1 > p_w_2);
		p_w_1 = parent_hpylm->compute_p_w(0, g0, d_m, theta_m);
		p_w_2 = parent_hpylm->compute_p_w(100, g0, d_m, theta_m);
		assert(p_w_1 > p_w_2);
		for(int n = 0;n < 100;n++){
			for(id token_id = 0;token_id < 100;token_id++){
				child_hpylm->remove_customer(token_id);
			}
		}
		assert(child_hpylm->_num_customers == 0);
		assert(child_hpylm->_num_tables == 0);
		assert(parent_hpylm->_num_customers == 0);
		assert(parent_hpylm->_num_tables == 0);
	}

}

int main(){
	test_add_customer();
	cout << "OK" << endl;
	test_remove_customer();
	cout << "OK" << endl;
	test_parent_child();
	cout << "OK" << endl;
	return 0;
}