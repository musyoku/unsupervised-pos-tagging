#include <iostream>
#include <cassert>
#include <string>
#include "../../src/ithmm/table.h"
using namespace ithmm;
using std::cout;
using std::flush;
using std::endl;
using std::vector;

void test_add_customer(){
	double concentration_parameter = 1;
	double g0 = 0.001;
	bool new_table_generated = false;
	Table* table = new Table();
	for(int i = 0;i < 100;i++){
		table->add_customer(concentration_parameter, g0, i * 2, new_table_generated);
	}
	assert(table->_num_customers == 100);
	delete table;
}

void test_remove_customer(){
	double concentration_parameter = 1;
	double g0 = 0.001;
	bool new_table_generated = false;
	Table* table = new Table();
	for(int i = 0;i < 100;i++){
		table->add_customer(concentration_parameter, g0, i * 2, new_table_generated);
	}
	for(int i = 0;i < 100;i++){
		table->remove_customer(new_table_generated);
	}
	assert(table->_num_customers == 0);
	assert(table->_arrangement.size() == 0);
	delete table;
}

int main(){
	test_add_customer();
	cout << "OK" << endl;
	test_remove_customer();
	cout << "OK" << endl;
	return 0;
}