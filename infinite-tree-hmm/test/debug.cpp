#include  <iostream>
#include  <chrono>
#include <algorithm>
#include "../src/ithmm/ithmm.h"
#include "../src/ithmm/cprintf.h"
#include "../src/python/model.h"
#include "../src/python/dataset.h"
#include "../src/python/dictionary.h"
#include "../src/python/trainer.h"
using namespace std;
using namespace ithmm;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->sample_node_in_htssb(model->_structure_tssb->_root->_transition_tssb);
		model->add_customer_to_htssb_node(node);
	}
}

void test1(iTHMM* model){
	add_customer(model, 100);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
	}
}

void test2(iTHMM* model){
	add_customer(model, 20);
	Node* target_in_structure = model->_structure_tssb->find_node_with_id(11);
	assert(target_in_structure != NULL);
	for(int n = 0;n < 100;n++){
		model->add_customer_to_htssb_node(target_in_structure);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = target_in_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}

	for(int n = 0;n < 100;n++){
		model->remove_customer_from_htssb_node(target_in_structure);
	}

	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	parent = target_in_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test6(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node_in_htssb(model->_structure_tssb->_root->_transition_tssb);
		nodes.push_back(node);
	}
	int rand_index = sampler::uniform_int(0, nodes.size());
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to_htssb_node(node);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to_htssb_node(node);
		}
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_customer_from_htssb_node(node);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_customer_from_htssb_node(node);
		}
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	c_printf("[*]%s\n", "transition");
	model->_structure_tssb->_root->_transition_tssb->dump();
}

void test7(iTHMM* model){
	add_customer(model, 20);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target_in_structure = model->_structure_tssb->find_node_with_id(21);
	assert(target_in_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_in_structure->dump();
	target_in_structure->_transition_tssb_myself->dump();
	for(int i = 0;i < 10;i++){
		model->add_customer_to_htssb_node(target_in_structure->_transition_tssb_myself);
	}
	c_printf("[*]%s\n", "transition");
	target_in_structure->_transition_tssb->dump();
	Node* parent_in_structure = target_in_structure->_parent;
	while(parent_in_structure){
		c_printf("[*]%d\n", parent_in_structure->_identifier);
		parent_in_structure->_transition_tssb->dump();
		parent_in_structure = parent_in_structure->_parent;
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_vertical_sbr_ratio(target_in_structure->_transition_tssb_myself, true);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;
}

void test8(iTHMM* model){
	add_customer(model, 20);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target_in_structure = model->_structure_tssb->find_node_with_id(21);
	assert(target_in_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_in_structure->dump();
	target_in_structure->_transition_tssb_myself->dump();
	for(int i = 0;i < 10;i++){
		model->add_customer_to_htssb_node(target_in_structure->_transition_tssb_myself);
	}
	c_printf("[*]%s\n", "transition");
	target_in_structure->_transition_tssb->dump();
	Node* parent_in_structure = target_in_structure->_parent;
	while(parent_in_structure){
		c_printf("[*]%d\n", parent_in_structure->_identifier);
		parent_in_structure->_transition_tssb->dump();
		parent_in_structure = parent_in_structure->_parent;
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_horizontal_sbr_ratio(target_in_structure->_transition_tssb_myself, true);
	}
	auto end_time = chrono::system_clock::now();
	auto duration = end_time - start_time;
	auto msec = chrono::duration_cast<chrono::milliseconds>(duration).count();
	cout << ratio << endl;
	cout << msec << " msec" << endl;
}


void test9(iTHMM* model){
	add_customer(model, 10);
	double uniform = 0;
	TSSB* tssb = model->_structure_tssb->_root->_transition_tssb;
	model->update_stick_length_of_tssb(tssb, 1.0, true);
	tssb->dump();
	vector<Node*> nodes_true;
	vector<Node*> nodes;
	int num_nodes_true = tssb->get_num_nodes();
	int num_nodes = 0;
	int prev_id = -1;
	for(int i = 0;i < 100000;i++){
		uniform = i / 100000.0;
		Node* node = model->retrospective_sampling(uniform, tssb, 1.0, true);
		assert(node != NULL);
		if(node->_identifier != prev_id){
			cout << uniform << ": " << node->_identifier << endl;
			prev_id = node->_identifier;
			num_nodes += 1;
			nodes.push_back(node);
		}
	}
	tssb->enumerate_nodes_from_left_to_right(nodes_true);
	int limit = (nodes_true.size() > nodes.size()) ? nodes_true.size() : nodes.size();
	for(int i = 0;i < limit;i++){
		int identifier = (i < nodes.size()) ? nodes[i]->_identifier : -1;
		int identifier_true = (i < nodes_true.size()) ? nodes_true[i]->_identifier : -1;
		cout << identifier << " : " << identifier_true << endl;
	}
	cout << num_nodes << " == " << num_nodes_true << endl;
	assert(num_nodes == num_nodes_true);
}

void test10(iTHMM* model){
	add_customer(model, 100);
	double uniform = 0.91;
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(7);
	model->compute_expectation_of_horizontal_sbr_ratio(target->_transition_tssb_myself, true);
	assert(target != NULL);
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
	model->update_stick_length_of_tssb(target->_transition_tssb, 1.0, true);
	Node* node = model->retrospective_sampling(uniform, target->_transition_tssb, 1.0, true);
	c_printf("[*]%s\n", "sampled");
	node->dump();
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test11(iTHMM* model){
	for(int n = 0;n < 10;n++){
		Node* node = model->sample_node_in_htssb(model->_structure_tssb->_root->_transition_tssb);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = model->_structure_tssb->find_node_with_id(12);
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		Node* myself = parent->_transition_tssb_myself;
		myself->dump();
		if(myself->_identifier != 1){
			assert(myself->_parent_transition_tssb_myself != NULL);
			myself->_parent_transition_tssb_myself->dump();
		}
		assert(myself->_structure_tssb_myself != NULL);
		myself->_structure_tssb_myself->dump();
		parent = parent->_parent;
	}
}

void test12(iTHMM* model){
	for(int n = 0;n < 10;n++){
		Node* node = model->sample_node_in_htssb(model->_structure_tssb->_root->_transition_tssb);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(12);
	model->delete_node_in_structure_if_needed(target);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = model->_structure_tssb->find_node_with_id(14);
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		Node* myself = parent->_transition_tssb_myself;
		myself->dump();
		if(myself->_identifier != 1){
			assert(myself->_parent_transition_tssb_myself != NULL);
			myself->_parent_transition_tssb_myself->dump();
		}
		assert(myself->_structure_tssb_myself != NULL);
		myself->_structure_tssb_myself->dump();
		parent = parent->_parent;
	}
}

void test13(iTHMM* model){
	add_customer(model, 1000);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(7);
	assert(target != NULL);
	Table* table = target->_table_v;
	int n = 1;
	double alpha = 10;
	double g0 = 0.5;
	cout << (n / (n + alpha)) / ((n + alpha * g0) / (n + alpha)) << " : " << (alpha * g0 / (n + alpha)) / ((n + alpha * g0) / (n + alpha)) << endl;
	int generated = 0;
	int total = 0;
	for(int i = 0;i < 100000;i++){
		bool new_table_generated = false;
		table->add_customer(alpha, g0, n, new_table_generated);
		total++;
		if(new_table_generated){
			generated++;
		}
	}
	cout << ((total - generated) / (double)total) << " : " << (generated / (double)total) << endl;
}

void test14(iTHMM* model){
	add_customer(model, 10);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(7);
	for(int i = 0;i < 10000;i++){
		model->add_customer_to_htssb_node(target->_transition_tssb_myself);
	}
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
	for(int i = 0;i < 10000;i++){
		model->remove_customer_from_htssb_node(target->_transition_tssb_myself);
	}
	parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test15(iTHMM* model){
	add_customer(model, 10000);
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		vector<Node*> nodes_in_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_in_htssb);
		for(const auto node_in_htssb: nodes_in_htssb){
			for(int i = 0;i < 1000;i++){
				model->add_customer_to_htssb_node(node_in_htssb);
			}
		}
	}
	for(const auto node: nodes){
		vector<Node*> nodes_in_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_in_htssb);
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
		for(const auto node_in_htssb: nodes_in_htssb){
			node_in_htssb->dump();
			double ratio_v = model->compute_expectation_of_vertical_sbr_ratio(node_in_htssb, true);
			double ratio_h = model->compute_expectation_of_horizontal_sbr_ratio(node_in_htssb, true);
			cout << ratio_v << ", " << ratio_h << endl;
		}
	}
}

void test16(iTHMM* model){
	add_customer(model, 1000);
	c_printf("[*]%s\n", "structure");
	model->set_word_g0(1.0 / 100.0);
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(60);
	for(int i = 0;i < 100000;i++){
		int token_id = sampler::uniform_int(0, 100);
		model->add_customer_to_hpylm(target, token_id);
	}
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		double pw = model->compute_Pw_given_s(5, parent);
		cout << pw << endl;
		assert(parent->_hpylm != NULL);
		parent->_hpylm->dump();
		parent = parent->_parent;
	}
}

void test17(iTHMM* model){
	string dir = "out";
	test16(model);
	model->save(dir);
	iTHMM* copy = new iTHMM();
	copy->load(dir);
	c_printf("[*]%s\n", "structure");
	copy->_structure_tssb->dump();
	Node* parent = copy->_structure_tssb->find_node_with_id(60);
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		double pw = copy->compute_Pw_given_s(5, parent);
		cout << pw << endl;
		assert(parent->_hpylm != NULL);
		parent->_hpylm->dump();
		parent = parent->_parent;
	}
	delete copy;
}

void test18(iTHMM* model){
	string dir = "out";
	test15(model);
	model->save(dir);
	iTHMM* copy = new iTHMM();
	copy->load(dir);

	c_printf("[*]%s\n", "structure");
	copy->_structure_tssb->dump();
	vector<Node*> nodes;
	copy->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		vector<Node*> nodes_in_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_in_htssb);
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
		for(const auto node_in_htssb: nodes_in_htssb){
			node_in_htssb->dump();
			double ratio_v = copy->compute_expectation_of_vertical_sbr_ratio(node_in_htssb, true);
			double ratio_h = copy->compute_expectation_of_horizontal_sbr_ratio(node_in_htssb, true);
			cout << ratio_v << ", " << ratio_h << endl;
		}
	}
	delete copy;
}
void _test19(iTHMM* model, TSSB* tssb){
	model->update_stick_length_of_tssb(tssb, 1.0, true);
	c_printf("[*]%d\n", tssb->_owner_id);
	tssb->dump();
	vector<Node*> nodes;
	tssb->enumerate_nodes_from_left_to_right(nodes);
	double prev_probability = 0;
	for(const auto node: nodes){
		double probability = node->_sum_probability;
		assert(probability > prev_probability);
		prev_probability = probability;
		cout << probability << endl;
	}
}

void test19(iTHMM* model){
	string dir = "out";
	test15(model);
	model->save(dir);
	iTHMM* copy = new iTHMM();
	copy->load(dir);
	c_printf("[*]%s\n", "structure");
	copy->_structure_tssb->dump();
	vector<Node*> nodes;
	copy->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		vector<Node*> nodes_in_htssb;
		_test19(model, node->_transition_tssb);
	}
	delete copy;
}

void test20(iTHMM* model){
	string dir = "out";

	add_customer(model, 10000);
	c_printf("[*]%s\n", "bos");
	model->_bos_tssb->dump();

	model->save(dir);
	iTHMM* copy = new iTHMM();
	copy->load(dir);

	c_printf("[*]%s\n", "bos");
	copy->_bos_tssb->dump();

	for(int n = 0;n < 10000;n++){
		Node* node = copy->sample_node_in_tssb(copy->_bos_tssb);
		copy->add_customer_to_tssb_node(node);
	}

	c_printf("[*]%s\n", "bos");
	copy->_bos_tssb->dump();
	delete copy;
}

void test21(){
	string filename = "../alice.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->get_dict();
	dictionary->save("ithmm.dict");
	Model* model = new Model();
	Trainer* trainer = new Trainer(dataset, model, dictionary);
	dataset->mark_low_frequency_words_as_unknown(1);
	model->update_hyperparameters();
	model->_ithmm->_structure_tssb->dump();
	delete dictionary;
	delete dataset;
	delete model;
	delete trainer;
}

void test22(){
	iTHMM* model = new iTHMM();
	Node* node = model->sample_node_in_tssb(model->_structure_tssb);
	model->add_customer_to_htssb_node(node->_transition_tssb_myself);

	c_printf("[*]%s\n", "before");
	model->_structure_tssb->dump();
	model->_bos_tssb->dump();
	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		node->_transition_tssb->dump();
	}

	model->remove_customer_from_htssb_node(node->_transition_tssb_myself);

	c_printf("[*]%s\n", "after");
	model->_structure_tssb->dump();
	model->_bos_tssb->dump();
	for(const auto node: nodes){
		node->_transition_tssb->dump();
	}
	delete model;
}

void test23(){
	iTHMM* model = new iTHMM();
	vector<Node*> nodes = {model->_structure_tssb->_root};
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node_in_htssb(nodes[0]->_transition_tssb);
		node->dump();
		assert(node->_structure_tssb_myself);
		nodes.push_back(node->_structure_tssb_myself);
		std::random_shuffle(nodes.begin(), nodes.end());
	}
	int rand_index = sampler::uniform_int(0, nodes.size());
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to_htssb_node(node->_transition_tssb_myself);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->add_customer_to_htssb_node(node->_transition_tssb_myself);
		}
	}
	c_printf("[*]%s\n", "before");
	model->_structure_tssb->dump();
	model->_bos_tssb->dump();

	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_customer_from_htssb_node(node->_transition_tssb_myself);
		}
	}
	std::random_shuffle(nodes.begin(), nodes.end());
	for(auto node: nodes){
		for(int i = 0;i < 1000;i++){
			model->remove_customer_from_htssb_node(node->_transition_tssb_myself);
		}
	}

	c_printf("[*]%s\n", "after");
	model->_structure_tssb->dump();
	model->_bos_tssb->dump();
	delete model;
}

void test24(iTHMM* model){
	test15(model);
	Node* target = model->_structure_tssb->find_node_with_id(7);
	model->update_stick_length_of_tssb(target->_transition_tssb, 1.0, true);
	Node* node = model->retrospective_sampling(0.9, target->_transition_tssb, 1.0, true);
	c_printf("[*]%s\n", "sampled");
	node->dump();
	target->_transition_tssb->dump();
}

void test25(iTHMM* model){
	test15(model);
	Node* target = model->_structure_tssb->find_node_with_id(55);
	TSSB* tssb = target->_transition_tssb;
	target = tssb->find_node_by_tracing_horizontal_indices(target);
	double p = model->compute_node_probability_in_tssb(tssb, target, 1);
	cout << p << endl;
	model->update_stick_length_of_tssb(tssb, 1.0, true);
	tssb->dump();
}

void test26(iTHMM* model){
	test15(model);
	vector<Node*> nodes;
	model->update_stick_length_of_tssb(model->_structure_tssb->_root->_transition_tssb, 1.0, true);
	model->_structure_tssb->_root->_transition_tssb->enumerate_nodes_from_left_to_right(nodes);
	cout << nodes.size() << endl;
	for(int i = 0;i < nodes.size();i++){
		for(int j = 0;j < nodes.size();j++){
			if(i == j){
				continue;
			}
			Node* left = nodes[i];
			Node* right = nodes[j];
			if(left->_identifier == right->_identifier){
				continue;
			}
			if(left->_depth_v == 0){
				continue;
			}
			if(right->_depth_v == 0){
				continue;
			}
			bool flag = model->is_node_to_the_left_of_node(left, right);
			assert((flag && (left->_sum_probability < right->_sum_probability)) || (flag == false && (left->_sum_probability > right->_sum_probability)));
		}
	}
}

void test27(){
	string filename = "../alice.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->get_dict();
	dictionary->save("ithmm.dict");
	Model* model = new Model();
	Trainer* trainer = new Trainer(dataset, model, dictionary);

	dataset->add_textfile(filename, 1000);
	dataset->mark_low_frequency_words_as_unknown(1);
	model->update_hyperparameters();
	c_printf("[*]%s\n", "structure");
	model->_ithmm->_structure_tssb->dump();

	Node* prev_state = model->_ithmm->_structure_tssb->find_node_with_id(20);
	Node* state = model->_ithmm->_structure_tssb->find_node_with_id(14);
	Node* next_state = model->_ithmm->_structure_tssb->find_node_with_id(31);
	Node* new_state = model->_ithmm->draw_state(prev_state, state, next_state, 10);
	new_state->dump();
	new_state = model->_ithmm->_draw_state_from_bos(state, next_state, 10);
	new_state->dump();
	new_state = model->_ithmm->_draw_state_to_eos(prev_state, state, 10);
	new_state->dump();

	trainer->remove_all_data();
	model->_ithmm->delete_invalid_children();
	c_printf("[*]%s\n", "structure");
	model->_ithmm->_structure_tssb->dump();

	delete dictionary;
	delete dataset;
	delete model;
	delete trainer;
}

void test28(){
	iTHMM* model = new iTHMM();
	for(int n = 0;n < 100;n++){
		Node* node = model->sample_node_in_htssb(model->_structure_tssb->_root->_transition_tssb);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();

	Node* target_in_structure = model->_structure_tssb->find_node_with_id(20);
	assert(target_in_structure != NULL);
	Node* target_in_htssb = target_in_structure->_transition_tssb_myself;
	assert(target_in_htssb != NULL);
	for(int n = 0;n < 1000;n++){
		model->add_customer_to_htssb_node(target_in_htssb);
	}
	c_printf("[*]%s\n", "htssb");
	target_in_structure->_transition_tssb->dump();

	vector<Node*> nodes;
	model->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	for(const auto node: nodes){
		int true_count = node->_transition_tssb->get_num_customers();
		int count = node->_transition_tssb->_num_customers;
		cout << count << " == " << true_count << endl;
		node->_transition_tssb->dump();
		assert(count == true_count);
	}
	delete model;
}

void test29(){
	string filename = "../alice.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->get_dict();
	dictionary->save("ithmm.dict");
	Model* model = new Model();
	Trainer* trainer = new Trainer(dataset, model, dictionary);

	model->set_depth_limit(1);
	dataset->add_textfile(filename, 0.8);

	string dir = "out";
	// model->mark_low_frequency_words_as_unknown(1);
	model->show_assigned_words_for_each_tag(dictionary, 20, false);
	model->set_metropolis_hastings_enabled(false);
	cout << dataset->get_num_words() << " words" << endl;
	for(int i = 0;i < 100;i++){
		trainer->perform_gibbs_sampling();
		trainer->update_hyperparameters();
		// model->save(dir);
		// c_printf("[*]%s\n", "structure");
		cout << "\repoch:" << i << flush;
		if(i % 100 == 0){
			cout << "epoch:" << i << endl;
			// model->_ithmm->_structure_tssb->dump();
			model->show_assigned_words_for_each_tag(dictionary, 20, false);
			model->save(dir);
			cout << "alpha: " << model->_ithmm->_alpha << endl;
			cout << "gamma: " << model->_ithmm->_gamma << endl;
			cout << "lambda_alpha: " << model->_ithmm->_lambda_alpha << endl;
			cout << "lambda_gamma: " << model->_ithmm->_lambda_gamma << endl;
			cout << "strength: " << model->_ithmm->_strength << endl;
			cout << "log_Pdata: " << trainer->compute_log_p_dataset_dev() << ", " << trainer->compute_log_p_dataset_train() << endl;
			cout << "PPL: " << trainer->compute_perplexity_dev() << ", " << trainer->compute_perplexity_train() << endl;
			cout << "MH: " << model->_ithmm->_num_mh_acceptance / (double)(model->_ithmm->_num_mh_acceptance + model->_ithmm->_num_mh_rejection) << endl;;
			for(int i = 0;i <= model->_ithmm->_current_max_depth;i++){
				cout << "d[" << i << "] = " << model->_ithmm->_hpylm_d_m[i] << endl;
				cout << "theta[" << i << "] = " << model->_ithmm->_hpylm_theta_m[i] << endl;
			}
		}
		// if(i == 4000){
		// 	model->set_depth_limit(2);
		// }
	}
	model->_ithmm->_structure_tssb->dump();
	trainer->remove_all_data();
	model->_ithmm->_structure_tssb->dump();

	delete model;
	delete trainer;
	delete dataset;
	delete dictionary;
}

void test30(iTHMM* model){
	add_customer(model, 100000);
	for(int i = 0;i < 1000;i++){
		double uniform = sampler::uniform(0, 1);
		Node* node = model->retrospective_sampling(uniform, model->_structure_tssb->_root->_transition_tssb, 1.0, true);
		model->add_customer_to_htssb_node(node);
	}
	model->update_stick_length_of_tssb(model->_structure_tssb->_root->_transition_tssb, 0.5, true);
	model->_structure_tssb->_root->_transition_tssb->dump();
	// Node* target = model->_structure_tssb->find_node_with_id(22);
	// cout << model->compute_node_probability_in_tssb(model->_structure_tssb, target, 1.0) << endl;
}

void test31(){
	string filename = "../alice.txt";
	string dir = "out";

	
	for(int i = 0;i < 5;i++){
		Dataset* dataset = new Dataset();
		dataset->add_textfile(filename, 0.8);
		Dictionary* dictionary = dataset->get_dict();
		dictionary->save("ithmm.dict");
		Model* model = new Model();
		Trainer* trainer = new Trainer(dataset, model, dictionary);

		dataset->add_textfile(filename, 0.8);
		dataset->mark_low_frequency_words_as_unknown(1);
		cout << dataset->get_num_words() << " words" << endl;
		for(int i = 0;i < 10;i++){
			trainer->perform_gibbs_sampling();
			trainer->update_hyperparameters();
		}
		trainer->remove_all_data();
		delete trainer;
		delete dataset;
		delete dictionary;
		delete model;
	}
}

void test32(){
	bool new_table_generated = true;
	for(int i = 0;i < 10;i++){
		Node* node = new Node();
		for(int j = 0;j < 10;j++){
			node->add_customer_to_vertical_crp(0.5, 0.5, new_table_generated);
			node->add_customer_to_horizontal_crp(0.5, 0.5, new_table_generated);
			node->increment_word_assignment(j);
		}
		delete node;
	}
}
void test33(){
	bool new_table_generated = true;
	for(int i = 0;i < 10;i++){
		Node* node = new Node();
		TSSB* tssb = new TSSB(node);
		for(int j = 0;j < 10;j++){
			node->add_customer_to_vertical_crp(0.5, 0.5, new_table_generated);
			node->add_customer_to_horizontal_crp(0.5, 0.5, new_table_generated);
			node->increment_word_assignment(j);
		}
		delete tssb;
	}
}
void test34(){
	bool new_table_generated = true;
	for(int i = 0;i < 10;i++){
		iTHMM* model = new iTHMM();
		for(int j = 0;j < 100;j++){
			Node* node = model->sample_node_in_tssb(model->_structure_tssb);
			model->add_customer_to_tssb_node(node);
			node = model->sample_node_in_tssb(model->_bos_tssb);
			model->add_customer_to_tssb_node(node);
		}
		delete model;
	}
}

void test35(){
	Dictionary* dictionary = new Dictionary();
	dictionary->load("ithmm.dict");
	Model* model = new Model();

	model->load("ithmm.model");
	model->show_hpylm_for_each_tag(dictionary);
	model->show_assigned_words_for_each_tag(dictionary, 20, true);
	model->show_sticks();
	for(int i = 0;i <= model->_ithmm->_current_max_depth;i++){
		cout << "d[" << i << "] = " << model->_ithmm->_hpylm_d_m[i] << endl;
		cout << "theta[" << i << "] = " << model->_ithmm->_hpylm_theta_m[i] << endl;
	}
	delete model;
	delete dictionary;
}

void test36(){
	Model* model = new Model();
	Dictionary* dictionary = new Dictionary();
	dictionary->load("ithmm.dict");
	model->_ithmm->set_word_g0(0.001);
	for(int i = 0;i < 10000;i++){
		double uniform = sampler::uniform(0, 1);
		Node* node = model->_ithmm->retrospective_sampling(uniform, model->_ithmm->_structure_tssb, 1.0, false);
		if(node->_depth_v == 0){
			continue;
		}
		model->_ithmm->add_customer_to_hpylm(node, 100);
	}
	c_printf("[*]%s\n", "sampled");
	model->show_hpylm_for_each_tag(dictionary);
	delete model;
}

void test37(){
	string filename = "../alice.txt";
	Dataset* dataset = new Dataset();
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->get_dict();
	Model* model = new Model();
	Trainer* trainer = new Trainer(dataset, model, dictionary);
	model->load("ithmm.model");
	model->show_assigned_words_for_each_tag(dictionary, 20, false);
	cout << dataset->get_num_words() << " words" << endl;
	cout << "alpha: " << model->_ithmm->_alpha << endl;
	cout << "gamma: " << model->_ithmm->_gamma << endl;
	cout << "lambda_alpha: " << model->_ithmm->_lambda_alpha << endl;
	cout << "strength: " << model->_ithmm->_strength << endl;
	cout << "log_Pdata: " << trainer->compute_log_p_dataset_dev() << endl;
	cout << "PPL: " << trainer->compute_perplexity_dev() << endl;
	for(int i = 0;i <= model->_ithmm->_current_max_depth;i++){
		cout << "d[" << i << "] = " << model->_ithmm->_hpylm_d_m[i] << endl;
		cout << "theta[" << i << "] = " << model->_ithmm->_hpylm_theta_m[i] << endl;
	}
	delete model;
	delete dictionary;
	delete trainer;
	delete dataset;
}

void test38(){
	string filename = "../alice.txt";
	Dataset* dataset = new Dataset();
	Model* model = new Model();
	model->set_depth_limit(1);
	dataset->add_textfile(filename, 0.8);
	Dictionary* dictionary = dataset->get_dict();
	Trainer* trainer = new Trainer(dataset, model, dictionary);

	// model->mark_low_frequency_words_as_unknown(1);
	// model->set_metropolis_hastings_enabled(true);
	for(int i = 0;i < 1000000;i++){
		trainer->perform_gibbs_sampling();
		trainer->update_hyperparameters();
		if(i % 100 == 0){
			cout << "log_Pdata: " << trainer->compute_log_p_dataset_train() << ", " << trainer->compute_log_p_dataset_dev() << endl;
			cout << "PPL: " << trainer->compute_perplexity_train() << ", " << trainer->compute_perplexity_dev() << endl;
		}
		// if(i == 4000){
		// 	model->set_depth_limit(2);
		// }
	}
}


int main(){
	iTHMM* hmm = new iTHMM();
	test1(hmm);
	test2(hmm);
	test6(hmm);
	test7(hmm);
	test8(hmm);
	test9(hmm);
	test10(hmm);
	test11(hmm);
	test12(hmm);
	test13(hmm);
	test14(hmm);
	test15(hmm);
	test16(hmm);
	test17(hmm);
	test18(hmm);
	test19(hmm);
	test20(hmm);
	test21();
	test22();
	test23();
	test24(hmm);
	test25(hmm);
	test26(hmm);
	test27();
	test28();
	test29();
	test30(hmm);
	test31();
	test32();
	test33();
	test34();
	test35();
	test36();
	test37();
	test38();
	return 0;
}