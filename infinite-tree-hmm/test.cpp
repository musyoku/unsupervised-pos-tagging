#include  <iostream>
#include  <chrono>
#include <algorithm>
#include "core/ithmm.h"
#include "core/cprintf.h"
#include "model.cpp"
using namespace std;

void add_customer(iTHMM* model, int count){
	for(int n = 0;n < count;n++){
		Node* node = model->sample_node_on_htssb(model->_structure_tssb->_root->_transition_tssb);
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
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(11);
	assert(target_on_structure != NULL);
	for(int n = 0;n < 100;n++){
		model->add_customer_to_htssb_node(target_on_structure);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* parent = target_on_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}

	for(int n = 0;n < 100;n++){
		model->remove_customer_from_htssb_node(target_on_structure);
	}

	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	parent = target_on_structure;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		parent->_transition_tssb->dump();
		parent = parent->_parent;
	}
}

void test6(iTHMM* model){
	vector<Node*> nodes;
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node_on_htssb(model->_structure_tssb->_root->_transition_tssb);
		nodes.push_back(node);
	}
	int rand_index = Sampler::uniform_int(0, nodes.size());
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
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(21);
	assert(target_on_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_on_structure->dump();
	target_on_structure->_transition_tssb_myself->dump();
	for(int i = 0;i < 10;i++){
		model->add_customer_to_htssb_node(target_on_structure->_transition_tssb_myself);
	}
	c_printf("[*]%s\n", "transition");
	target_on_structure->_transition_tssb->dump();
	Node* parent_on_structure = target_on_structure->_parent;
	while(parent_on_structure){
		c_printf("[*]%d\n", parent_on_structure->_identifier);
		parent_on_structure->_transition_tssb->dump();
		parent_on_structure = parent_on_structure->_parent;
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_vertical_sbr_ratio(target_on_structure->_transition_tssb_myself, true);
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
	Node* target_on_structure = model->_structure_tssb->find_node_with_id(21);
	assert(target_on_structure != NULL);
	c_printf("[*]%s\n", "target");
	target_on_structure->dump();
	target_on_structure->_transition_tssb_myself->dump();
	for(int i = 0;i < 10;i++){
		model->add_customer_to_htssb_node(target_on_structure->_transition_tssb_myself);
	}
	c_printf("[*]%s\n", "transition");
	target_on_structure->_transition_tssb->dump();
	Node* parent_on_structure = target_on_structure->_parent;
	while(parent_on_structure){
		c_printf("[*]%d\n", parent_on_structure->_identifier);
		parent_on_structure->_transition_tssb->dump();
		parent_on_structure = parent_on_structure->_parent;
	}
	double ratio = 0;
	auto start_time = chrono::system_clock::now();
	for(int i = 0;i < 100000;i++){
		ratio = model->compute_expectation_of_horizontal_sbr_ratio(target_on_structure->_transition_tssb_myself, true);
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
	model->update_stick_length_of_tssb(tssb);
	tssb->dump();
	vector<Node*> nodes_true;
	vector<Node*> nodes;
	int num_nodes_true = tssb->get_num_nodes();
	int num_nodes = 0;
	int prev_id = -1;
	for(int i = 0;i < 100000;i++){
		uniform = i / 100000.0;
		Node* node = model->retrospective_sampling_on_htssb(uniform, tssb);
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
	model->update_stick_length_of_tssb(target->_transition_tssb);
	Node* node = model->retrospective_sampling_on_htssb(uniform, target->_transition_tssb);
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
		Node* node = model->sample_node_on_htssb(model->_structure_tssb->_root->_transition_tssb);
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
		Node* node = model->sample_node_on_htssb(model->_structure_tssb->_root->_transition_tssb);
	}
	c_printf("[*]%s\n", "structure");
	model->_structure_tssb->dump();
	Node* target = model->_structure_tssb->find_node_with_id(12);
	model->delete_node_if_needed(target);
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
		vector<Node*> nodes_on_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_on_htssb);
		for(const auto node_on_htssb: nodes_on_htssb){
			for(int i = 0;i < 1000;i++){
				model->add_customer_to_htssb_node(node_on_htssb);
			}
		}
	}
	for(const auto node: nodes){
		vector<Node*> nodes_on_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_on_htssb);
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
		for(const auto node_on_htssb: nodes_on_htssb){
			node_on_htssb->dump();
			double ratio_v = model->compute_expectation_of_vertical_sbr_ratio(node_on_htssb, true);
			double ratio_h = model->compute_expectation_of_horizontal_sbr_ratio(node_on_htssb, true);
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
		int token_id = Sampler::uniform_int(0, 100);
		model->add_customer_to_hpylm(target, token_id);
	}
	Node* parent = target;
	while(parent){
		c_printf("[*]%d\n", parent->_identifier);
		double pw = model->compute_word_probability_given_node(5, parent);
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
		double pw = copy->compute_word_probability_given_node(5, parent);
		cout << pw << endl;
		assert(parent->_hpylm != NULL);
		parent->_hpylm->dump();
		parent = parent->_parent;
	}
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
		vector<Node*> nodes_on_htssb;
		node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_on_htssb);
		c_printf("[*]%d\n", node->_identifier);
		node->_transition_tssb->dump();
		for(const auto node_on_htssb: nodes_on_htssb){
			node_on_htssb->dump();
			double ratio_v = copy->compute_expectation_of_vertical_sbr_ratio(node_on_htssb, true);
			double ratio_h = copy->compute_expectation_of_horizontal_sbr_ratio(node_on_htssb, true);
			cout << ratio_v << ", " << ratio_h << endl;
		}
	}
}
void _test19(iTHMM* model, TSSB* tssb){
	model->update_stick_length_of_tssb(tssb);
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
		vector<Node*> nodes_on_htssb;
		_test19(model, node->_transition_tssb);
	}
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
		Node* node = copy->sample_node_on_tssb(copy->_bos_tssb);
		copy->add_customer_to_tssb_node(node);
	}

	c_printf("[*]%s\n", "bos");
	copy->_bos_tssb->dump();
}

void test21(){
	string filename = "../alice.txt";
	PyInfiniteTreeHMM* model = new PyInfiniteTreeHMM();
	model->load_textfile(filename);

	model->mark_low_frequency_words_as_unknown(1);
	model->compile();
	model->update_hyperparameters();
	model->_ithmm->_structure_tssb->dump();
}

void test22(){
	iTHMM* model = new iTHMM();
	Node* node = model->sample_node_on_tssb(model->_structure_tssb);
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
}

void test23(){
	iTHMM* model = new iTHMM();
	vector<Node*> nodes = {model->_structure_tssb->_root};
	for(int n = 0;n < 1000;n++){
		Node* node = model->sample_node_on_htssb(nodes[0]->_transition_tssb);
		node->dump();
		assert(node->_structure_tssb_myself);
		nodes.push_back(node->_structure_tssb_myself);
		std::random_shuffle(nodes.begin(), nodes.end());
	}
	int rand_index = Sampler::uniform_int(0, nodes.size());
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
}

int main(){
	iTHMM* model = new iTHMM();
	test21();
	return 0;
}