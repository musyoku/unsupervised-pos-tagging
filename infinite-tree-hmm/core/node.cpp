#include <boost/format.hpp>
#include <cassert>
#include "cprintf.h"
#include "sampler.h"
#include "util.h"
#include "hpylm.hpp"
#include "tssb.hpp"
#include "node.hpp"

Node::Node(Node* parent){
	_identifier = _auto_increment;
	_auto_increment++;
	_parent = parent;
	init();
}
Node::Node(Node* parent, int identifier){
	_identifier = identifier;
	_parent = parent;
	init();
}
void Node::init(){
	_depth_v = (_parent != NULL) ? _parent->_depth_v + 1 : 0;
	_depth_h = (_parent != NULL) ? _parent->_children.size() : 0;
	_htssb_owner_id = (_parent != NULL) ? _parent->_htssb_owner_id : 0;
	_stick_length = -1;
	_children_stick_length = -1;
	_pass_count_v = 0;
	_stop_count_v = 0;
	_pass_count_h = 0;
	_stop_count_h = 0;
	_probability = -1;
	_transition_tssb = NULL;
	_transition_tssb_myself = NULL;
	_parent_transition_tssb_myself = NULL;
	_htssb_owner = (_parent != NULL) ? _parent->_htssb_owner : 0;
	_stop_probability_v_over_parent = new double[_depth_v + 1];
	_stop_ratio_v_over_parent = new double[_depth_v + 1];
	_nodes_from_root_to_myself = new Node*[_depth_v + 1];
	_stop_probability_h_over_parent = new double[_depth_h + 1];
	_stop_ratio_h_over_parent = new double[_depth_h + 1];
	_horizontal_indices_from_root = new int[_depth_v];
	_table_v = new Table();
	_table_h = new Table();
	_hpylm = NULL;
	if(_htssb_owner_id == 0){	// HPYLMは木構造上のノードにだけあればよい
		_hpylm = new HPYLM(this);
	}
	set_horizontal_indices();
	set_pointers_from_root_to_myself();
}
Node::~Node(){
	delete _table_v;
	delete _table_h;
	delete[] _stop_probability_v_over_parent;
	delete[] _stop_ratio_v_over_parent;
	delete[] _nodes_from_root_to_myself;
	delete[] _horizontal_indices_from_root;
}
Node* Node::generate_child(){
	Node* child = new Node(this);
	add_child(child);
	return child;
}
void Node::copy_transition_tssb_from_structure(TSSB* structure){
	_transition_tssb = structure->generate_transition_tssb_belonging_to(this);
}
void Node::set_horizontal_indices(){
	Node* iterator = this;
	for(int i = 0;i < _depth_v;i++){
		_horizontal_indices_from_root[_depth_v - i - 1] = iterator->_depth_h;
		iterator = iterator->_parent;
	}
}
void Node::set_pointers_from_root_to_myself(){
	Node* iterator = this;
	_nodes_from_root_to_myself[_depth_v] = iterator;
	for(int n = 0;n < _depth_v;n++){
		iterator = iterator->_parent;
		assert(iterator != NULL);
		_nodes_from_root_to_myself[_depth_v - n - 1] = iterator;
	}
}
void Node::add_child(Node* node){
	assert(node != NULL);
	_children.push_back(node);
}
Node* Node::find_same_node_on_transition_tssb(){
	return _transition_tssb->find_node_by_tracing_horizontal_indices(this);
}
int Node::get_vertical_stop_count(){
	return _stop_count_v;
}
int Node::get_vertical_pass_count(){
	return _pass_count_v;
}
int Node::get_horizontal_stop_count(){
	return _stop_count_h;
}
int Node::get_horizontal_pass_count(){
	return _pass_count_h;
}
Table* Node::get_vertical_table(){
	return _table_v;
}
Table* Node::get_horizontal_table(){
	return _table_h;
}
bool Node::has_child(){
	return _children.size() != 0;
}
// 遷移確率TSSBに客を追加
void Node::add_customer_to_vertical_crp(double concentration, double g0, bool &new_table_generated){
	Table* table = get_vertical_table();
	Node* parent = _parent;
	assert(table != NULL);
	int total_count = _pass_count_v + _stop_count_v;
	table->add_customer(concentration, g0, total_count, new_table_generated);
	// 停止回数・通過回数を更新
	increment_vertical_stop_count();
	while(parent){
		parent->increment_vertical_pass_count();
		parent = parent->_parent;
	}
}
void Node::increment_vertical_stop_count(){
	_stop_count_v += 1;
}
void Node::decrement_vertical_stop_count(){
	_stop_count_v -= 1;
	assert(_stop_count_v >= 0);
}
void Node::increment_vertical_pass_count(){
	_pass_count_v += 1;
}
void Node::decrement_vertical_pass_count(){
	_pass_count_v -= 1;
	assert(_pass_count_v >= 0);
}
void Node::add_customer_to_horizontal_crp(double concentration, double g0, bool &new_table_generated){
	Table* table = get_horizontal_table();
	assert(table != NULL);
	int total_count = _pass_count_h + _stop_count_h;
	table->add_customer(concentration, g0, total_count, new_table_generated);
	// 停止回数・通過回数を更新
	Node* stopped_child = this;
	Node* parent = _parent;
	while(parent){
		for(int i = 0;i < parent->_children.size();i++){
			Node* child = parent->_children[i];
			if(child == stopped_child){
				child->increment_horizontal_stop_count();
				break;
			}
			child->increment_horizontal_pass_count();
		}
		stopped_child = parent;
		parent = parent->_parent;
	}
	stopped_child->increment_horizontal_stop_count();
}
void Node::increment_horizontal_stop_count(){
	_stop_count_h += 1;
}
void Node::decrement_horizontal_stop_count(){
	_stop_count_h -= 1;
	assert(_stop_count_h >= 0);
}
void Node::increment_horizontal_pass_count(){
	_pass_count_h += 1;
}
void Node::decrement_horizontal_pass_count(){
	_pass_count_h -= 1;
	assert(_pass_count_h >= 0);
}
// 客を除去
void Node::remove_customer_from_vertical_crp(bool &empty_table_deleted){
	// cout << "remove_customer_from_vertical_crp: " << tssb_identifier << ", " << node->_identifier << endl;
	Table* table = get_vertical_table();
	assert(table != NULL);
	table->remove_customer(empty_table_deleted);
	// 停止回数・通過回数を更新
	decrement_vertical_stop_count();
	Node* parent = _parent;
	while(parent){
		parent->decrement_vertical_pass_count();
		parent = parent->_parent;
	}
}
void Node::remove_customer_from_horizontal_crp(bool &empty_table_deleted){
	// cout << "remove_customer_from_horizontal_crp: " << _identifier << "," << node->_identifier << endl;
	Table* table = get_horizontal_table();
	assert(table != NULL);
	table->remove_customer(empty_table_deleted);
	// 停止回数・通過回数を更新
	Node* stopped_child = this;
	Node* parent = _parent;
	while(parent){
		bool found = false;
		for(int i = parent->_children.size() - 1;i >= 0;i--){	// 逆向きに辿らないと通過ノードが先に消えてしまう
			Node* child = parent->_children[i];
			// cout << "foreach: " << child->_identifier << endl;
			// cout << "foreach: " << child->_identifier << endl;

			if(child == stopped_child){
				found = true;
				child->decrement_horizontal_stop_count();
				continue;
			}
			if(found){
				child->decrement_horizontal_pass_count();
			}
		}
		stopped_child = parent;
		parent = parent->_parent;
	}
	// ルートノードのカウントを減らす
	stopped_child->decrement_horizontal_stop_count();
}
bool Node::delete_node_if_needed(){
	int pass_count_v = get_vertical_pass_count();
	int stop_count_v = get_vertical_stop_count();
	int pass_count_h = get_horizontal_pass_count();
	int stop_count_h = get_horizontal_stop_count();
	if(pass_count_v + stop_count_v + pass_count_h + stop_count_h == 0 && _parent != NULL){
		// cout << pass_count_v << "," << stop_count_v << "," << pass_count_h << "," << stop_count_h << endl;
		// cout << "requesting parent " << _parent->_identifier << ", me = " << _identifier << endl;
		_parent->delete_child_node(_identifier);
		return true;
	}
	return false;
}
Node* Node::delete_child_node(int node_id){
	Node* return_node = NULL;
	for(int i = 0;i < _children.size();i++){
		Node* &target = _children[i];
		if(target->_identifier == node_id){
			assert(target->get_vertical_pass_count() == 0);
			assert(target->get_vertical_stop_count() == 0);
			assert(target->get_horizontal_pass_count() == 0);
			assert(target->get_horizontal_stop_count() == 0);
			return_node = target;
			continue;
		}
		if(return_node != NULL){
			target->_depth_h -= 1;
			target->_horizontal_indices_from_root[target->_depth_v - 1] -= 1;
		}
	}
	if(return_node != NULL){
		_children.erase(_children.begin() + return_node->_depth_h);
	}
	return return_node;
}
void Node::dump(){
	string indices_str = "";
	for(int i = 0;i < _depth_v;i++){
		indices_str += std::to_string(_horizontal_indices_from_root[i]);
		indices_str += ",";
	}
	cout << (boost::format("%d [vp:%d,vs:%d,hp:%d,hs:%d][len:%f,self:%f,ch:%f,p:%f][ow:%d,dv:%d,dh:%d][%s]") % _identifier % _pass_count_v % _stop_count_v % _pass_count_h % _stop_count_h % _stick_length % (_stick_length - _children_stick_length) % _children_stick_length % _probability % _htssb_owner_id % _depth_v % _depth_h % indices_str.c_str()).str() << endl;
}

int Node::_auto_increment = 1;