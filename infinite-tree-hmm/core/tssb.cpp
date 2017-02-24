#include <boost/format.hpp>
#include "tssb.hpp"
#include "node.hpp"

TSSB::TSSB(){
	_root = NULL;
	_owner_id = 0;
	_owner = NULL;
	_alpha = 0;
	_gamma = 0;
	_lambda = 0;
}
TSSB::TSSB(double alpha, double gamma, double lambda){
	_root = new Node(NULL);
	_root->_stick_length = 1;
	_owner_id = 0;
	_owner = NULL;
	_alpha = alpha;
	_gamma = gamma;
	_lambda = lambda;
}
TSSB::TSSB(Node* root, double alpha, double gamma, double lambda){
	_root = root;
	_owner_id = 0;
	_owner = NULL;
	_alpha = alpha;
	_gamma = gamma;
	_lambda = lambda;
}
void TSSB::enumerate_nodes_from_left_to_right(vector<Node*> &nodes){
	_enumerate_nodes_from_left_to_right(_root, nodes);
}
void TSSB::_enumerate_nodes_from_left_to_right(Node* node, vector<Node*> &nodes){
	nodes.push_back(node);
	for(const auto child: node->_children){
		_enumerate_nodes_from_left_to_right(child, nodes);
	}
}
Node* TSSB::find_node_by_tracing_horizontal_indices(Node* base){
	assert(base != NULL);
	int* indices = base->_horizontal_indices_from_root;
	int depth_v = base->_depth_v;
	Node* iterator = _root;
	for(int i = 0;i < depth_v;i++){
		int index = indices[i];
		assert(index >= 0);
		assert(index < iterator->_children.size());
		iterator = iterator->_children[index];
	}
	return iterator;
}
Node* TSSB::find_node_with_id(int identifier){
	if(_root->_identifier == identifier){
		return _root;
	}
	return _find_node_with_id(identifier, _root);
}
Node* TSSB::_find_node_with_id(int identifier,  Node* node){
	for(const auto child: node->_children){
		if(child->_identifier == identifier){
			return child;
		}
		Node* found = _find_node_with_id(identifier, child);
		if(found != NULL){
			return found;
		}
	}
	return NULL;
}
int TSSB::get_num_nodes(){
	return _get_num_children(_root) + 1;
}
int TSSB::_get_num_children(Node* node){
	int sum = node->_children.size();
	for(const auto child: node->_children){
		sum += _get_num_children(child);
	}
	return sum;
}
int TSSB::get_max_depth(){
	return _get_max_depth(_root);
}
int TSSB::_get_max_depth(Node* node){
	int max_depth = node->_depth_v;
	for(const auto child: node->_children){
		int depth = _get_max_depth(child);
		if(depth > max_depth){
			max_depth = depth;
		}
	}
	return max_depth;
}
void TSSB::increment_num_customers(){
	_num_customers += 1;
}
void TSSB::decrement_num_customers(){
	_num_customers -= 1;
	assert(_num_customers >= 0);
}
void TSSB::dump(){
	_dump(_root);
}
void TSSB::_dump(Node* node){
	string tab = "";
	for(int i = 0;i < node->_depth_v;i++){
		tab += "	";
	}
	cout << tab;
	cout << node->_dump() << endl;
	for(int i = 0;i < node->_children.size();i++){
		Node* child = node->_children[i];
		_dump(child);
	}
}