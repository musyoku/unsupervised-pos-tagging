#include "tssb.hpp"
#include "node.hpp"

TSSB::TSSB(double alpha, double gamma, double lambda){
	_root = new Node(NULL);
	_root->_stick_length = 1;
	_owner_id = 0;
	_alpha = alpha;
	_gamma = gamma;
	_lambda = lambda;
}
TSSB::TSSB(Node* root, double alpha, double gamma, double lambda){
	_root = root;
	_owner_id = 0;
	_alpha = alpha;
	_gamma = gamma;
	_lambda = lambda;
}
void TSSB::update_stick_length(){
	double ratio_v = _root->compute_expectation_of_clustering_vertical_sbr_ratio(_gamma);
	double sum_probability = ratio_v;
	_root->_stick_length = 1;
	_root->_children_stick_length = 1.0 - ratio_v;
	_root->_probability = ratio_v;
	_update_stick_length(sum_probability, _root);
	cout << "sum_probability: " << sum_probability << endl;
}
TSSB* TSSB::generate_transition_tssb_belonging_to(int owner_node_id_on_structure){
	Node* root = new Node(NULL, _root->_identifier);
	root->_htssb_owner_id = owner_node_id_on_structure;
	copy_children(_root, root, owner_node_id_on_structure);
	TSSB* target = new TSSB(root, _alpha, _gamma, _lambda);
	target->_owner_id = owner_node_id_on_structure;
	return target;
}
void TSSB::copy_children(Node* source, Node* target, int owner_node_id_on_structure){
	for(const auto source_child: source->_children){
		Node* child = new Node(target, source_child->_identifier);
		child->_htssb_owner_id = owner_node_id_on_structure;
		target->add_child(child);
		copy_children(source_child, child, owner_node_id_on_structure);
	}
}
void TSSB::_update_stick_length(double &sum_probability, Node* node){
	assert(node->_children_stick_length > 0);
	double rest_stick_length = node->_children_stick_length;
	for(int i = 0;i < node->_children.size();i++){
		Node* child = node->_children[i];
		double ratio_h = child->compute_expectation_of_clustering_horizontal_sbr_ratio(_gamma);
		child->_stick_length = rest_stick_length * ratio_h;
		double ratio_v = child->compute_expectation_of_clustering_vertical_sbr_ratio(_gamma);
		child->_probability = child->_stick_length * ratio_v;
		sum_probability += child->_probability;
		rest_stick_length *= 1.0 - ratio_h;
		double alpha = _alpha * pow(_lambda, child->_depth_v);
		child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
		if(child->has_child()){
			_update_stick_length(sum_probability, child);
		}
	}
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
void TSSB::dump(){
	_dump(_root);
}
void TSSB::_dump(Node* node){
	string tab = "";
	for(int i = 0;i < node->_depth_v;i++){
		tab += "	";
	}
	cout << tab;
	int pass_count_v = node->get_vertical_pass_count();
	int stop_count_v = node->get_vertical_stop_count();
	int pass_count_h = node->get_horizontal_pass_count();
	int stop_count_h = node->get_horizontal_stop_count();
	string indices_str = "";
	for(int i = 0;i < node->_depth_v;i++){
		indices_str += node->_horizontal_indices_from_root[i];
		indices_str += ",";
	}
	cout << (boost::format("%d [vp:%d,vs:%d,hp:%d,hs:%d][len:%f,self:%f,ch:%f,p:%f][ow:%d,dv:%d,dh:%d][%s]") % node->_identifier % pass_count_v % stop_count_v % pass_count_h % stop_count_h % node->_stick_length % (node->_stick_length - node->_children_stick_length) % node->_children_stick_length % node->_probability % node->_htssb_owner_id % node->_depth_v % node->_depth_h % indices_str.c_str()).str() << endl;
	for(int i = 0;i < node->_children.size();i++){
		Node* child = node->_children[i];
		_dump(child);
	}
}