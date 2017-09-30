#include <boost/format.hpp>
#include "tssb.h"
#include "node.h"

namespace ithmm {
	using std::cout;
	using std::endl;

	TSSB::TSSB(){
		_root = new Node(NULL);
		_root->_stick_length = 1;
		_owner = NULL;
		_num_customers = 0;
		_is_bos = false;
		_is_structure = false;
		_is_htssb = false;
	}
	TSSB::TSSB(Node* root){
		_root = root;
		_owner = NULL;
		_num_customers = 0;
		_is_bos = false;
		_is_structure = false;
		_is_htssb = false;
	}
	TSSB::~TSSB(){
		_delete_children(_root);
		if(_root->is_structure_node()){
			TSSB* transition_tssb = _root->get_transition_tssb();
			if(transition_tssb != NULL){
				delete transition_tssb;
			}
		}
		delete _root;
	}
	int TSSB::get_owner_node_id(){
		assert(_owner != NULL);
		return _owner->_identifier;
	}
	template <class Archive>
	void TSSB::serialize(Archive & ar, unsigned int version)
	{
		ar & _root;
		ar & _owner;
		ar & _num_customers;
	}
	template void TSSB::serialize(boost::archive::binary_iarchive &ar, unsigned int version);
	template void TSSB::serialize(boost::archive::binary_oarchive &ar, unsigned int version);
	void TSSB::_delete_children(Node* node){
		for(auto &child: node->_children){
			if(child->has_child()){
				_delete_children(child);
			}
			if(child->is_structure_node()){
				TSSB* transition_tssb = child->get_transition_tssb();
				if(transition_tssb != NULL){
					delete transition_tssb;
				}
			}
			delete child;
		}
	}
	void TSSB::enumerate_nodes_from_left_to_right(std::vector<Node*> &nodes){
		_enumerate_nodes_from_left_to_right(_root, nodes);
	}
	void TSSB::_enumerate_nodes_from_left_to_right(Node* node, std::vector<Node*> &nodes){
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
	int TSSB::get_num_customers(){
		return _get_num_customers(_root);
	}
	int TSSB::_get_num_customers(Node* node){
		int count = node->_table_v->_num_customers + node->_table_h->_num_customers;
		for(const auto &child: node->_children){
			count += _get_num_customers(child);
		}
		return count;
	}
	void TSSB::increment_num_customers(){
		_num_customers += 1;
	}
	void TSSB::decrement_num_customers(){
		_num_customers -= 1;
		if(_num_customers < 0){
			dump();
			cout << _num_customers << endl;
			cout << get_num_customers() << endl;
		}
		assert(_num_customers >= 0);
	}
	void TSSB::dump(){
		cout << (boost::format("TSSB[ow:$%d,#c:%d,#ct:%d]") % get_owner_node_id() % _num_customers % get_num_customers()).str() << endl;
		_dump(_root);
	}
	void TSSB::_dump(Node* node){
		std::string tab = "";
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
}