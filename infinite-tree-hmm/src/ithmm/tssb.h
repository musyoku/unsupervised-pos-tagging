#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <vector>
#include <iostream>

namespace ithmm {
	class Node;
	class TSSB {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, unsigned int version);
	public:
		Node* _root;
		int _owner_id;		// HTSSBの場合は木構造上でこのTSSBを持っているノードのID
		Node* _owner;
		int _num_customers;
		TSSB();
		TSSB(Node* root);
		~TSSB();
		void _delete_children(Node* node);
		void enumerate_nodes_from_left_to_right(std::vector<Node*> &nodes);
		void _enumerate_nodes_from_left_to_right(Node* node, std::vector<Node*> &nodes);
		Node* find_node_by_tracing_horizontal_indices(Node* base);
		Node* find_node_with_id(int identifier);
		Node* _find_node_with_id(int identifier,  Node* node);
		int get_num_nodes();
		int _get_num_children(Node* node);
		int get_max_depth();
		int _get_max_depth(Node* node);
		int get_num_customers();
		int _get_num_customers(Node* node);
		void increment_num_customers();
		void decrement_num_customers();
		void dump();
		void _dump(Node* node);
	};
}