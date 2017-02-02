#ifndef _hpylm_
#define _hpylm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>
#include <random>
#include <unordered_map> 
#include <cstdlib>
#include <cassert>
#include "cprintf.h"
#include "node.h"
#include "sampler.h"

class CRP{
private:
	friend class boost::serialization::access;
	template <class Archive>
	// モデルの保存
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version); // No use
		archive & _root;
		archive & _g0;
		archive & _d_m;
		archive & _theta_m;
		archive & _a_m;
		archive & _b_m;
		archive & _alpha_m;
		archive & _beta_m;
	}
public:
	Node* _root;				// 文脈木のルートノード
	double _alpha;
	double _beta;
	double _gamma;
	CRP(){
		_root = new Node();
		_root->_depth = 0;		// ルートは深さ0
		_alpha = 1;
		_beta = 1;
		_gamma = 1;
	}
	bool add_customer(int context_token_id, int token_id){
		Node* node = find_context_node(context_token_id, true);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "客を追加できません. ノードが見つかりません.");
			exit(1);
		}
		int added_to_table_k;
		node->add_customer(token_id, _g0, _alpha, _beta, _gamma, true, added_to_table_k);
		return true;
	}
	bool remove_customer(int context_token_id, int token_id){
		Node* node = find_context_node(context_token_id, false);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. ノードが見つかりません.");
			exit(1);
		}
		int removed_from_table_k;
		node->remove_customer(token_id, _alpha, _beta, _gamma, true, removed_from_table_k);
		// 客が一人もいなくなったらノードを削除する
		if(node->need_to_remove_from_parent()){
			node->remove_from_parent();
		}
		return true;
	}
	Node* find_context_node(int context_token_id, bool generate_node_if_needed = false, bool return_middle_node = false){
		Node* child = _root->find_child_node(context_token_id, generate_node_if_needed);
		if(return_middle_node && child == NULL){
			return _root;
		}
		return child;
	}
	double compute_Pw_h(int context_token_id, int token_id){
		Node* node = find_context_node(context_token_id, false, true);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. node == NULL");
			exit(1);
		}
		return node->compute_Pw(token_id, _g0, _concent_m);
	}
	bool save(string filename = "hpylm.model"){
		std::ofstream ofs(filename);
		boost::archive::binary_oarchive oarchive(ofs);
		oarchive << static_cast<const HPYLM&>(*this);
		ofs.close();
		return true;
	}
	bool load(string filename = "hpylm.model"){
		std::ifstream ifs(filename);
		if(ifs.good() == false){
			return false;
		}
		boost::archive::binary_iarchive iarchive(ifs);
		iarchive >> *this;
		ifs.close();
		return true;
	}
};

#endif