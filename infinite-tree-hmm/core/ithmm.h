#ifndef _ithmm_
#define _ithmm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <vector>
#include <fstream>
#include "tssb.hpp"
#include "node.hpp"
#include "hpylm.hpp"
#include "sampler.h"
#include "cprintf.h"
#include "util.h"
#include "hyperparameters.h"

// 10以下ならなんでもいい
#define TSSB_STRUCTURE_ID 8
#define TSSB_BOS_ID 7

typedef struct Word {
	id id;
	Node* state;
} Word;

class iTHMM{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _structure_tssb;
		archive & _bos_tssb;
		archive & _alpha;
		archive & _gamma;
		archive & _lambda;
		archive & _tau0;
		archive & _tau1;
		archive & _word_g0;
		archive & _max_depth;
		archive & _hpylm_d_m;
		archive & _hpylm_theta_m;
		archive & _hpylm_a_m;
		archive & _hpylm_b_m;
		archive & _hpylm_alpha_m;
		archive & _hpylm_beta_m;
	}
public:
	TSSB* _structure_tssb;	// 木構造を表すためだけのTSSB。HTSSBは全てこのノードを基準に成形する
	TSSB* _bos_tssb;		// <bos>からの遷移を表すTSSB
	double _alpha;
	double _gamma;
	double _lambda;
	double _tau0;
	double _tau1;
	double _word_g0;
	int _max_depth;
	vector<double> _hpylm_d_m;		// HPYLMのハイパーパラメータ（ディスカウント係数）
	vector<double> _hpylm_theta_m;	// HPYLMのハイパーパラメータ（集中度）
	vector<double> _hpylm_a_m;		// ベータ分布のパラメータ	dの推定用
	vector<double> _hpylm_b_m;		// ベータ分布のパラメータ	dの推定用
	vector<double> _hpylm_alpha_m;	// ガンマ分布のパラメータ	θの推定用
	vector<double> _hpylm_beta_m;	// ガンマ分布のパラメータ	θの推定用
	iTHMM(){
		_alpha = iTHMM_ALPHA;
		_gamma = iTHMM_GAMMA;
		_lambda = iTHMM_LAMBDA;
		_tau0 = iTHMM_TAU_0;
		_tau1 = iTHMM_TAU_1;
		_max_depth = 0;
		_word_g0 = -1;

		_structure_tssb = new TSSB(_alpha, _gamma, _lambda);
		_structure_tssb->_root->_owner_id_on_structure = TSSB_STRUCTURE_ID;
		_structure_tssb->_owner_id = TSSB_STRUCTURE_ID;
		Node* root_on_structure = _structure_tssb->_root;
		root_on_structure->init_hpylm();
		Node* root_on_htssb = new Node(NULL, root_on_structure->_identifier);
		root_on_htssb->_owner_id_on_structure = root_on_structure->_identifier;
		root_on_htssb->_owner_on_structure = root_on_structure;
		root_on_htssb->_structure_tssb_myself = root_on_structure;
		root_on_structure->_transition_tssb = new TSSB(root_on_htssb, _alpha, _gamma, _lambda);
		root_on_structure->_transition_tssb->_owner_id = root_on_structure->_identifier;
		root_on_structure->_transition_tssb_myself = root_on_htssb;

		Node* root_on_bos = new Node(NULL, root_on_structure->_identifier);
		root_on_bos->_owner_id_on_structure = TSSB_BOS_ID;		// そもそも木構造上に所有者がいないが気にしない
		root_on_bos->_structure_tssb_myself = root_on_structure;
		_bos_tssb = new TSSB(root_on_bos, _alpha, _gamma, _lambda);
		_bos_tssb->_owner_id = TSSB_BOS_ID;

		_hpylm_d_m.push_back(HPYLM_D);
		_hpylm_theta_m.push_back(HPYLM_THETA);
		_hpylm_a_m.push_back(HPYLM_A);
		_hpylm_b_m.push_back(HPYLM_B);
		_hpylm_alpha_m.push_back(HPYLM_ALPHA);
		_hpylm_beta_m.push_back(HPYLM_BETA);
	}
	void initialize_data(vector<vector<Word*>> &dataset){
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			if(line.size() == 0){
				continue;
			}
			// 状態路ランダムに設定
			for(int i = 0;i < line.size();i++){
				Word* word = line[i];
				Node* state = sample_node_on_tssb(_structure_tssb);
				assert(state != NULL);
				word->state = state;
			}
			Node* prev_state = NULL;
			Node* next_state = line.size() == 1 ? NULL : line[1]->state;
			for(int i = 0;i < line.size();i++){
				Word* word = line[i];
				Node* state = word->state;
				add_parameters(prev_state, state, next_state, word->id);
				prev_state = state;
				next_state = i < line.size() - 2 ? line[i + 2]->state : NULL;
			}
		}
	}
	// デバッグ用
	// これを呼んで全パラメータが消えなかったらバグっている
	void remove_all_data(vector<vector<Word*>> &dataset){
		for(int data_index = 0;data_index < dataset.size();data_index++){
			vector<Word*> &line = dataset[data_index];
			if(line.size() == 0){
				continue;
			}
			Node* prev_state = NULL;
			Node* next_state = line.size() == 1 ? NULL : line[1]->state;
			for(int i = 0;i < line.size();i++){
				Word* word = line[i];
				Node* state = word->state;
				remove_parameters(prev_state, state, next_state, word->id);
				prev_state = state;
				next_state = i < line.size() - 2 ? line[i + 2]->state : NULL;
			}
		}
	}
	void set_word_g0(double g0){
		_word_g0 = g0;
	}
	bool is_node_on_bos_tssb(Node* node){
		assert(node != NULL);
		return node->_owner_id_on_structure == TSSB_BOS_ID;
	}
	bool is_node_on_structure_tssb(Node* node){
		assert(node != NULL);
		return node->_owner_id_on_structure == TSSB_STRUCTURE_ID;
	}
	bool is_node_on_htssb(Node* node){
		assert(node != NULL);
		return is_node_on_bos_tssb(node) == false && is_node_on_structure_tssb(node) == false;
	}
	bool is_tssb_bos(TSSB* tssb){
		assert(tssb != NULL);
		return tssb->_owner_id == TSSB_BOS_ID;
	}
	bool is_tssb_structure(TSSB* tssb){
		assert(tssb != NULL);
		return tssb->_owner_id == TSSB_STRUCTURE_ID;
	}
	bool is_tssb_htssb(TSSB* tssb){
		return is_tssb_bos(tssb) == false && is_tssb_structure(tssb) == false;
	}
	bool is_node_to_the_left_of_node(Node* left, Node* right){
		assert(left->_identifier != right->_identifier);
		int limit = std::min(left->_depth_v, right->_depth_v);
		for(int i = 0;i < limit;i++){
			if(left->_horizontal_indices_from_root[i] != right->_horizontal_indices_from_root[i]){
				return left->_horizontal_indices_from_root[i] < right->_horizontal_indices_from_root[i];
			}
		}
		return left->_depth_v < right->_depth_v;
	}
	// 木構造で子ノードを生成した際に全てのHTSSBの同じ位置に子ノードを生成する
	Node* generate_and_add_new_child_to(Node* parent){
		assert(parent != NULL);
		// まず木構造上で子ノードを作る
		Node* generated_child_on_structure = NULL;
		if(is_node_on_structure_tssb(parent)){	// parentが木構造上のノードの場合
			generated_child_on_structure = parent->generate_child();
		}else{	// parentが別のノードの遷移確率用TSSB上のノードだった場合
			Node* parent_on_structure = _structure_tssb->find_node_by_tracing_horizontal_indices(parent);
			generated_child_on_structure = parent_on_structure->generate_child();
		}
		assert(generated_child_on_structure != NULL);
		// HTSSBをセット
		generated_child_on_structure->_transition_tssb = generate_transition_tssb_belonging_to(generated_child_on_structure);
		Node* myself_on_htssb = generated_child_on_structure->find_same_node_on_transition_tssb();
		assert(myself_on_htssb != NULL);
		generated_child_on_structure->_transition_tssb_myself = myself_on_htssb;
		// HPYLM
		generated_child_on_structure->init_hpylm();

		Node* return_child = generated_child_on_structure;	// 実際に返すノード
		// <bos>TSSB上で子ノードを作成
		Node* generated_child_on_bos = _generate_and_add_new_child_to_bos_tssb(generated_child_on_structure);
		if(is_node_on_bos_tssb(parent)){
			return_child = generated_child_on_bos;
		}
		// 木構造上の全ノードのHTSSBにノードを追加
		_generate_and_add_new_child_to_all_htssb(_structure_tssb->_root, parent, generated_child_on_structure, return_child);
		// ポインタを張る
		Node* generated_child_on_htssb = generated_child_on_structure->_transition_tssb_myself;
		assert(generated_child_on_htssb != NULL);
		// generated_child_on_htssb->_structure_tssb_myself = generated_child_on_structure;
		//// 木構造上の親ノードのHTSSBの自分と同じ位置のノードへのポインタ
		Node* iterator_on_structure = generated_child_on_structure;
		Node* parent_on_structure = iterator_on_structure->_parent;
		Node* iterator_on_htssb = generated_child_on_htssb;
		Node* iterator_on_parent_htssb = NULL;
		while(parent_on_structure != NULL){
			assert(iterator_on_structure->_transition_tssb_myself != NULL);
			// 木構造上での親ノードが持つHTSSBにある対応するノードを取る
			iterator_on_parent_htssb = parent_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(generated_child_on_structure);
			assert(iterator_on_parent_htssb != NULL);
			// ポインタを張る
			iterator_on_htssb->_parent_transition_tssb_myself = iterator_on_parent_htssb;
			iterator_on_htssb->_structure_tssb_myself = generated_child_on_structure;
			assert(iterator_on_htssb->_structure_tssb_myself->_identifier == generated_child_on_structure->_identifier);
			// 木構造上で次の親ノードへ
			iterator_on_structure = parent_on_structure;
			parent_on_structure = iterator_on_structure->_parent;
			iterator_on_htssb = iterator_on_parent_htssb;
		}
		// HPYLM用のハイパーパラメータを追加
		if(return_child->_depth_v > _max_depth){
			_max_depth = return_child->_depth_v;		
			while(_max_depth >= _hpylm_d_m.size()){
				_hpylm_d_m.push_back(HPYLM_D);
			}			
			while(_max_depth >= _hpylm_theta_m.size()){
				_hpylm_theta_m.push_back(HPYLM_THETA);
			}		
			while(_max_depth >= _hpylm_a_m.size()){
				_hpylm_a_m.push_back(HPYLM_A);
			}			
			while(_max_depth >= _hpylm_b_m.size()){
				_hpylm_b_m.push_back(HPYLM_B);
			}		
			while(_max_depth >= _hpylm_alpha_m.size()){
				_hpylm_alpha_m.push_back(HPYLM_ALPHA);
			}			
			while(_max_depth >= _hpylm_beta_m.size()){
				_hpylm_beta_m.push_back(HPYLM_BETA);
			}	
		}
		return return_child;
	}
	void _generate_and_add_new_child_to_all_htssb(Node* iterator_on_structure, Node* parent, Node* generated_child_on_structure, Node* &return_child){
		// iteratorとgenerated_childが同一の場合はすでに追加されているのでスキップ
		if(iterator_on_structure->_identifier != generated_child_on_structure->_identifier){
			assert(iterator_on_structure->_transition_tssb != NULL);
			int owner_id_of_htssb_parent_belongs = parent->_owner_id_on_structure;
			int child_id_to_generate = generated_child_on_structure->_identifier;
			// 遷移確率用TSSBの同じ位置に子ノードを挿入
			Node* parent_on_htssb = iterator_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(parent);
			assert(parent_on_htssb != NULL);
			assert(parent_on_htssb->_identifier == generated_child_on_structure->_parent->_identifier);
			Node* child_on_htssb = new Node(parent_on_htssb, child_id_to_generate);
			child_on_htssb->_structure_tssb_myself = generated_child_on_structure;
			if(child_on_htssb->_owner_id_on_structure == owner_id_of_htssb_parent_belongs){	// 親と同じTSSB上の子ノードを返す
				return_child = child_on_htssb;
			}
			parent_on_htssb->add_child(child_on_htssb);
		}
		for(const auto &child: iterator_on_structure->_children){
			_generate_and_add_new_child_to_all_htssb(child, parent, generated_child_on_structure, return_child);
		}
	}
	Node* _generate_and_add_new_child_to_bos_tssb(Node* generated_child_on_structure){
		// 木構造上での親ノードが<bos>TSSBのどのノードに対応するかを調べる
		Node* parent = _bos_tssb->find_node_by_tracing_horizontal_indices(generated_child_on_structure->_parent);
		assert(parent != NULL);
		Node* child = new Node(parent, generated_child_on_structure->_identifier);
		child->_owner_id_on_structure = TSSB_BOS_ID;
		parent->add_child(child);
		// ポインタを張る
		generated_child_on_structure->_bos_tssb_myself = child;
		child->_structure_tssb_myself = generated_child_on_structure;
		return child;
	}
	// 木構造上のノードにHTSSBを追加
	TSSB* generate_transition_tssb_belonging_to(Node* owner_on_structure){
		assert(is_node_on_structure_tssb(owner_on_structure));
		Node* root_on_structure = _structure_tssb->_root;
		Node* root_on_htssb = new Node(NULL, root_on_structure->_identifier);
		root_on_htssb->_owner_id_on_structure = owner_on_structure->_identifier;
		root_on_htssb->_owner_on_structure = owner_on_structure;
		root_on_htssb->_parent_transition_tssb_myself = NULL;
		if(owner_on_structure->_parent != NULL){
			root_on_htssb->_parent_transition_tssb_myself = owner_on_structure->_parent->_transition_tssb->_root;
		}
		root_on_htssb->_structure_tssb_myself = root_on_structure;
		copy_children_on_structure_to_transition_tssb(root_on_structure, root_on_htssb, owner_on_structure);
		TSSB* target = new TSSB(root_on_htssb, _alpha, _gamma, _lambda);
		target->_owner_id = owner_on_structure->_identifier;
		target->_owner = owner_on_structure;
		return target;
	}
	// 生成したHTSSBを木構造と同一の形状にするために子ノードを生成・追加
	void copy_children_on_structure_to_transition_tssb(Node* source_on_structure, Node* target_on_htssb, Node* owner_on_structure){
		for(const auto source_child_on_structure: source_on_structure->_children){
			Node* child = new Node(target_on_htssb, source_child_on_structure->_identifier);
			child->_owner_id_on_structure = owner_on_structure->_identifier;
			child->_owner_on_structure = owner_on_structure;
			child->_structure_tssb_myself = source_child_on_structure;
			// child->_owner_id_on_structure = owner_on_structure;
			target_on_htssb->add_child(child);
			copy_children_on_structure_to_transition_tssb(source_child_on_structure, child, owner_on_structure);
		}
	}
	// コインを投げる操作を繰り返して到達したノードを返す
	Node* sample_node_on_tssb(TSSB* tssb){
		assert(is_tssb_structure(tssb));
		Node* node = _sample_node_on_tssb_by_iterating_node(tssb->_root, false);
		return node;
	}
	// HTSSB上でノードをサンプリング
	Node* sample_node_on_htssb(TSSB* tssb){
		assert(is_tssb_htssb(tssb));
		Node* node = _sample_node_on_tssb_by_iterating_node(tssb->_root, true);
		assert(node->_owner_id_on_structure == tssb->_owner_id);
		return node;
	}
	// 止まるノードを決定する
	// htssb_modeがtrueの場合、停止確率は親のHTSSBから生成する
	// htssb_modeがfalseの場合は普通のTSSBによるクラスタリング
	Node* _sample_node_on_tssb_by_iterating_node(Node* iterator, bool htssb_mode){
		assert(iterator != NULL);
		double head = compute_expectation_of_vertical_sbr_ratio(iterator, htssb_mode);
		iterator->_children_stick_length = iterator->_stick_length * (1 - head);
		double bernoulli = Sampler::uniform(0, 1);
		if(bernoulli <= head){			// 表が出たらこのノードに降りる
			return iterator;
		}
		// 子ノードがある場合
		for(int i = 0;i < iterator->_children.size();i++){
			Node* child = iterator->_children[i];
			assert(child != NULL);
			double head = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _sample_node_on_tssb_by_iterating_node(child, htssb_mode);
			}
		}
		// ない場合生成しながらコインを投げる
		while(true){
			Node* child = generate_and_add_new_child_to(iterator);
			double head = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double bernoulli = Sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _sample_node_on_tssb_by_iterating_node(child, htssb_mode);
			}
		}
	}
	// [0, 1)の一様分布からノードをサンプリング
	// HTSSBの場合（木構造に対してはそもそも行わない）
	Node* retrospective_sampling(double uniform, TSSB* tssb, double total_stick_length, bool htssb_mode){
		assert(uniform < total_stick_length);
		assert( (htssb_mode && is_tssb_htssb(tssb)) || (htssb_mode == false && is_tssb_htssb(tssb) == false) );
		Node* root = tssb->_root;
		double ratio_v = compute_expectation_of_vertical_sbr_ratio(root, htssb_mode);
		double sum_probability = total_stick_length * ratio_v;
		root->_stick_length = total_stick_length;
		root->_children_stick_length = total_stick_length * (1.0 - ratio_v);
		Node* node =  _retrospective_sampling_by_iterating_node(uniform, sum_probability, root, htssb_mode);
		assert(node != NULL);
		return node;
	}
	Node* _retrospective_sampling_by_iterating_node(double uniform, double &sum_probability, Node* iterator, bool htssb_mode){
		if(uniform <= sum_probability){
			return iterator;
		}
		// 棒の長さとノードの確率の関係に気をつける
		// [<------------------- 棒の長さ -------------------]
		// [<--親ノードの確率--><---子ノードに割り当てる長さ --->]
		//					  [<---子1の棒---><---子2の棒--->]
		assert(iterator->_children_stick_length > 0);
		double rest_stick_length = iterator->_children_stick_length;	// 子ノードに割り当てる棒の長さの総和
		assert(rest_stick_length > 0);
		double sum_stick_length_over_children = 0;			// 子ノードを走査する時の走査済みの棒の長さ
		// Node* last_node = NULL;
		for(int i = 0;i < iterator->_children.size();i++){
			Node* child = iterator->_children[i];
			double alpha = _alpha * pow(_lambda, child->_depth_v);
			double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
			child->_stick_length = rest_stick_length * ratio_h;
			child->_probability = child->_stick_length * ratio_v;
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
			if(uniform <= sum_probability + sum_stick_length_over_children + child->_stick_length){
				// child->_stick_lengthだけだとこのノードの棒の長さ（つまりこのノード+子ノードに割り当てる棒）なので
				// ratio_vも掛けてこのノードで止まる確率にする必要がある
				sum_probability += sum_stick_length_over_children + child->_stick_length * ratio_v;
				if(uniform <= sum_probability){
					return child;
				}
				// if(child->has_child()){
					return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);
				// }
				// 子ノード領域に当たった場合、uniformを超えるまで棒を折り続ける
				// Node* _child = generate_and_add_new_child_to(child);
				// double ratio_h = compute_expectation_of_horizontal_sbr_ratio(_child, htssb_mode);
				// _child->_stick_length = child->_children_stick_length * ratio_h;
				// double alpha = _alpha * pow(_lambda, _child->_depth_v);
				// double ratio_v = compute_expectation_of_vertical_sbr_ratio(_child, htssb_mode);
				// _child->_probability = _child->_stick_length * ratio_v;
				// _child->_children_stick_length = _child->_stick_length * (1.0 - ratio_v);
				// assert(child->has_child());
				// return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);
			}
			sum_stick_length_over_children += child->_stick_length;
			rest_stick_length *= 1.0 - ratio_h;
			// last_node = child;
		}

		// 全ての存在する子ノードを通り過ぎた場合は生成
		while(true){
			Node* child = generate_and_add_new_child_to(iterator);
			double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
			child->_stick_length = rest_stick_length * ratio_h;
			child->_probability = child->_stick_length * ratio_v;
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
			if(uniform <= sum_probability + sum_stick_length_over_children + child->_stick_length){
				sum_probability += sum_stick_length_over_children + child->_probability;
				if(uniform <= sum_probability){
					return child;
				}
				return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);
			}
			sum_stick_length_over_children += child->_stick_length;
			rest_stick_length *= 1.0 - ratio_h;
		}
		
		// double alpha = _alpha * pow(_lambda, _child->_depth_v);
		// double ratio_v = compute_expectation_of_vertical_sbr_ratio(_child, htssb_mode);
		// _child->_probability = _child->_stick_length * ratio_v;
		// _child->_children_stick_length = _child->_stick_length * (1.0 - ratio_v);
		// assert(child->has_child());
		// return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);

		// if(is_node_on_htssb(iterator)){
		// 	Node* owner = _structure_tssb->find_node_with_id(iterator->_owner_id_on_structure);
		// 	iterator->dump();
		// 	assert(owner != NULL);
		// 	double p = owner->compute_transition_probability_to_eos(_tau0, _tau1);
		// 	update_stick_length_of_tssb(owner->_transition_tssb, 1 - p);
		// 	owner->_transition_tssb->dump();
		// 	last_node->dump();
		// 	cout << uniform << ", ";
		// 	cout << sum_probability << endl;
		// }
		// // 見つからなかったら一番右端のノードを返す
		// if(last_node != NULL){
		// 	if(last_node->has_child()){
		// 		return _retrospective_sampling_by_iterating_node(uniform, sum_probability, last_node, htssb_mode);
		// 	}
		// 	return last_node;
		// }
		// return NULL;
	}
	void perform_gibbs_sampling_line(vector<Word*> &line){
		assert(line.size() > 0);
		Node* prev_state = NULL;
		Node* next_state = line.size() == 1 ? NULL : line[1]->state;
		for(int i = 0;i < line.size();i++){
			Word* word = line[i];
			Node* state = word->state;
			remove_parameters(prev_state, state, next_state, word->id);
			state = draw_state(prev_state, state, next_state, word->id);
			add_parameters(prev_state, state, next_state, word->id);
			prev_state = state;
			next_state = i < line.size() - 2 ? line[i + 2]->state : NULL;
			delete_invalid_children_of_node_on_structure(_structure_tssb->_root);
			word->state = state;
		}
	}
	void add_parameters(Node* prev_state_on_structure, Node* state_on_structure, Node* next_state_on_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_on_structure == NULL){
			assert(state_on_structure != NULL);
			assert(state_on_structure->_transition_tssb != NULL);
			assert(next_state_on_structure != NULL);
			assert(is_node_on_structure_tssb(state_on_structure));
			assert(is_node_on_structure_tssb(next_state_on_structure));
			Node* state_on_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_bos);
			add_customer_to_tssb_node(state_on_bos);
			add_customer_to_tssb_node(state_on_structure);			// 参照カウント用
			Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
			assert(next_state_on_htssb != NULL);
			add_customer_to_htssb_node(next_state_on_htssb);
			add_customer_to_tssb_node(next_state_on_structure);		// 参照カウント用
			add_customer_to_hpylm(state_on_structure, word_id);
			state_on_structure->increment_transition_count_to_other();
			return;
		}
		// <eos>への遷移を含む場合
		if(next_state_on_structure == NULL){
			assert(prev_state_on_structure != NULL);
			assert(prev_state_on_structure->_transition_tssb != NULL);
			assert(state_on_structure != NULL);
			assert(is_node_on_structure_tssb(prev_state_on_structure));
			assert(is_node_on_structure_tssb(state_on_structure));
			Node* state_on_htssb = prev_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_htssb != NULL);
			add_customer_to_htssb_node(state_on_htssb);
			add_customer_to_tssb_node(state_on_structure);		// 参照カウント用
			add_customer_to_hpylm(state_on_structure, word_id);
			prev_state_on_structure->increment_transition_count_to_other();
			state_on_structure->increment_transition_count_to_eos();
			return;
		}
		// <bos>と<eos>両方を含む場合
		if(prev_state_on_structure == NULL && next_state_on_structure == NULL){
			assert(state_on_structure != NULL);
			assert(is_node_on_structure_tssb(state_on_structure));
			Node* state_on_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_bos);
			add_customer_to_tssb_node(state_on_bos);
			add_customer_to_tssb_node(state_on_structure);			// 参照カウント用
			add_customer_to_hpylm(state_on_structure, word_id);
			state_on_structure->increment_transition_count_to_eos();
			return;
		}
		assert(prev_state_on_structure != NULL);
		assert(prev_state_on_structure->_transition_tssb != NULL);
		assert(state_on_structure != NULL);
		assert(state_on_structure->_transition_tssb != NULL);
		assert(next_state_on_structure != NULL);
		assert(is_node_on_structure_tssb(prev_state_on_structure));
		assert(is_node_on_structure_tssb(state_on_structure));
		assert(is_node_on_structure_tssb(next_state_on_structure));

		Node* state_on_htssb = prev_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
		assert(state_on_htssb != NULL);
		add_customer_to_htssb_node(state_on_htssb);
		add_customer_to_tssb_node(state_on_structure);			// 参照カウント用

		Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
		assert(next_state_on_htssb != NULL);
		add_customer_to_htssb_node(next_state_on_htssb);
		add_customer_to_tssb_node(next_state_on_structure);		// 参照カウント用

		add_customer_to_hpylm(state_on_structure, word_id);
		prev_state_on_structure->increment_transition_count_to_other();
		state_on_structure->increment_transition_count_to_other();
	}
	void remove_parameters(Node* prev_state_on_structure, Node* state_on_structure, Node* next_state_on_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_on_structure == NULL){
			assert(state_on_structure != NULL);
			assert(state_on_structure->_transition_tssb != NULL);
			assert(next_state_on_structure != NULL);
			assert(is_node_on_structure_tssb(state_on_structure));
			assert(is_node_on_structure_tssb(next_state_on_structure));
			Node* state_on_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_bos);
			remove_customer_from_tssb_node(state_on_bos);
			remove_customer_from_tssb_node(state_on_structure);			// 参照カウント用
			Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
			assert(next_state_on_htssb != NULL);
			remove_customer_from_htssb_node(next_state_on_htssb);
			remove_customer_from_tssb_node(next_state_on_structure);		// 参照カウント用
			remove_customer_from_hpylm(state_on_structure, word_id);
			state_on_structure->decrement_transition_count_to_other();
			return;
		}
		// <eos>への遷移を含む場合
		if(next_state_on_structure == NULL){
			assert(prev_state_on_structure != NULL);
			assert(prev_state_on_structure->_transition_tssb != NULL);
			assert(state_on_structure != NULL);
			assert(is_node_on_structure_tssb(prev_state_on_structure));
			assert(is_node_on_structure_tssb(state_on_structure));
			Node* state_on_htssb = prev_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_htssb != NULL);
			remove_customer_from_htssb_node(state_on_htssb);
			remove_customer_from_tssb_node(state_on_structure);		// 参照カウント用
			remove_customer_from_hpylm(state_on_structure, word_id);
			prev_state_on_structure->decrement_transition_count_to_other();
			state_on_structure->decrement_transition_count_to_eos();
			return;
		}
		// <bos>と<eos>両方を含む場合
		if(prev_state_on_structure == NULL && next_state_on_structure == NULL){
			assert(state_on_structure != NULL);
			assert(is_node_on_structure_tssb(state_on_structure));
			Node* state_on_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
			assert(state_on_bos);
			remove_customer_from_tssb_node(state_on_bos);
			remove_customer_from_tssb_node(state_on_structure);			// 参照カウント用
			remove_customer_from_hpylm(state_on_structure, word_id);
			state_on_structure->decrement_transition_count_to_eos();
			return;
		}
		assert(prev_state_on_structure != NULL);
		assert(prev_state_on_structure->_transition_tssb != NULL);
		assert(state_on_structure != NULL);
		assert(state_on_structure->_transition_tssb != NULL);
		assert(next_state_on_structure != NULL);
		assert(is_node_on_structure_tssb(prev_state_on_structure));
		assert(is_node_on_structure_tssb(state_on_structure));
		assert(is_node_on_structure_tssb(next_state_on_structure));

		Node* state_on_htssb = prev_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_on_structure);
		assert(state_on_htssb != NULL);
		remove_customer_from_htssb_node(state_on_htssb);
		remove_customer_from_tssb_node(state_on_structure);			// 参照カウント用

		Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
		assert(next_state_on_htssb != NULL);
		remove_customer_from_htssb_node(next_state_on_htssb);
		remove_customer_from_tssb_node(next_state_on_structure);		// 参照カウント用

		remove_customer_from_hpylm(state_on_structure, word_id);
		prev_state_on_structure->decrement_transition_count_to_other();
		state_on_structure->decrement_transition_count_to_other();
	}
	// 新しい状態のギブスサンプリング
	// なるべく論文の記号を使う
	Node* draw_state(Node* prev_state_on_structure, Node* state_on_structure, Node* next_state_on_structure, id word_id){
		if(prev_state_on_structure == NULL){
			return _draw_state_from_bos(state_on_structure, next_state_on_structure, word_id);
		}
		if(next_state_on_structure == NULL){
			return _draw_state_to_eos(prev_state_on_structure, state_on_structure, word_id);
		}
		return _draw_state(prev_state_on_structure, state_on_structure, next_state_on_structure, word_id);
	}
	Node* _draw_state(Node* prev_state_on_structure, Node* state_on_structure, Node* next_state_on_structure, id word_id){
		assert(is_node_on_structure_tssb(state_on_structure));
		assert(is_node_on_structure_tssb(prev_state_on_structure));
		assert(is_node_on_structure_tssb(next_state_on_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_on_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		//// s_{t}から<eos>へ接続する確率
		double Peos_given_s = state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double stick_length = 1.0 - Peos_given_s;
		assert(state_on_structure->_transition_tssb != NULL);
		Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
		assert(next_state_on_htssb != NULL);
		double Pt_given_s = compute_node_probability_on_tssb(state_on_structure->_transition_tssb, next_state_on_htssb, stick_length);
		assert(0 < Pt_given_s && Pt_given_s <= 1);

		// スライス
		double slice = Pw_given_s * Pt_given_s * Sampler::uniform(0, 1);

		// s_{t-1}から<eos>へ接続する確率
		double Peos_given_ps = prev_state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double total_stick_length_of_prev_tssb = 1.0 - Peos_given_ps;
		double st = 0;
		double ed = total_stick_length_of_prev_tssb;


		// cout << "Pw_given_s: " << Pw_given_s << endl;
		// cout << "Peos_given_s: " << Peos_given_s << endl;
		// cout << "stick_length: " << stick_length << endl;
		// next_state_on_htssb->dump();
		// cout << "Pt_given_s: " << Pt_given_s << endl;
		// cout << "slice: " << slice << endl;
		// cout << "Peos_given_ps: " << Peos_given_ps << endl;
		// cout << "total_stick_length_of_prev_tssb: " << total_stick_length_of_prev_tssb << endl;
		
		while(true){
			double u = Sampler::uniform(st, ed);	// 最大値は棒の長さなので補正する
			Node* new_state_on_htssb = retrospective_sampling(u, prev_state_on_structure->_transition_tssb, total_stick_length_of_prev_tssb, true);
			assert(new_state_on_htssb != NULL);
			Node* new_state_on_structure = new_state_on_htssb->_structure_tssb_myself;
			assert(new_state_on_structure != NULL);
			assert(new_state_on_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_on_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			// 遷移確率
			//// s_{new}からs_{t+1}へ接続する確率
			Node* next_state_on_new_state_htssb = new_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
			//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
			double Peos_given_ns = new_state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
			double total_stick_length_of_new_tssb = 1.0 - Peos_given_ns;
			double new_Pt_given_s = compute_node_probability_on_tssb(new_state_on_structure->_transition_tssb, next_state_on_new_state_htssb, total_stick_length_of_new_tssb);
			assert(0 < new_Pt_given_s && new_Pt_given_s <= 1);
			// 尤度を計算
			double likelihoood = new_Pw_given_s * new_Pt_given_s;

			// cout << "u: " << u << endl;
			// new_state_on_htssb->dump();
			// new_state_on_structure->dump();
			// cout << "new_Pw_given_s: " << new_Pw_given_s << endl;
			// next_state_on_new_state_htssb->dump();
			// cout << "total_stick_length_of_new_tssb: " << total_stick_length_of_new_tssb << endl;
			// cout << "new_Pt_given_s: " << new_Pt_given_s << endl;
			// cout << "likelihoood: " << likelihoood << endl;

			// cout << likelihoood << ", " << slice << endl;
			// update_stick_length_of_tssb(prev_state_on_structure->_transition_tssb, total_stick_length_of_new_tssb, true);
			// prev_state_on_structure->_transition_tssb->dump();
			// state_on_structure->dump();
			// new_state_on_structure->dump();

			if(likelihoood > slice){
				return new_state_on_structure;
			}
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_on_structure, state_on_structure)){
				st = u;
			}else{
				ed = u;
			}
		}
	}
	Node* _draw_state_from_bos(Node* state_on_structure, Node* next_state_on_structure, id word_id){
		assert(is_node_on_structure_tssb(state_on_structure));
		assert(is_node_on_structure_tssb(next_state_on_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_on_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		//// s_{t}から<eos>へ接続する確率
		double Peos_given_s = state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double stick_length = 1.0 - Peos_given_s;
		assert(state_on_structure->_transition_tssb != NULL);
		Node* next_state_on_htssb = state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
		assert(next_state_on_htssb != NULL);
		double Pt_given_s = compute_node_probability_on_tssb(state_on_structure->_transition_tssb, next_state_on_htssb, stick_length);
		assert(0 < Pt_given_s && Pt_given_s <= 1);
		// スライス
		double slice = Pw_given_s * Pt_given_s * Sampler::uniform(0, 1);
		double st = 0;
		double ed = 1;	// ここは補正なし
		while(true){
			double u = Sampler::uniform(st, ed);
			Node* new_state_on_bos = retrospective_sampling(u, _bos_tssb, 1.0, false);
			assert(is_node_on_bos_tssb(new_state_on_bos));
			assert(new_state_on_bos != NULL);
			Node* new_state_on_structure = new_state_on_bos->_structure_tssb_myself;
			assert(new_state_on_structure != NULL);
			assert(new_state_on_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_on_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			// 遷移確率
			//// s_{new}からs_{t+1}へ接続する確率
			Node* next_state_on_new_state_htssb = new_state_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_on_structure);
			//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
			double Peos_given_ns = new_state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
			double total_stick_length_of_new_tssb = 1.0 - Peos_given_ns;
			double new_Pt_given_s = compute_node_probability_on_tssb(new_state_on_structure->_transition_tssb, next_state_on_new_state_htssb, total_stick_length_of_new_tssb);
			assert(0 < new_Pt_given_s && new_Pt_given_s <= 1);
			// 尤度を計算
			double likelihoood = new_Pw_given_s * new_Pt_given_s;
			if(likelihoood > slice){
				return new_state_on_structure;
			}
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_on_structure, state_on_structure)){
				st = u;
			}else{
				ed = u;
			}
		}
	}
	Node* _draw_state_to_eos(Node* prev_state_on_structure, Node* state_on_structure, id word_id){
		assert(is_node_on_structure_tssb(prev_state_on_structure));
		assert(is_node_on_structure_tssb(state_on_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_on_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		//// s_{t}から<eos>へ接続する確率
		double Peos_given_s = state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		assert(0 < Peos_given_s && Peos_given_s <= 1);
		// スライス
		double slice = Pw_given_s * Peos_given_s * Sampler::uniform(0, 1);
		// s_{t-1}から<eos>へ接続する確率
		double Peos_given_ps = prev_state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double total_stick_length_of_prev_tssb = 1.0 - Peos_given_ps;
		double st = 0;
		double ed = total_stick_length_of_prev_tssb;
		while(true){
			double u = Sampler::uniform(st, ed);	// 最大値は棒の長さなので補正する
			Node* new_state_on_htssb = retrospective_sampling(u, prev_state_on_structure->_transition_tssb, total_stick_length_of_prev_tssb, true);
			assert(new_state_on_htssb != NULL);
			Node* new_state_on_structure = new_state_on_htssb->_structure_tssb_myself;
			assert(new_state_on_structure != NULL);
			assert(new_state_on_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_on_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			// 遷移確率
			//// s_{new}から<eos>へ接続する確率
			double Peos_given_ns = new_state_on_structure->compute_transition_probability_to_eos(_tau0, _tau1);
			assert(0 < Peos_given_ns && Peos_given_ns <= 1);
			// 尤度を計算
			double likelihoood = new_Pw_given_s * Peos_given_ns;
			if(likelihoood > slice){
				return new_state_on_structure;
			}
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_on_structure, state_on_structure)){
				st = u;
			}else{
				ed = u;
			}
		}
	}
	void add_customer_to_hpylm(Node* target_on_structure, id token_id){
		assert(target_on_structure != NULL);
		assert(is_node_on_structure_tssb(target_on_structure));	// 木構造上のノードのみ
		assert(target_on_structure->_depth_v == target_on_structure->_hpylm->_depth);
		assert(_hpylm_d_m.size() > target_on_structure->_depth_v);
		assert(_hpylm_theta_m.size() > target_on_structure->_depth_v);
		assert(_word_g0 > 0);
		target_on_structure->_hpylm->add_customer(token_id, _word_g0, _hpylm_d_m, _hpylm_theta_m);
	}
	void add_customer_to_tssb_node(Node* target_on_tssb){
		assert(target_on_tssb != NULL);
		assert(is_node_on_htssb(target_on_tssb) == false);
		double alpha = _alpha * pow(_lambda, target_on_tssb->_depth_v);
		bool new_table_generated = false;
		// double ratio_v = compute_expectation_of_vertical_tssb_sbr_ratio(target_on_tssb);		// 特に計算しても意味はない
		target_on_tssb->add_customer_to_vertical_crp(alpha, 0, new_table_generated);
		// double ratio_h = compute_expectation_of_horizontal_tssb_sbr_ratio(target_on_tssb);		// 特に計算しても意味はない
		target_on_tssb->add_customer_to_horizontal_crp(_gamma, 0, new_table_generated);
		// 総客数のインクリメント
		if(is_node_on_structure_tssb(target_on_tssb)){
			_structure_tssb->increment_num_customers();
		}else if(is_node_on_bos_tssb(target_on_tssb)){
			_bos_tssb->increment_num_customers();
		}
		// 参照カウントのインクリメント
		//// <bos>からの接続のカウント
		if(is_node_on_bos_tssb(target_on_tssb)){
			Node* target_on_structure = target_on_tssb->_structure_tssb_myself;
			assert(target_on_structure != NULL);
			target_on_structure->increment_ref_count();
		}
	}
	void add_customer_to_htssb_node(Node* target_on_htssb){
		assert(target_on_htssb != NULL);
		assert(is_node_on_htssb(target_on_htssb));
		double alpha = _alpha * pow(_lambda, target_on_htssb->_depth_v);
		_add_customer_to_htssb_vertical_crp(alpha, target_on_htssb);
		_add_customer_to_htssb_horizontal_crp(_gamma, target_on_htssb);
	}
	void _add_customer_to_htssb_vertical_crp(double alpha, Node* iterator){
		assert(iterator != NULL);
		assert(is_node_on_htssb(iterator));
		bool new_table_generated = false;
		double ratio_v = compute_expectation_of_vertical_htssb_sbr_ratio(iterator);
		iterator->add_customer_to_vertical_crp(alpha, ratio_v, new_table_generated);
		// 総客数のインクリメント
		Node* owner_on_structure = iterator->_owner_on_structure;
		assert(owner_on_structure != NULL);
		TSSB* htssb = owner_on_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_on_structure);
		htssb->increment_num_customers();
		// 参照カウントのインクリメント
		Node* iterator_on_structure = iterator->_structure_tssb_myself;
		assert(iterator_on_structure != NULL);
		iterator_on_structure->increment_ref_count();
		// 親TSSBに代理客を追加
		Node* iterator_on_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(new_table_generated && iterator_on_parent_htssb != NULL){
			_add_customer_to_htssb_vertical_crp(alpha, iterator_on_parent_htssb);
		}
	}
	void _add_customer_to_htssb_horizontal_crp(double gamma, Node* iterator){
		assert(iterator != NULL);
		assert(is_node_on_htssb(iterator));
		bool new_table_generated = false;
		double ratio_h = compute_expectation_of_horizontal_htssb_sbr_ratio(iterator);
		iterator->add_customer_to_horizontal_crp(gamma, ratio_h, new_table_generated);
		// 総客数のインクリメント
		Node* owner_on_structure = iterator->_owner_on_structure;
		assert(owner_on_structure != NULL);
		TSSB* htssb = owner_on_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_on_structure);
		for(int d = 0;d <= iterator->_depth_v;d++){	// 水平方向には深さの数だけ別のSBRがあり、別の客が追加されることに注意
			htssb->increment_num_customers();
		}
		// 参照カウントのインクリメント
		Node* iterator_on_structure = iterator->_structure_tssb_myself;
		assert(iterator_on_structure != NULL);
		iterator_on_structure->increment_ref_count();
		// 親TSSBに代理客を追加
		Node* iterator_on_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(new_table_generated && iterator_on_parent_htssb != NULL){
			_add_customer_to_htssb_horizontal_crp(gamma, iterator_on_parent_htssb);
		}
	}
	void remove_customer_from_hpylm(Node* target_on_structure, id token_id){
		assert(target_on_structure != NULL);
		assert(is_node_on_structure_tssb(target_on_structure));	// 木構造上のノードのみ
		assert(target_on_structure->_depth_v == target_on_structure->_hpylm->_depth);
		target_on_structure->_hpylm->remove_customer(token_id);
	}
	void remove_customer_from_tssb_node(Node* target_on_tssb){
		assert(target_on_tssb != NULL);
		assert(is_node_on_htssb(target_on_tssb) == false);
		bool empty_table_deleted = false;
		target_on_tssb->remove_customer_from_vertical_crp(empty_table_deleted);
		target_on_tssb->remove_customer_from_horizontal_crp(empty_table_deleted);
		// 総客数のインクリメント
		if(is_node_on_structure_tssb(target_on_tssb)){
			_structure_tssb->decrement_num_customers();
		}else if(is_node_on_bos_tssb(target_on_tssb)){
			_bos_tssb->decrement_num_customers();
		}
		// 参照カウントのインクリメント
		//// <bos>からの接続のカウント
		if(is_node_on_bos_tssb(target_on_tssb)){
			Node* target_on_structure = target_on_tssb->_structure_tssb_myself;
			assert(target_on_structure != NULL);
			target_on_structure->decrement_ref_count();
		}
	}
	void remove_customer_from_htssb_node(Node* target_on_htssb){
		assert(target_on_htssb != NULL);
		assert(is_node_on_htssb(target_on_htssb));
		_remove_customer_from_htssb_vertical_crp(target_on_htssb);
		_remove_customer_from_htssb_horizontal_crp(target_on_htssb);
	}
	void _remove_customer_from_htssb_vertical_crp(Node* iterator){
		assert(iterator != NULL);
		assert(is_node_on_htssb(iterator));
		bool empty_table_deleted = false;
		iterator->remove_customer_from_vertical_crp(empty_table_deleted);
		// 総客数のインクリメント
		Node* owner_on_structure = iterator->_owner_on_structure;
		assert(owner_on_structure != NULL);
		TSSB* htssb = owner_on_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_on_structure);
		htssb->decrement_num_customers();
		// 参照カウントのインクリメント
		Node* iterator_on_structure = iterator->_structure_tssb_myself;
		assert(iterator_on_structure != NULL);
		iterator_on_structure->decrement_ref_count();
		// 親TSSBから代理客を削除
		Node* iterator_on_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(empty_table_deleted && iterator_on_parent_htssb != NULL){
			_remove_customer_from_htssb_vertical_crp(iterator_on_parent_htssb);
		}
	}
	void _remove_customer_from_htssb_horizontal_crp(Node* iterator){
		assert(iterator != NULL);
		bool empty_table_deleted = false;
		iterator->remove_customer_from_horizontal_crp(empty_table_deleted);
		// 総客数のインクリメント
		Node* owner_on_structure = iterator->_owner_on_structure;
		assert(owner_on_structure != NULL);
		TSSB* htssb = owner_on_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_on_structure);
		for(int d = 0;d <= iterator->_depth_v;d++){	// 水平方向には深さの数だけ別のSBRがあり、別の客が追加されることに注意
			htssb->decrement_num_customers();
		}
		// 参照カウントのインクリメント
		Node* iterator_on_structure = iterator->_structure_tssb_myself;
		assert(iterator_on_structure != NULL);
		iterator_on_structure->decrement_ref_count();
		// 親TSSBから代理客を削除
		Node* iterator_on_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(empty_table_deleted && iterator_on_parent_htssb != NULL){
			_remove_customer_from_htssb_horizontal_crp(iterator_on_parent_htssb);
		}
	}
	// update_stick_length_of_tssbは全ノードを更新するのに対しこっちは対象ノードのみ正確に計算する
	double compute_node_probability_on_tssb(TSSB* tssb, Node* node, double total_stick_length){
		assert(tssb->_owner_id == node->_owner_id_on_structure);
		bool htssb_mode = is_tssb_htssb(tssb);
		Node* iterator = tssb->_root;
		iterator->_stick_length = total_stick_length;
		double ratio_v = compute_expectation_of_vertical_sbr_ratio(iterator, htssb_mode);
		iterator->_probability = iterator->_stick_length * ratio_v;
		iterator->_children_stick_length = iterator->_stick_length * (1.0 - ratio_v);
		for(int n = 0;n < node->_depth_v;n++){
			int depth_h = node->_horizontal_indices_from_root[n];
			double rest_stick_length = iterator->_children_stick_length;
			for(int m = 0;m <= depth_h;m++){
				Node* child = iterator->_children[m];
				double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
				child->_stick_length = rest_stick_length * ratio_h;
				double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
				child->_probability = child->_stick_length * ratio_v;
				child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
				rest_stick_length *= (1.0 - ratio_h);
			}
			iterator = iterator->_children[depth_h];
		}
		return iterator->_probability;
	}
	double compute_expectation_of_vertical_sbr_ratio(Node* iterator, bool htssb_mode){
		if(htssb_mode){
			assert(is_node_on_htssb(iterator));
			return compute_expectation_of_vertical_htssb_sbr_ratio(iterator);
		}
		assert(is_node_on_htssb(iterator) == false);
		return compute_expectation_of_vertical_tssb_sbr_ratio(iterator);
	}
	double compute_expectation_of_horizontal_sbr_ratio(Node* iterator, bool htssb_mode){
		if(htssb_mode){
			assert(is_node_on_htssb(iterator));
			return compute_expectation_of_horizontal_htssb_sbr_ratio(iterator);
		}
		assert(is_node_on_htssb(iterator) == false);
		return compute_expectation_of_horizontal_tssb_sbr_ratio(iterator);
	}
	double compute_expectation_of_vertical_tssb_sbr_ratio(Node* target_on_tssb){
		int pass_count = target_on_tssb->_pass_count_v;
		int stop_count = target_on_tssb->_stop_count_v;
		double alpha = _alpha * pow(_lambda, target_on_tssb->_depth_v);
		return (1.0 + stop_count) / (1.0 + alpha + stop_count + target_on_tssb->_pass_count_v);
	}
	double compute_expectation_of_horizontal_tssb_sbr_ratio(Node* target_on_tssb){
		int pass_count = target_on_tssb->_pass_count_h;
		int stop_count = target_on_tssb->_stop_count_h;
		return (1.0 + stop_count) / (1.0 + _gamma + stop_count + pass_count);
	}
	// 縦の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double compute_expectation_of_vertical_htssb_sbr_ratio(Node* target_on_htssb){
		// c_printf("[*]%s\n", "compute_expectation_of_vertical_htssb_sbr_ratio");
		assert(target_on_htssb != NULL);
		assert(is_node_on_htssb(target_on_htssb));	// 木構造上のノードだった場合は計算できない
		Node* owner_on_structure = target_on_htssb->_owner_on_structure;
		assert(owner_on_structure != NULL);
		double sbr_ratio = -1;
		// target_on_htssb->dump();

		// 木構造上での基準となるノードを選ぶ
		assert(owner_on_structure != NULL);
		// 階層TSSBなので親TSSBのSBPを全て睿珊しないと次ノードのSBPを計算できない
		int num_itr_on_structure = owner_on_structure->_depth_v + 1;
		int num_itr_on_htssb = target_on_htssb->_depth_v + 1;
		// キャッシュ用配列
		double* stop_ratio_over_parent = target_on_htssb->_stop_ratio_v_over_parent;
		double* stop_probability_over_parent = target_on_htssb->_stop_probability_v_over_parent;
		// 計算
		for(int n = 0;n < num_itr_on_structure;n++){
			// cout << "n = " << n << endl;
			// 木構造上での基準となるノードを選ぶ
			// nが増えるごとに木構造を下に降りていく
			Node* iterator_on_structure = owner_on_structure->_nodes_from_root_to_myself[n];
			assert(iterator_on_structure != NULL);
			assert(iterator_on_structure->_transition_tssb != NULL);
			// iterator_on_structure->dump();
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			Node* iterator_on_htssb = iterator_on_structure->_transition_tssb->_root;
			assert(iterator_on_htssb != NULL);
			for(int m = 0;m < num_itr_on_htssb;m++){
				assert(iterator_on_htssb != NULL);
				// cout << "m = " << m << endl;
				// iterator_on_htssb->dump();
				if(n == 0){	// 木構造の親ノードの場合
					int pass_count = iterator_on_htssb->_pass_count_v;
					int stop_count = iterator_on_htssb->_stop_count_v;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					// cout << "depth = " << iterator_on_htssb->_depth_v << endl;
					double alpha = _alpha * pow(_lambda, iterator_on_htssb->_depth_v);
					// cout << "alpha = " << alpha << endl;
					double ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + iterator_on_htssb->_pass_count_v + EPS);
					// cout << "ratio_v = " << ratio_v << endl;
					assert(ratio_v < 1);
					stop_ratio_over_parent[m] = ratio_v;
					sbr_ratio = ratio_v;
				}else{	// 親の遷移確率用HTSSBから生成
					int pass_count = iterator_on_htssb->_pass_count_v;
					int stop_count = iterator_on_htssb->_stop_count_v;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					// cout << "depth = " << iterator_on_htssb->_depth_v << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					double alpha = _alpha * pow(_lambda, iterator_on_htssb->_depth_v);
					// cout << "alpha = " << alpha << endl;
					double ratio_v = (alpha * parent_stop_probability + stop_count) / (alpha * (1.0 - sum_parent_stop_probability) + stop_count + pass_count + EPS);
					// cout << "ratio_v = " << ratio_v << endl;
					assert(ratio_v < 1);
					stop_ratio_over_parent[m] = ratio_v;
					sum_parent_stop_probability += parent_stop_probability;
					sbr_ratio = ratio_v;
				}
				// HTSSB上で親から子へ降りていく
				if(m < num_itr_on_htssb - 1){
					int index_h = target_on_htssb->_horizontal_indices_from_root[m];
					assert(index_h < iterator_on_htssb->_children.size());
					iterator_on_htssb = iterator_on_htssb->_children[index_h];
				}
			}
			// 計算した棒を折る比率から確率を計算
			if(n < num_itr_on_structure - 1){
				double rest_stick_length = 1;
				for(int m = 0;m < target_on_htssb->_depth_v + 1;m++){
					// cout << "m = " << m << endl;
					double ratio_v = stop_ratio_over_parent[m];
					double stop_probability = rest_stick_length * ratio_v;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_v;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					stop_probability_over_parent[m] = stop_probability;
				}
			}
		}
		assert(sbr_ratio > 0);
		return sbr_ratio;
	}
	// 横の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double compute_expectation_of_horizontal_htssb_sbr_ratio(Node* target_on_htssb){
		// c_printf("[*]%s\n", "compute_expectation_of_horizontal_htssb_sbr_ratio");
		assert(target_on_htssb != NULL);
		if(target_on_htssb->_depth_v == 0){	// ルートノードなら必ず止まる
			return 1;
		}
		Node* owner_on_structure = target_on_htssb->_owner_on_structure;
		assert(owner_on_structure != NULL);
		int num_itr_on_structure = owner_on_structure->_depth_v + 1;
		// int depth_v = owner_on_structure->_depth_v;
		int num_itr_horizontal = target_on_htssb->_depth_h + 1;
		double sbr_ratio = 0;
		double* stop_ratio_over_parent = target_on_htssb->_stop_ratio_h_over_parent;
		double* stop_probability_over_parent = target_on_htssb->_stop_probability_h_over_parent;
		for(int n = 0;n < num_itr_on_structure;n++){
			// cout << "n = " << n << endl;
			// トップレベルのノードから順に停止確率を計算
			Node* iterator_on_structure = owner_on_structure->_nodes_from_root_to_myself[n];
			Node* iterator_on_htssb = iterator_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(target_on_htssb);
			assert(iterator_on_htssb != NULL);
			double sum_parent_stop_probability = 0;
			Node* parent_contains_target_on_htssb = iterator_on_htssb->_parent;
			assert(parent_contains_target_on_htssb != NULL);
			// cout << "parent contains" << endl << "	";
			// parent_contains_target_on_htssb->dump();
			for(int m = 0;m < num_itr_horizontal;m++){		// 自分自身も含めるので+1
				Node* child_on_htssb = parent_contains_target_on_htssb->_children[m];
				assert(child_on_htssb);
				// cout << "m = " << m << endl << "	";
				if(n == 0){		// 親ノードの場合
					int pass_count = child_on_htssb->_pass_count_h;
					int stop_count = child_on_htssb->_stop_count_h;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double ratio_h = (1.0 + stop_count) / (1.0 + _gamma + stop_count + pass_count + EPS);
					assert(ratio_h < 1);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
					sbr_ratio = ratio_h;
				}else{
					int pass_count = child_on_htssb->_pass_count_h;
					int stop_count = child_on_htssb->_stop_count_h;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					double ratio_h = (_gamma * parent_stop_probability + stop_count) / (_gamma * (1.0 - sum_parent_stop_probability) + stop_count + pass_count + EPS);
					assert(ratio_h < 1);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
					sum_parent_stop_probability += parent_stop_probability;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					sbr_ratio = ratio_h;
				}

			}
			if(n < num_itr_on_structure - 1){
				double rest_stick_length = 1;
				for(int m = 0;m < num_itr_horizontal;m++){
					double ratio_h = stop_ratio_over_parent[m];
					double stop_probability = rest_stick_length * ratio_h;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_h;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					stop_probability_over_parent[m] = stop_probability;
				}
			}
		}
		assert(sbr_ratio > 0);
		return sbr_ratio;
	}
	double compute_Pw_given_s(id token_id, Node* node_on_structure){
		assert(node_on_structure != NULL);
		assert(node_on_structure->_hpylm != NULL);
		assert(is_node_on_structure_tssb(node_on_structure));
		assert(_hpylm_d_m.size() > node_on_structure->_depth_v);
		assert(_hpylm_theta_m.size() > node_on_structure->_depth_v);
		assert(_word_g0 > 0);
		return node_on_structure->_hpylm->compute_Pw(token_id, _word_g0, _hpylm_d_m, _hpylm_theta_m);
	}
	void update_stick_length_of_tssb(TSSB* tssb, double total_stick_length, bool htssb_mode){
		// assert(tssb->_owner_id != 0);	// 木構造の場合は計算しない
		Node* root = tssb->_root;
		double ratio_v = compute_expectation_of_vertical_sbr_ratio(root, htssb_mode);
		double sum_probability = total_stick_length * ratio_v;
		root->_stick_length = total_stick_length;
		root->_children_stick_length = total_stick_length * (1.0 - ratio_v);
		root->_probability = ratio_v * total_stick_length;
		root->_sum_probability = root->_probability;
		_update_stick_length_of_parent_node(sum_probability, root, htssb_mode);
	}
	void _update_stick_length_of_parent_node(double &sum_probability, Node* parent, bool htssb_mode){
		assert(parent->_children_stick_length > 0);
		double rest_stick_length = parent->_children_stick_length;	// 親ノードが持っている子ノードに割り当てる棒の長さ
		double sum_stick_length_over_children = sum_probability;
		for(int i = 0;i < parent->_children.size();i++){
			Node* child = parent->_children[i];
			double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			child->_stick_length = rest_stick_length * ratio_h;		// このノードかこのノードの子ノードに止まる確率
			double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
			child->_probability = child->_stick_length * ratio_v;	// このノードに止まる確率
			sum_probability += child->_probability;
			sum_stick_length_over_children += child->_stick_length;
			child->_sum_probability = sum_stick_length_over_children - child->_children_stick_length;				// このノードより左側の全ての棒の長さの総和
			rest_stick_length *= 1.0 - ratio_h;
			double alpha = _alpha * pow(_lambda, child->_depth_v);
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);	// 子ノードに割り当てる長さはこのノードに降りない確率
			if(child->has_child()){
				_update_stick_length_of_parent_node(child->_sum_probability, child, htssb_mode);
			}
		}
	}
	void delete_invalid_children_on_structure_tssb(TSSB* tssb){
		assert(is_tssb_structure(tssb));
		delete_invalid_children_of_node_on_structure(tssb->_root);
	}
	void delete_invalid_children_of_node_on_structure(Node* parent){
		assert(is_node_on_structure_tssb(parent));
		vector<Node*> &children = parent->_children;
		for(int i = children.size() - 1;i >= 0;i--){
			Node* child = children[i];
			delete_invalid_children_of_node_on_structure(child);
			bool success = delete_node_on_structure_if_needed(child);
			if(success == false){	// 失敗したらそれ以上は消さない
				// break;
			}
		}
	}
	bool delete_node_on_structure_if_needed(Node* target_on_structure){
		assert(target_on_structure != NULL);
		assert(is_node_on_structure_tssb(target_on_structure));
		if(target_on_structure->_depth_v == 0){
			return false;
		}
		assert(target_on_structure->_parent != NULL);
		int delete_id = target_on_structure->_identifier;
		if(target_on_structure->_pass_count_v != 0){
			return false;
		}
		if(target_on_structure->_stop_count_v != 0){
			return false;
		}
		if(target_on_structure->_pass_count_h != 0){
			return false;
		}
		if(target_on_structure->_stop_count_h != 0){
			return false;
		}
		if(target_on_structure->_ref_count != 0){
			return false;
		}
		if(target_on_structure->_num_transitions_to_eos != 0){
			return false;
		}
		if(target_on_structure->_num_transitions_to_other != 0){
			return false;
		}
		Node* parent_on_structure = target_on_structure->_parent;
		Node* delete_node = target_on_structure->_parent->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		// 全てのHTSSBから削除
		_delete_node_on_all_htssb(delete_id, _structure_tssb->_root, parent_on_structure);
		// <bos>から削除
		Node* parent_on_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(parent_on_structure);
		assert(parent_on_bos != NULL);
		assert(is_node_on_bos_tssb(parent_on_bos));
		delete_node = parent_on_bos->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		return true;
	}
	void _delete_node_on_all_htssb(int delete_id, Node* iterator_on_structure, Node* target_parent_on_structure){
		assert(target_parent_on_structure != NULL);
		assert(iterator_on_structure->_transition_tssb != NULL);
		// 遷移確率用TSSBでの同じ位置の子ノードを削除
		Node* parent_on_htssb = iterator_on_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(target_parent_on_structure);
		assert(parent_on_htssb != NULL);
		assert(is_node_on_htssb(parent_on_htssb));
		Node* delete_node = parent_on_htssb->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		for(const auto &child: iterator_on_structure->_children){
			_delete_node_on_all_htssb(delete_id, child, target_parent_on_structure);
		}
	}
	// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
	// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
	void sum_auxiliary_variables_recursively_for_hpylm(Node* parent, vector<double> &sum_log_x_u_m, vector<double> &sum_y_ui_m, vector<double> &sum_1_y_ui_m, vector<double> &sum_1_z_uwkj_m){
		for(const auto &child: parent->_children){
			int depth = child->_depth_v;
			assert(depth < _hpylm_d_m.size());
			assert(depth < _hpylm_theta_m.size());

			double d = _hpylm_d_m[depth];
			double theta = _hpylm_theta_m[depth];
			sum_log_x_u_m[depth] += child->_hpylm->auxiliary_log_x_u(theta);	// log(x_u)
			sum_y_ui_m[depth] += child->_hpylm->auxiliary_y_ui(d, theta);		// y_ui
			sum_1_y_ui_m[depth] += child->_hpylm->auxiliary_1_y_ui(d, theta);	// 1 - y_ui
			sum_1_z_uwkj_m[depth] += child->_hpylm->auxiliary_1_z_uwkj(d);		// 1 - z_uwkj

			sum_auxiliary_variables_recursively_for_hpylm(child, sum_log_x_u_m, sum_y_ui_m, sum_1_y_ui_m, sum_1_z_uwkj_m);
		}
	}
	// dとθの推定
	void sample_hpylm_hyperparameters(){
		assert(_max_depth < _hpylm_d_m.size());
		assert(_max_depth < _hpylm_theta_m.size());
		assert(_max_depth < _hpylm_a_m.size());
		assert(_max_depth < _hpylm_b_m.size());
		assert(_max_depth < _hpylm_alpha_m.size());
		assert(_max_depth < _hpylm_beta_m.size());

		// 親ノードの深さが0であることに注意
		vector<double> sum_log_x_u_m(_max_depth + 1, 0.0);
		vector<double> sum_y_ui_m(_max_depth + 1, 0.0);
		vector<double> sum_1_y_ui_m(_max_depth + 1, 0.0);
		vector<double> sum_1_z_uwkj_m(_max_depth + 1, 0.0);

		// _root
		HPYLM* root = _structure_tssb->_root->_hpylm;
		sum_log_x_u_m[0] = root->auxiliary_log_x_u(_hpylm_theta_m[0]);			// log(x_u)
		sum_y_ui_m[0] = root->auxiliary_y_ui(_hpylm_d_m[0], _hpylm_theta_m[0]);			// y_ui
		sum_1_y_ui_m[0] = root->auxiliary_1_y_ui(_hpylm_d_m[0], _hpylm_theta_m[0]);		// 1 - y_ui
		sum_1_z_uwkj_m[0] = root->auxiliary_1_z_uwkj(_hpylm_d_m[0]);				// 1 - z_uwkj

		// それ以外
		sum_auxiliary_variables_recursively_for_hpylm(_structure_tssb->_root, sum_log_x_u_m, sum_y_ui_m, sum_1_y_ui_m, sum_1_z_uwkj_m);

		for(int u = 0;u <= _max_depth;u++){
			_hpylm_d_m[u] = Sampler::beta(_hpylm_a_m[u] + sum_1_y_ui_m[u], _hpylm_b_m[u] + sum_1_z_uwkj_m[u]);
			_hpylm_theta_m[u] = Sampler::gamma(_hpylm_alpha_m[u] + sum_y_ui_m[u], _hpylm_beta_m[u] - sum_log_x_u_m[u]);
		}
		// 不要な深さのハイパーパラメータを削除
		int num_remove = _hpylm_d_m.size() - _max_depth - 1;
		for(int n = 0;n < num_remove;n++){
			_hpylm_d_m.pop_back();
			_hpylm_theta_m.pop_back();
			_hpylm_a_m.pop_back();
			_hpylm_b_m.pop_back();
			_hpylm_alpha_m.pop_back();
			_hpylm_beta_m.pop_back();
		}
	}
	bool save(string dir = "out"){
		bool success = false;
		ofstream ofs(dir + "/ithmm.model");
		if(ofs.good()){
			boost::archive::binary_oarchive oarchive(ofs);
			oarchive << static_cast<const iTHMM&>(*this);
			success = true;
		}
		ofs.close();
		return success;
	}
	bool load(string dir = "out"){
		bool success = false;
		ifstream ifs(dir + "/ithmm.model");
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> *this;
			assert(_structure_tssb != NULL);
			assert(_structure_tssb->_root != NULL);
			vector<Node*> nodes;
			_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
			for(auto node: nodes){
				// 配列を確保
				node->init_arrays();
				node->init_horizontal_indices();
				node->init_pointers_from_root_to_myself();
				vector<Node*> nodes_on_htssb;
				node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_on_htssb);
				for(auto node_on_htssb: nodes_on_htssb){
					// 配列を確保
					node_on_htssb->init_arrays();
					node_on_htssb->init_horizontal_indices();
					node_on_htssb->init_pointers_from_root_to_myself();
				}
				vector<Node*>().swap(nodes_on_htssb);	// 解放
			}
			nodes.clear();
			vector<Node*>().swap(nodes);				// 解放
			_bos_tssb->enumerate_nodes_from_left_to_right(nodes);
			for(auto node: nodes){
				// 配列を確保
				node->init_arrays();
				node->init_horizontal_indices();
				node->init_pointers_from_root_to_myself();
			}
			vector<Node*>().swap(nodes);				// 解放
			success = true;
		}
		ifs.close();

		return success;
	}
};

#endif