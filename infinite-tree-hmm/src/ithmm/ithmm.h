#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <vector>
#include <set>
#include <fstream>
#include "tssb.h"
#include "node.h"
#include "hpylm.h"
#include "sampler.h"
#include "cprintf.h"
#include "utils.h"
#include "common.h"
#include "hyperparameters.h"

// 10以下ならなんでもいい
#define TSSB_STRUCTURE_ID 8
#define TSSB_BOS_ID 7

namespace ithmm {
	struct multiset_value_comparator {
		bool operator()(const std::pair<id, double> &a, const std::pair<id, double> &b) {
			return a.second > b.second;
		}   
	};

	class iTHMM {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version);
	public:
		TSSB* _structure_tssb;	// 木構造を表すためだけのTSSB。HTSSBは全てこのノードを基準に成形する
		TSSB* _bos_tssb;		// <bos>からの遷移を表すTSSB
		double _alpha;
		double _gamma;
		double _lambda_alpha;	// 縦のTSSBの減衰率. 論文のlambdaに該当
		double _lambda_gamma;	// 横のTSSBの減衰率
		double _strength;		// HTSSBにおいて親の情報をどの程度受け継ぐかをコントロール
		double _tau0;
		double _tau1;
		double _word_g0;
		int _current_max_depth;
		int _depth_limit;		// -1なら無限大
		std::vector<double> _hpylm_d_m;		// HPYLMのハイパーパラメータ（ディスカウント係数）
		std::vector<double> _hpylm_theta_m;	// HPYLMのハイパーパラメータ（集中度）
		std::vector<double> _hpylm_a_m;		// ベータ分布のパラメータ	dの推定用
		std::vector<double> _hpylm_b_m;		// ベータ分布のパラメータ	dの推定用
		std::vector<double> _hpylm_alpha_m;	// ガンマ分布のパラメータ	θの推定用
		std::vector<double> _hpylm_beta_m;	// ガンマ分布のパラメータ	θの推定用
		bool _mh_enabled;				// メトロポリス・ヘイスティングス法による補正を行うかどうか
		// 統計
		int _num_mh_acceptance;
		int _num_mh_rejection;
		iTHMM();
		~iTHMM();
		void initialize_with_training_dataset(std::vector<std::vector<Word*>> &dataset);
		void remove_all_data(std::vector<std::vector<Word*>> &dataset);
		void set_depth_limit(int limit);
		void set_word_g0(double g0);
		bool is_node_in_bos_tssb(Node* node);
		bool is_node_in_structure_tssb(Node* node);
		bool is_node_in_htssb(Node* node);
		bool is_node_root(Node* node);
		bool is_tssb_bos(TSSB* tssb);
		bool is_tssb_structure(TSSB* tssb);
		bool is_tssb_htssb(TSSB* tssb);
		bool is_node_to_the_left_of_node(Node* left, Node* right);
		Node* generate_and_add_new_child_to(Node* parent);
		void _generate_and_add_new_child_to_all_htssb(Node* iterator_in_structure, Node* parent, Node* generated_child_in_structure, Node* &return_child);
		Node* _generate_and_add_new_child_to_bos_tssb(Node* generated_child_in_structure);
		TSSB* generate_transition_tssb_belonging_to(Node* owner_in_structure);
		void copy_children_in_structure_to_transition_tssb(Node* source_in_structure, Node* target_in_htssb, Node* owner_in_structure);
		Node* sample_node_in_tssb(TSSB* tssb, bool ignore_root = false);
		Node* sample_node_in_htssb(TSSB* tssb, bool ignore_root = false);
		Node* _sample_node_in_tssb_by_iterating_node(Node* iterator, bool htssb_mode, bool ignore_root = false);
		Node* retrospective_sampling(double uniform, TSSB* tssb, double total_stick_length, bool htssb_mode);
		Node* _retrospective_sampling_by_iterating_node(double uniform, double &sum_probability, Node* iterator, bool htssb_mode);
		void perform_gibbs_sampling_data(std::vector<Word*> &data);
		void add_initial_parameters(Node* prev_state_in_structure, Node* state_in_structure, id word_id);
		void add_temporal_parameters(Node* prev_state_in_structure, Node* state_in_structure);
		void add_parameters(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id);
		void remove_initial_parameters(Node* prev_state_in_structure, Node* state_in_structure, id word_id);
		void remove_temporal_parameters(Node* prev_state_in_structure, Node* state_in_structure);
		void remove_parameters(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id);
		Node* draw_state(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id);
		Node* _draw_state(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id);
		Node* _draw_state_from_bos(Node* state_in_structure, Node* next_state_in_structure, id word_id);
		Node* _draw_state_to_eos(Node* prev_state_in_structure, Node* state_in_structure, id word_id);
		void add_customer_to_hpylm(Node* target_in_structure, id token_id);
		void add_customer_to_tssb_node(Node* target_in_tssb);
		void add_customer_to_htssb_node(Node* target_in_htssb);
		void _add_customer_to_htssb_vertical_crp(double alpha, Node* iterator);
		void _add_customer_to_htssb_horizontal_crp(double gamma, Node* iterator);
		void remove_customer_from_hpylm(Node* target_in_structure, id token_id);
		void remove_customer_from_tssb_node(Node* target_in_tssb);
		void remove_customer_from_htssb_node(Node* target_in_htssb, bool remove_last_customer = false);
		void _remove_customer_from_htssb_vertical_crp(Node* iterator, bool remove_last_customer = false);
		void _remove_customer_from_htssb_horizontal_crp(Node* iterator, bool remove_last_customer = false);
		double compute_node_probability_in_tssb(TSSB* tssb, Node* node, double total_stick_length);
		double compute_expectation_of_vertical_sbr_ratio(Node* iterator, bool htssb_mode);
		double compute_expectation_of_horizontal_sbr_ratio(Node* iterator, bool htssb_mode);
		double compute_expectation_of_vertical_tssb_sbr_ratio(Node* target_in_tssb);
		double compute_expectation_of_horizontal_tssb_sbr_ratio(Node* target_in_tssb);
		double compute_expectation_of_vertical_htssb_sbr_ratio(Node* target_in_htssb);
		double compute_expectation_of_horizontal_htssb_sbr_ratio(Node* target_in_htssb);
		double compute_Pw_given_s(id token_id, Node* node_in_structure);
		void update_stick_length_of_tssb(TSSB* tssb, double total_stick_length, bool htssb_mode);
		void _update_stick_length_of_parent_node(double &sum_probability, Node* parent, bool htssb_mode);
		void delete_invalid_children();
		void _delete_invalid_children_in_structure_tssb(TSSB* tssb);
		void _delete_invalid_children_of_node_in_structure(Node* parent);
		bool delete_node_in_structure_if_needed(Node* target_in_structure);
		void _delete_node_in_all_htssb(int delete_id, Node* iterator_in_structure, Node* target_parent_in_structure);
		void sum_auxiliary_variables_recursively_for_hpylm(Node* parent, std::vector<double> &sum_log_x_u_m, std::vector<double> &sum_y_ui_m, std::vector<double> &sum_1_y_ui_m, std::vector<double> &sum_1_z_uwkj_m);
		void sample_hpylm_hyperparameters();
		void geneerate_word_ranking_of_node(Node* node_in_structure, std::multiset<std::pair<id, double>, multiset_value_comparator> &ranking);
		bool save(std::string filename);
		bool load(std::string filename);
	};
}