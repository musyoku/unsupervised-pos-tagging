#pragma once
#include <boost/serialization/serialization.hpp>
#include <unordered_map>
#include <map>
#include "common.h"
#include "table.h"

namespace ithmm {
	class HPYLM;
	class TSSB;
	class Node {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, unsigned int version);
	public:
		static int _auto_increment;
		int _identifier;		// ノードID
		Node* _parent;			// 親ノード
		int _depth_v;			// 縦の深さ。 論文中の|s|に相当
		int _depth_h;			// 横の深さ。 論文中のkに相当
		int _pass_count_v;		// 通過回数。 縦方向のCDP
		int _stop_count_v;		// 停止回数。 縦方向のCDP
		int _pass_count_h;		// 通過回数。 横方向のCDP
		int _stop_count_h;		// 停止回数。 横方向のCDP
		int _num_transitions_to_eos;	// EOSへの遷移回数
		int _num_transitions_to_other;	// EOS以外への遷移回数
		Table* _table_v;		// 客を管理するテーブル。 縦方向のCRP
		Table* _table_h;		// 客を管理するテーブル。 横方向のCRP
		std::vector<Node*> _children;
		double _stick_length;				// 自分の棒の木全体に対する長さ
		double _children_stick_length;		// 自分の棒の子ノードに割り当てる長さ
		double _probability;				// このノードの確率
		double _sum_probability;			// 自分より左側にある全ての棒の長さの合計
		HPYLM* _hpylm;						// 出力分布
		std::unordered_map<int, int> _num_word_assignment;	// 単語がこのノードに割り当てられた回数。結果表示用でiTHMMとは無関係。
		// 計算時に使う配列
		Node** _nodes_from_root_to_myself;
		double* _stop_probability_v_over_parent;
		double* _stop_ratio_v_over_parent;
		int* _horizontal_indices_from_root;		// 親から辿ってこのノードに到達するための水平方向のインデックス
		double* _stop_probability_h_over_parent;
		double* _stop_ratio_h_over_parent;
		int _ref_count;							// 参照カウント
		// ポインタを張る
		//// 木構造上のノードの場合は以下のみ有効
		TSSB* _transition_tssb;					// 遷移確率を表すTSSBのルート
		Node* _myself_in_transition_tssb;			// 遷移確率を表すTSSBの自分と同じ位置のノード
		Node* _myself_in_bos_tssb;					// <bos>TSSBでの自分と同じ位置のノード
		//// HTSSB上のノードの場合は以下のみ有効
		Node* _myself_in_parent_transition_tssb;	// 木構造上の親ノードが持つ遷移確率TSSBの自分と同じ位置のノード
		Node* _myself_in_structure_tssb;			// 木構造上の自分と同じ位置のノード
		bool _is_structure_node;
		bool _is_htssb_node;
		bool _is_bos_tssb_node;
		Node* _htssb_owner_node_in_structure;	// このHTSSBを木構造のどのノードが持っているか
		Node();
		Node(Node* parent);
		Node(Node* parent, int identifier);
		void init();
		void init_arrays();
		void init_horizontal_indices();
		void init_pointers_from_root_to_myself();
		void init_hpylm();
		~Node();
		Node* generate_child();
		void add_child(Node* node);
		Node* find_same_node_in_transition_tssb();
		int get_htssb_owner_node_id();
		int get_vertical_stop_count();
		int get_vertical_pass_count();
		int get_horizontal_stop_count();
		int get_horizontal_pass_count();
		Table* get_vertical_table();
		Table* get_horizontal_table();
		TSSB* get_transition_tssb();
		Node* get_myself_in_transition_tssb();
		Node* get_myself_in_bos_tssb();
		Node* get_myself_in_parent_transition_tssb();
		Node* get_myself_in_structure_tssb();
		Node* get_htssb_owner_node_in_structure();
		void set_transition_tssb(TSSB* tssb);
		void set_myself_in_transition_tssb(Node* node);
		void set_myself_in_bos_tssb(Node* node);
		void set_myself_in_parent_transition_tssb(Node* node);
		void set_myself_in_structure_tssb(Node* node);
		void set_htssb_owner_node_in_structure(Node* node);
		bool is_structure_node();
		bool is_htssb_node();
		bool is_bos_tssb_node();
		void set_as_structure_node();
		void set_as_htssb_node();
		void set_as_bos_tssb_node();
		double compute_transition_probability_to_eos(double tau0, double tau1);
		bool has_child();
		bool has_parent();
		void add_customer_to_vertical_crp(double concentration, double g0, bool &new_table_generated);
		void increment_vertical_stop_count();
		void decrement_vertical_stop_count();
		void increment_vertical_pass_count();
		void decrement_vertical_pass_count();
		void add_customer_to_horizontal_crp(double concentration, double g0, bool &new_table_generated);
		void increment_horizontal_stop_count();
		void decrement_horizontal_stop_count();
		void increment_horizontal_pass_count();
		void decrement_horizontal_pass_count();
		void increment_transition_count_to_eos();
		void decrement_transition_count_to_eos();
		void increment_transition_count_to_other();
		void decrement_transition_count_to_other();
		void increment_ref_count();
		void decrement_ref_count();
		void increment_word_assignment(int word_id);
		void decrement_word_assignment(int word_id);
		// 客を除去
		void remove_customer_from_vertical_crp(bool &empty_table_deleted);
		void remove_last_customer_from_vertical_crp(bool &empty_table_deleted);
		void _remove_customer_from_vertical_crp(bool remove_last_customer, bool &empty_table_deleted);
		void remove_customer_from_horizontal_crp(bool &empty_table_deleted);
		void remove_last_customer_from_horizontal_crp(bool &empty_table_deleted);
		void _remove_customer_from_horizontal_crp(bool remove_last_customer, bool &empty_table_deleted);
		bool delete_node_if_needed();
		Node* delete_child_node(int node_id);
		static void _delete_all_children(Node* parent);
		void dump();
		std::string _dump_indices();
		std::wstring _wdump_indices();
		std::string _dump();
	};
}