#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <unordered_map>
#include <map>
#include "table.hpp"
#include "hpylm.hpp"

class HPYLM;
class TSSB;
class Node{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _auto_increment;;
		archive & _identifier;
		archive & _owner_id_on_structure;
		archive & _owner_on_structure;
		archive & _parent;
		archive & _depth_v;
		archive & _depth_h;
		archive & _pass_count_v;
		archive & _stop_count_v;
		archive & _pass_count_h;
		archive & _stop_count_h;
		archive & _num_word_assignment;
		archive & _num_transitions_to_eos;
		archive & _num_transitions_to_other;
		archive & _table_v;
		archive & _table_h;
		archive & _children;
		archive & _stick_length;
		archive & _children_stick_length;
		archive & _probability;
		archive & _sum_probability;
		archive & _hpylm;
		archive & _transition_tssb;
		archive & _transition_tssb_myself;
		archive & _parent_transition_tssb_myself;
		archive & _ref_count;
		archive & _structure_tssb_myself;
		archive & _bos_tssb_myself;
	}
public:
	static int _auto_increment;
	int _identifier;		// ノードID
	int _owner_id_on_structure;
	Node* _owner_on_structure;		// このHTSSBは木構造のどのノードが持っているか
	Node* _parent;			// 親ノード
	int _depth_v;			// 縦の深さ。 論文中の|s|に相当
	int _depth_h;			// 横の深さ。 論文中のkに相当
	// 各ノードの遷移確率TSSBは自己同型になっている必要があるため、構造を共有する
	// カウントは各ノードのIDごとに管理
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
	std::map<id, int> _num_word_assignment;	// 単語がこのノードに割り当てられた回数。結果表示用でiTHMMとは無関係。
	// 計算時に使う配列
	Node** _nodes_from_root_to_myself;
	double* _stop_probability_v_over_parent;
	double* _stop_ratio_v_over_parent;
	int* _horizontal_indices_from_root;		// 親から辿ってこのノードに到達するための水平方向のインデックス
	double* _stop_probability_h_over_parent;
	double* _stop_ratio_h_over_parent;
	// ポインタを張る
	//// 木構造上のノードの場合は以下のみ有効
	TSSB* _transition_tssb;					// 遷移確率を表すTSSBのルート
	Node* _transition_tssb_myself;			// 遷移確率を表すTSSBの自分と同じ位置のノード
	Node* _bos_tssb_myself;					// <bos>TSSBでの自分と同じ位置のノード
	int _ref_count;							// 参照カウント
	//// HTSSB上のノードの場合は以下のみ有効
	Node* _parent_transition_tssb_myself;	// 木構造上の親ノードが持つ遷移確率TSSBの自分と同じ位置のノード
	Node* _structure_tssb_myself;			// 木構造上の自分と同じ位置のノード
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
	Node* find_same_node_on_transition_tssb();
	int get_vertical_stop_count();
	int get_vertical_pass_count();
	int get_horizontal_stop_count();
	int get_horizontal_pass_count();
	Table* get_vertical_table();
	Table* get_horizontal_table();
	double compute_transition_probability_to_eos(double tau0, double tau1);
	bool has_child();
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
	void increment_word_assignment(id word_id);
	void decrement_word_assignment(id word_id);
	// 客を除去
	void remove_customer_from_vertical_crp(bool &empty_table_deleted);
	void remove_last_customer_from_vertical_crp(bool &empty_table_deleted);
	void _remove_customer_from_vertical_crp(bool remove_last_customer, bool &empty_table_deleted);
	void remove_customer_from_horizontal_crp(bool &empty_table_deleted);
	void remove_last_customer_from_horizontal_crp(bool &empty_table_deleted);
	void _remove_customer_from_horizontal_crp(bool remove_last_customer, bool &empty_table_deleted);
	bool delete_node_if_needed();
	Node* delete_child_node(int node_id);
	void dump();
	std::string _dump_indices();
	std::wstring _wdump_indices();
	std::string _dump();
};
