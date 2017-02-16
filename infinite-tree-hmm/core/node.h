#ifndef _node_
#define _node_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include "cprintf.h"
#include "sampler.h"
#include "util.h"
using namespace std;
using namespace boost;
#define CLUSTERING_TSSB_ID 0

// 中華料理店過程のテーブル
// 通常CRPではテーブルが各クラスタを表すが、TSSBでは全テーブルが同じクラスタに属する
class Table{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _arrangement;
		archive & _num_customers;
		archive & _token_id;
	}
public:
	vector<int> _arrangement;
	int _num_customers;
	int _token_id;
	Table(){
		_num_customers = 0;
		_token_id = 0;
	}
	Table(int token_id){
		_num_customers = 0;
		_token_id = token_id;
	}
	bool is_empty(){
		return _arrangement.size() == 0;
	}
	void add_customer(double concentration_parameter, bool &new_table_generated){
		_num_customers += 1;
		if(_arrangement.size() == 0){
			_arrangement.push_back(1);
			new_table_generated = true;
			return;
		}
		new_table_generated = false;
		double sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0) + concentration_parameter;
		double normalizer = 1 / sum;
		double bernoulli = Sampler::uniform(0, 1);
		sum = 0;
		for(int i = 0;i < _arrangement.size();i++){
			sum += _arrangement[i] * normalizer;
			if(bernoulli <= sum){
				_arrangement[i] += 1;
				return;
			}
		}
		_arrangement.push_back(1);
		new_table_generated = true;
	}
	void remove_customer(bool &empty_table_deleted){
		assert(_arrangement.size() > 0);
		empty_table_deleted = false;
		_num_customers -= 1;
		int sum = std::accumulate(_arrangement.begin(), _arrangement.end(), 0);
		int bernoulli = Sampler::uniform_int(0, sum);
		sum = 0;
		int target_index = _arrangement.size() - 1;
		for(int i = 0;i < _arrangement.size();i++){
			sum += _arrangement[i];
			if(bernoulli <= sum){
				target_index = i;
				break;
			}
		}
		_arrangement[target_index] -= 1;
		if(_arrangement[target_index] == 0){
			_arrangement.erase(_arrangement.begin() + target_index);
			empty_table_deleted = true;
		}
	}
};

class Node{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _identifier;
		archive & _parent;
		archive & _depth_v;
		archive & _depth_h;
		archive & _pass_count_v;
		archive & _stop_count_v;
		archive & _pass_count_h;
		archive & _stop_count_h;
		archive & _table_v;
		archive & _table_h;
	}
public:
	static int _auto_increment;
	int _identifier;	// ノードID
	Node* _parent;		// 親ノード
	int _depth_v;		// 縦の深さ。 論文中の|s|に相当
	int _depth_h;		// 横の深さ。 論文中のkに相当
	// 各ノードの遷移確率TSSBは自己同型になっている必要があるため、構造を共有する
	// カウントは各ノードのIDごとに管理
	unordered_map<int, int> _pass_count_v;	// 通過回数。 縦方向のCDP
	unordered_map<int, int> _stop_count_v;	// 停止回数。 縦方向のCDP
	unordered_map<int, int> _pass_count_h;	// 通過回数。 横方向のCDP
	unordered_map<int, int> _stop_count_h;	// 停止回数。 横方向のCDP
	unordered_map<int, Table*> _table_v;	// 客を管理するテーブル。 縦方向のCRP
	unordered_map<int, Table*> _table_h;	// 客を管理するテーブル。 横方向のCRP
	vector<Node*> _children;
	double _stick_length;					// 自分の棒の木全体に対する長さ
	double _children_stick_length;			// 自分の棒の子ノードに割り当てる長さ
	Node(Node* parent){
		_identifier = _auto_increment;
		_auto_increment++;
		_parent = parent;
		_depth_v = (parent != NULL) ? parent->_depth_v + 1 : 0;
		_stick_length = -1;
		_children_stick_length = -1;
	}
	~Node(){
		for(auto &elem: _table_v){
			delete elem.second;
		}
		for(auto &elem: _table_h){
			delete elem.second;
		}
	}
	void init_table(){
		Node* parent = _parent;
		while(parent != NULL){
			Table* table_v = new Table();
			_table_v[parent->_identifier] = table_v;
			Table* table_h = new Table();
			_table_h[parent->_identifier] = table_h;
			parent = parent->_parent;
		}
	}
	Node* generate_child(){
		Node* child = new Node(this);
		_children.push_back(child);
		return child;
	}
	int get_vertical_stop_count_with_id(int identifier){
		auto itr = _stop_count_v.find(identifier);
		if(itr == _stop_count_v.end()){
			return 0;
		}
		return itr->second;
	}
	int get_vertical_pass_count_with_id(int identifier){
		auto itr = _pass_count_v.find(identifier);
		if(itr == _pass_count_v.end()){
			return 0;
		}
		return itr->second;
	}
	int get_horizontal_stop_count_with_id(int identifier){
		auto itr = _stop_count_h.find(identifier);
		if(itr == _stop_count_h.end()){
			return 0;
		}
		return itr->second;
	}
	int get_horizontal_pass_count_with_id(int identifier){
		auto itr = _pass_count_h.find(identifier);
		if(itr == _pass_count_h.end()){
			return 0;
		}
		return itr->second;
	}
	Table* get_vertical_table(int identifier, bool generate_if_not_exist){
		auto itr = _table_v.find(identifier);
		if(itr == _table_v.end()){
			if(generate_if_not_exist){
				Table* table = new Table();
				_table_v[identifier] = table;
				return table;
			}
			return NULL;
		}
		return itr->second;
	}
	Table* get_htssb_horizontal_table(int identifier, bool generate_if_not_exist){
		auto itr = _table_h.find(identifier);
		if(itr == _table_h.end()){
			if(generate_if_not_exist){
				Table* table = new Table();
				_table_h[identifier] = table;
				return table;
			}
			return NULL;
		}
		return itr->second;
	}
	// 縦の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_clustering_vertical_sbr_ratio(double alpha){
		int pass_count = get_vertical_pass_count_with_id(CLUSTERING_TSSB_ID);
		int stop_count = get_vertical_stop_count_with_id(CLUSTERING_TSSB_ID);
		return (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
	}
	double compute_stop_probability_of_vertical_sbr_ratio(double alpha){
		vector<double> ratio_over_parents_reverse;	// 縦のSBRでの棒を折る比率を後ろから順に格納
		int num_parents = _depth_v;
		Node* target = this;
		for(int n = 0;n < num_parents + 1;n++){
			double expectation = target->compute_expectation_of_clustering_vertical_sbr_ratio(alpha);
			target = target->_parent;
			ratio_over_parents_reverse.push_back(expectation);
		}
		return compute_sbr_probability_given_params_reverse(ratio_over_parents_reverse);
	}
	// 縦の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_htssb_vertical_sbr_ratio(double alpha){
		vector<double> ratio_over_parents;	// 縦のSBRでの棒を折る比率をトップレベルから順に格納
		// トップレベルのノードから順に下りながら計算すると効率が良い
		int num_parents = _depth_v;
		for(int n = 0;n < num_parents;n++){
			Node* target = this;
			for(int step = 0;step < num_parents - n;step++){
				target = target->_parent;
				assert(target != NULL);
			}
			// 親のTSSBにおけるこのノードの期待値なので自分のメソッドを呼ぶ
			double expectation = _compute_expectation_of_htssb_vertical_sbr_ratio(alpha, target, ratio_over_parents);
			ratio_over_parents.push_back(expectation);
		}
		return _compute_expectation_of_htssb_vertical_sbr_ratio(alpha, this, ratio_over_parents);
	}
	double _compute_expectation_of_htssb_vertical_sbr_ratio(double alpha, Node* target, vector<double> &expectation_over_parents){
		int pass_count = get_vertical_pass_count_with_id(target->_identifier);
		int stop_count = get_vertical_stop_count_with_id(target->_identifier);
		if(target->_parent == NULL){
			return (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
		}
		double v_parent = expectation_over_parents.back();
		double sum_v_parents = std::accumulate(expectation_over_parents.begin(), expectation_over_parents.end(), 0.0);
		return (alpha * v_parent + stop_count) / (alpha * (1.0 - sum_v_parents) + stop_count + pass_count);
	}
	// 横の棒折り過程における、棒を折る比率を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_clustering_horizontal_sbr_ratio(double gamma){
		int pass_count = get_horizontal_pass_count_with_id(CLUSTERING_TSSB_ID);
		int stop_count = get_horizontal_stop_count_with_id(CLUSTERING_TSSB_ID);
		return (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
	}
	// 横の棒折り過程における、棒を折る比率を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_htssb_horizontal_sbr_ratio(double gamma){
		return 0.6;
	}
	// クラスタリングTSSBでこのノードに止まる確率。
	double compute_stop_probability(){
		return 0.6;
	}
	// クラスタリングTSSBに客を追加
	void add_customer_to_clustering_vertical_crp(double concentration){
		add_customer_to_htssb_vertical_crp(concentration, CLUSTERING_TSSB_ID, this);
	}
	// 遷移確率TSSBに客を追加
	void add_customer_to_htssb_vertical_crp(double concentration, int tssb_identifier, Node* node){
		Table* table = get_vertical_table(tssb_identifier, true);
		assert(table != NULL);
		bool new_table_generated = false;
		table->add_customer(concentration, new_table_generated);
		if(tssb_identifier != CLUSTERING_TSSB_ID && new_table_generated){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			assert(node != NULL);
			if(node->_parent != NULL){
				add_customer_to_htssb_vertical_crp(concentration, node->_parent->_identifier, node->_parent);		// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
			}
		}
		// 停止回数・通過回数を更新
		increment_vertical_stop_count(tssb_identifier);
		Node* parent = _parent;
		while(parent){
			parent->increment_vertical_pass_count(tssb_identifier);
			parent = parent->_parent;
		}
	}
	void increment_vertical_stop_count(int identifier){
		_stop_count_v[identifier] += 1;
	}
	void decrement_vertical_stop_count(int identifier){
		auto itr = _stop_count_v.find(identifier);
		assert(itr != _stop_count_v.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_stop_count_v.erase(itr);
		}
		// delete_node_if_needed(identifier);
	}
	void increment_vertical_pass_count(int identifier){
		_pass_count_v[identifier] += 1;
	}
	void decrement_vertical_pass_count(int identifier){
		// cout << "decrement_vertical_pass_count: " << _identifier << endl;
		auto itr = _pass_count_v.find(identifier);
		assert(itr != _pass_count_v.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_pass_count_v.erase(itr);
		}
		// delete_node_if_needed(identifier);
	}
	void add_customer_to_clustering_horizontal_crp(double concentration){
		add_customer_to_htssb_horizontal_crp(concentration, CLUSTERING_TSSB_ID, this);
	}
	void add_customer_to_htssb_horizontal_crp(double concentration, int tssb_identifier, Node* node){
		Table* table = get_htssb_horizontal_table(tssb_identifier, true);
		assert(table != NULL);
		bool new_table_generated = false;
		table->add_customer(concentration, new_table_generated);
		if(tssb_identifier != CLUSTERING_TSSB_ID && new_table_generated){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			assert(node != NULL);
			if(node->_parent != NULL){
				add_customer_to_htssb_horizontal_crp(concentration, node->_parent->_identifier, node->_parent);	// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
			}
		}
		// 停止回数・通過回数を更新
		Node* stopped_child = this;
		Node* parent = _parent;
		while(parent){
			for(int i = 0;i < parent->_children.size();i++){
				Node* child = parent->_children[i];
				if(child == stopped_child){
					child->increment_horizontal_stop_count(tssb_identifier);
					break;
				}
				child->increment_horizontal_pass_count(tssb_identifier);
			}
			stopped_child = parent;
			parent = parent->_parent;
		}
		stopped_child->increment_horizontal_stop_count(tssb_identifier);
	}
	void increment_horizontal_stop_count(int identifier){
		_stop_count_h[identifier] += 1;
	}
	void decrement_horizontal_stop_count(int identifier){
		// cout << "decrement_horizontal_stop_count: " << _identifier << endl;
		auto itr = _stop_count_h.find(identifier);
		assert(itr != _stop_count_h.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_stop_count_h.erase(itr);
		}
		// delete_node_if_needed(identifier);
	}
	void increment_horizontal_pass_count(int identifier){
		_pass_count_h[identifier] += 1;
	}
	void decrement_horizontal_pass_count(int identifier){
		auto itr = _pass_count_h.find(identifier);
		assert(itr != _pass_count_h.end());
		itr->second -= 1;
		assert(itr->second >= 0);
		if(itr->second == 0){
			_pass_count_h.erase(itr);
		}
		// delete_node_if_needed(identifier);
	}
	// 客を除去
	void remove_customer_from_clustering_vertical_crp(){
		remove_customer_from_htssb_vertical_crp(CLUSTERING_TSSB_ID, this);
	}
	void remove_customer_from_htssb_vertical_crp(int tssb_identifier, Node* node){
		// cout << "remove_customer_from_htssb_vertical_crp: " << tssb_identifier << ", " << node->_identifier << endl;
		Table* table = get_vertical_table(tssb_identifier, false);
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(tssb_identifier != CLUSTERING_TSSB_ID && empty_table_deleted && node->_parent != NULL){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			remove_customer_from_htssb_vertical_crp(tssb_identifier, node->_parent);	// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
		}
		// 停止回数・通過回数を更新
		decrement_vertical_stop_count(tssb_identifier);
		delete_node_if_needed(tssb_identifier);
		Node* parent = _parent;
		while(parent){
			parent->decrement_vertical_pass_count(tssb_identifier);
			parent->delete_node_if_needed(tssb_identifier);
			parent = parent->_parent;
		}
	}
	void remove_customer_from_clustering_horizontal_crp(){
		remove_customer_from_htssb_horizontal_crp(CLUSTERING_TSSB_ID, this);
	}
	void remove_customer_from_htssb_horizontal_crp(int tssb_identifier, Node* node){
		// cout << "remove_customer_from_htssb_horizontal_crp: " << _identifier << "," << node->_identifier << endl;
		Table* table = get_htssb_horizontal_table(tssb_identifier, true);
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(tssb_identifier != CLUSTERING_TSSB_ID && empty_table_deleted && node->_parent != NULL){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			remove_customer_from_htssb_horizontal_crp(tssb_identifier, node->_parent);	// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
		}
		// 停止回数・通過回数を更新
		Node* stopped_child = node;
		Node* parent = node->_parent;
		while(parent){
			bool found = false;
			for(int i = parent->_children.size() - 1;i >= 0;i--){	// 逆向きに辿らないと通過ノードが先に消えてしまう
				Node* child = parent->_children[i];
				// cout << "foreach: " << child->_identifier << endl;
				// cout << "foreach: " << child->_identifier << endl;

				if(child == stopped_child){
					found = true;
					child->decrement_horizontal_stop_count(tssb_identifier);
					child->delete_node_if_needed(tssb_identifier);
					continue;
				}
				if(found){
					child->decrement_horizontal_pass_count(tssb_identifier);
					child->delete_node_if_needed(tssb_identifier);
				}
			}
			stopped_child = parent;
			parent = parent->_parent;
		}
		stopped_child->decrement_horizontal_stop_count(tssb_identifier);
		stopped_child->delete_node_if_needed(tssb_identifier);
	}
	void delete_node_if_needed(int tssb_identifier){
		int pass_count_v = get_vertical_pass_count_with_id(tssb_identifier);
		int stop_count_v = get_vertical_stop_count_with_id(tssb_identifier);
		int pass_count_h = get_horizontal_pass_count_with_id(tssb_identifier);
		int stop_count_h = get_horizontal_stop_count_with_id(tssb_identifier);
		if(pass_count_v + stop_count_v + pass_count_h + stop_count_h == 0 && _parent != NULL){
			// cout << pass_count_v << "," << stop_count_v << "," << pass_count_h << "," << stop_count_h << endl;
			// cout << "requesting parent " << _parent->_identifier << ", me = " << _identifier << endl;
			_parent->delete_child_node(tssb_identifier, _identifier);
		}
	}
	void delete_child_node(int tssb_identifier,  int node_id){
		for(int i = 0;i < _children.size();i++){
			Node* target = _children[i];
			if(target->_identifier == node_id){
				_children.erase(_children.begin() + i);
				break;
			}
			target->delete_node_if_needed(tssb_identifier);
		}
		// if(_children.size() == 0 && _parent != NULL){
		// 	_parent->delete_child_node(tssb_identifier, _identifier);
		// }
	}
};
int Node::_auto_increment = CLUSTERING_TSSB_ID + 1;
#endif