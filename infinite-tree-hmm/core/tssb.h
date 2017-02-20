#ifndef _tssb_
#define _tssb_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <unordered_map>
#include <map>
#include <algorithm>
#include "cprintf.h"
#include "sampler.h"
#include "util.h"
using namespace std;
using namespace boost;

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
class TSSB;
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
	int _htssb_owner_id;
	Node* _parent;		// 親ノード
	int _depth_v;		// 縦の深さ。 論文中の|s|に相当
	int _depth_h;		// 横の深さ。 論文中のkに相当
	// 各ノードの遷移確率TSSBは自己同型になっている必要があるため、構造を共有する
	// カウントは各ノードのIDごとに管理
	int _pass_count_v;	// 通過回数。 縦方向のCDP
	int _stop_count_v;	// 停止回数。 縦方向のCDP
	int _pass_count_h;	// 通過回数。 横方向のCDP
	int _stop_count_h;	// 停止回数。 横方向のCDP
	Table* _table_v;	// 客を管理するテーブル。 縦方向のCRP
	Table* _table_h;	// 客を管理するテーブル。 横方向のCRP
	TSSB* _transition_tssb;		// 遷移確率を表すTSSBのルート
	Node* _transition_tssb_myself;		// 遷移確率を表すTSSBの自分と同じ位置のノード
	vector<Node*> _children;
	double _stick_length;					// 自分の棒の木全体に対する長さ
	double _children_stick_length;			// 自分の棒の子ノードに割り当てる長さ
	double _probability;					// このノードの確率
	// 計算時に使う配列
	double* _stop_probability_v_over_parent;
	double* _stop_ratio_v_over_parent;
	Node** _pointer_nodes_v;
	double* _stop_probability_h_over_parent;
	double* _stop_ratio_h_over_parent;
	Node(Node* parent){
		_identifier = _auto_increment;
		_auto_increment++;
		_parent = parent;
		init();
	}
	Node(Node* parent, int identifier){
		_identifier = identifier;
		_parent = parent;
		init();
	}
	void init(){
		_depth_v = (_parent != NULL) ? _parent->_depth_v + 1 : 0;
		_depth_h = (_parent != NULL) ? _parent->_children.size() : 0;
		_htssb_owner_id = 0;
		_stick_length = -1;
		_children_stick_length = -1;
		_pass_count_v = 0;
		_stop_count_v = 0;
		_pass_count_h = 0;
		_stop_count_h = 0;
		_probability = -1;
		_transition_tssb = NULL;
		_transition_tssb_myself = NULL;
		_stop_probability_v_over_parent = new double[_depth_v + 1];
		_stop_ratio_v_over_parent = new double[_depth_v + 1];
		_pointer_nodes_v = new Node*[_depth_v + 1];
		_stop_probability_h_over_parent = new double[_depth_h + 1];
		_stop_ratio_h_over_parent = new double[_depth_h + 1];
		_table_v = new Table();
		_table_h = new Table();
	}
	~Node(){
		delete _table_v;
		delete _table_h;
		delete[] _stop_probability_v_over_parent;
		delete[] _stop_ratio_v_over_parent;
		delete[] _pointer_nodes_v;
	}
	Node* generate_child(){
		Node* child = new Node(this);
		add_child(child);
		return child;
	}
	void add_child(Node* node){
		assert(node != NULL);
		_children.push_back(node);
	}
	int get_vertical_stop_count(){
		return _stop_count_v;
	}
	int get_vertical_pass_count(){
		return _pass_count_v;
	}
	int get_horizontal_stop_count(){
		return _stop_count_h;
	}
	int get_horizontal_pass_count(){
		return _pass_count_h;
	}
	Table* get_vertical_table(){
		return _table_v;
	}
	Table* get_horizontal_table(){
		return _table_h;
	}
	bool has_child(){
		return _children.size() != 0;
	}
	// 縦の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_clustering_vertical_sbr_ratio(double alpha){
		int pass_count = get_vertical_pass_count();
		int stop_count = get_vertical_stop_count();
		return (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
	}
	// 縦の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_htssb_vertical_sbr_ratio(double alpha, double lambda){
		double sbr_ratio = 0;
		int num_parents = _depth_v;
		// トップレベルのノードから順に下りながら計算する
		Node* parent = this;
		_pointer_nodes_v[_depth_v] = parent;
		for(int n = 0;n < num_parents;n++){
			parent = parent->_parent;
			_pointer_nodes_v[_depth_v - n - 1] = parent;
		}

		for(int n = 0;n < num_parents + 1;n++){
			// cout << "n = " << n << endl;
			Node* target = this;
			for(int step = 0;step < num_parents - n;step++){
				target = target->_parent;
				assert(target != NULL);
			}
			// cout << "target = " << target->_identifier << endl;
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			for(int m = 0;m < num_parents + 1;m++){
				// cout << "m = " << m << endl;
				Node* pointer = _pointer_nodes_v[m];
				// cout << "pointer = " << pointer->_identifier << endl;
				if(target->_depth_v == 0){	// 親ノードの場合
					int pass_count = pointer->get_vertical_pass_count();
					int stop_count = pointer->get_vertical_stop_count();
					double ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count);
					// cout << "ratio_v = " << ratio_v << endl;
					// pointer->set_vertical_stop_probability_with_id(stop_probability, target->_identifier);
					_stop_ratio_v_over_parent[m] = ratio_v;
				}else{
					int pass_count = pointer->get_vertical_pass_count();
					int stop_count = pointer->get_vertical_stop_count();
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					assert(target->_parent != NULL);
					double parent_probability_to_stop_pointer_node = _stop_probability_v_over_parent[m];
					// cout << "parent_probability_to_stop_pointer_node = " << parent_probability_to_stop_pointer_node << endl;
					double _alpha = alpha * pow(lambda, target->_depth_v);
					// cout << "_alpha = " << _alpha << endl;
					double ratio_v = (_alpha * parent_probability_to_stop_pointer_node + stop_count) / (_alpha * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
					// cout << "ratio_v = " << ratio_v << endl;
					_stop_ratio_v_over_parent[m] = ratio_v;
					// double stop_probability = rest_stick_length * ratio_v;
					// cout << "stop_probability = " << stop_probability << endl;
					// rest_stick_length *= 1.0 - ratio_v;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					sum_parent_stop_probability += parent_probability_to_stop_pointer_node;
					// pointer->set_vertical_stop_probability_with_id(stop_probability, target->_identifier);
					sbr_ratio = ratio_v;
					// cout << "sbr_ratio = " << sbr_ratio << endl;
				}
			}
			if(n < num_parents){
				double rest_stick_length = 1;
				for(int m = 0;m < num_parents + 1;m++){
					double ratio_v = _stop_ratio_v_over_parent[m];
					double stop_probability = rest_stick_length * ratio_v;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_v;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					_stop_probability_v_over_parent[m] = stop_probability;
				}
			}
		}
		return sbr_ratio;
	}
	// 横の棒折り過程における、棒を折る比率を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_clustering_horizontal_sbr_ratio(double gamma){
		int pass_count = get_horizontal_pass_count();
		int stop_count = get_horizontal_stop_count();
		return (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
	}
	// 横の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double compute_expectation_of_htssb_horizontal_sbr_ratio(double gamma){
		if(_depth_v == 0){
			return 1;
		}
		double sbr_ratio = 0;
		int num_parents = _depth_v;

		for(int n = 0;n < num_parents;n++){
			cout << "n = " << n << endl;
			Node* target = this;
			for(int step = 0;step < num_parents - n;step++){
				target = target->_parent;
				assert(target != NULL);
			}
			// cout << "target = " << target->_identifier << endl;
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			for(int m = 0;m < _depth_h + 1;m++){
				Node* pointer = _parent->_children[m];
				// cout << "m = " << m << endl;
				// cout << "size = " << _parent->_children.size() << endl;
				// cout << "pointer = " << pointer->_identifier << endl;
				if(target->_depth_v == 0){	// 親ノードの場合
					int pass_count = pointer->get_horizontal_pass_count();
					int stop_count = pointer->get_horizontal_stop_count();
					double ratio_h = (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
					// cout << "ratio_h = " << ratio_h << endl;
					// pointer->set_horizontal_stop_probability_with_id(stop_probability, );
					_stop_ratio_h_over_parent[m] = ratio_h;
				}else{
					int pass_count = pointer->get_horizontal_pass_count();
					int stop_count = pointer->get_horizontal_stop_count();
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					assert(target->_parent != NULL);
					double parent_probability_to_stop_pointer_node = _stop_probability_h_over_parent[m];
					// cout << "parent_probability_to_stop_pointer_node = " << parent_probability_to_stop_pointer_node << endl;
					double ratio_h = (gamma * parent_probability_to_stop_pointer_node + stop_count) / (gamma * (1.0 - sum_parent_stop_probability) + stop_count + pass_count);
					// cout << "ratio_h = " << ratio_h << endl;
					_stop_ratio_h_over_parent[m] = ratio_h;
					sum_parent_stop_probability += parent_probability_to_stop_pointer_node;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					sbr_ratio = ratio_h;
					// cout << "sbr_ratio = " << sbr_ratio << endl;
				}
			}
			if(n < num_parents - 1){
				double rest_stick_length = 1;
				for(int m = 0;m < num_parents + 1;m++){
					double ratio_h = _stop_ratio_h_over_parent[m];
					double stop_probability = rest_stick_length * ratio_h;
					// cout << "stop_probability = " << stop_probability << endl;
					rest_stick_length *= 1.0 - ratio_h;
					// cout << "rest_stick_length = " << rest_stick_length << endl;
					_stop_probability_h_over_parent[m] = stop_probability;
				}
			}
		}
		return sbr_ratio;
	}
	// クラスタリングTSSBでこのノードに止まる確率。
	double compute_stop_probability(){
		return 0.6;
	}
	// 遷移確率TSSBに客を追加
	void add_customer_to_vertical_crp(double concentration, bool &new_table_generated){
		Table* table = get_vertical_table();
		Node* parent = _parent;
		assert(table != NULL);
		table->add_customer(concentration, new_table_generated);
		// 停止回数・通過回数を更新
		increment_vertical_stop_count();
		while(parent){
			parent->increment_vertical_pass_count();
			parent = parent->_parent;
		}
	}
	void increment_vertical_stop_count(){
		_stop_count_v += 1;
	}
	void decrement_vertical_stop_count(){
		_stop_count_v -= 1;
		assert(_stop_count_v >= 0);
	}
	void increment_vertical_pass_count(){
		_pass_count_v += 1;
	}
	void decrement_vertical_pass_count(){
		_pass_count_v -= 1;
		assert(_pass_count_v >= 0);
	}
	void add_customer_to_horizontal_crp(double concentration, bool &new_table_generated){
		Table* table = get_horizontal_table();
		assert(table != NULL);
		table->add_customer(concentration, new_table_generated);
		// 停止回数・通過回数を更新
		Node* stopped_child = this;
		Node* parent = _parent;
		while(parent){
			for(int i = 0;i < parent->_children.size();i++){
				Node* child = parent->_children[i];
				if(child == stopped_child){
					child->increment_horizontal_stop_count();
					break;
				}
				child->increment_horizontal_pass_count();
			}
			stopped_child = parent;
			parent = parent->_parent;
		}
		stopped_child->increment_horizontal_stop_count();
	}
	void increment_horizontal_stop_count(){
		_stop_count_h += 1;
	}
	void decrement_horizontal_stop_count(){
		_stop_count_h -= 1;
		assert(_stop_count_h >= 0);
	}
	void increment_horizontal_pass_count(){
		_pass_count_h += 1;
	}
	void decrement_horizontal_pass_count(){
		_pass_count_h -= 1;
		assert(_pass_count_h >= 0);
	}
	// 客を除去
	void remove_customer_from_vertical_crp(bool remove_proxy_customer = false){
		remove_customer_from_htssb_vertical_crp(this, remove_proxy_customer);
	}
	void remove_customer_from_htssb_vertical_crp(Node* node, bool remove_proxy_customer = false){
		// cout << "remove_customer_from_htssb_vertical_crp: " << tssb_identifier << ", " << node->_identifier << endl;
		Table* table = get_vertical_table();
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(remove_proxy_customer && empty_table_deleted){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			assert(node != NULL);
			if(node->_parent != NULL){
				remove_customer_from_htssb_vertical_crp(node->_parent, true);	// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
			}
		}
		// 停止回数・通過回数を更新
		decrement_vertical_stop_count();
		delete_node_if_needed();
		Node* parent = _parent;
		while(parent){
			parent->decrement_vertical_pass_count();
			parent->delete_node_if_needed();
			parent = parent->_parent;
		}
	}
	void remove_customer_from_horizontal_crp(bool remove_proxy_customer = false){
		remove_customer_from_horizontal_crp(this, remove_proxy_customer);
	}
	void remove_customer_from_horizontal_crp(Node* node, bool remove_proxy_customer = false){
		// cout << "remove_customer_from_horizontal_crp: " << _identifier << "," << node->_identifier << endl;
		Table* table = get_horizontal_table();
		assert(table != NULL);
		bool empty_table_deleted = false;
		table->remove_customer(empty_table_deleted);
		if(remove_proxy_customer && empty_table_deleted){	// 新しいテーブルが作られたら親のTSSBのこのノードに代理客を追加
			assert(node != NULL);
			if(node->_parent != NULL){
				remove_customer_from_horizontal_crp(node->_parent, true);	// 親のTSSBにおけるこのノードに客を追加するので自分のメソッドを呼ぶ
			}
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
					child->decrement_horizontal_stop_count();
					child->delete_node_if_needed();
					continue;
				}
				if(found){
					child->decrement_horizontal_pass_count();
					child->delete_node_if_needed();
				}
			}
			stopped_child = parent;
			parent = parent->_parent;
		}
		stopped_child->decrement_horizontal_stop_count();
		stopped_child->delete_node_if_needed();
	}
	void delete_node_if_needed(){
		int pass_count_v = get_vertical_pass_count();
		int stop_count_v = get_vertical_stop_count();
		int pass_count_h = get_horizontal_pass_count();
		int stop_count_h = get_horizontal_stop_count();
		if(pass_count_v + stop_count_v + pass_count_h + stop_count_h == 0 && _parent != NULL){
			// cout << pass_count_v << "," << stop_count_v << "," << pass_count_h << "," << stop_count_h << endl;
			// cout << "requesting parent " << _parent->_identifier << ", me = " << _identifier << endl;
			_parent->delete_child_node(_identifier);
		}
	}
	void delete_child_node(int node_id){
		for(int i = 0;i < _children.size();i++){
			Node* target = _children[i];
			if(target->_identifier == node_id){
				_children.erase(_children.begin() + i);
				break;
			}
			target->delete_node_if_needed();
		}
	}
};
int Node::_auto_increment = 1;

class TSSB{
public:
	Node* _root;
	double _alpha;
	double _gamma;
	double _lambda;
	int _owner_id;
	TSSB(double alpha, double gamma, double lambda){
		_root = new Node(NULL);
		_root->_stick_length = 1;
		_owner_id = 0;
		_alpha = alpha;
		_gamma = gamma;
		_lambda = lambda;
	}
	TSSB(Node* root, double alpha, double gamma, double lambda){
		_root = root;
		_owner_id = 0;
		_alpha = alpha;
		_gamma = gamma;
		_lambda = lambda;
	}
	TSSB* copy(int owner_id){
		Node* root = new Node(NULL, _root->_identifier);
		root->_htssb_owner_id = owner_id;
		copy_children(_root, root, owner_id);
		TSSB* target = new TSSB(root, _alpha, _gamma, _lambda);
		target->_owner_id = owner_id;
		return target;
	}
	void copy_children(Node* source, Node* target, int owner_id){
		for(const auto source_child: source->_children){
			Node* child = new Node(target, source_child->_identifier);
			child->_htssb_owner_id = owner_id;
			target->add_child(child);
			copy_children(source_child, child, owner_id);
		}
	}
	void update_stick_length(){
		double ratio_v = _root->compute_expectation_of_clustering_vertical_sbr_ratio(_gamma);
		double sum_probability = ratio_v;
		_root->_stick_length = 1;
		_root->_children_stick_length = 1.0 - ratio_v;
		_root->_probability = ratio_v;
		_update_stick_length(sum_probability, _root);
		cout << "sum_probability: " << sum_probability << endl;
	}
	void _update_stick_length(double &sum_probability, Node* node){
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
	double compute_expectation_of_htssb_vertical_sbr_ratio(Node* node){
		return node->compute_expectation_of_htssb_vertical_sbr_ratio(_alpha, _lambda);
	}
	double compute_expectation_of_htssb_horizontal_sbr_ratio(Node* node){
		return node->compute_expectation_of_htssb_horizontal_sbr_ratio(_gamma);
	}
	void enumerate_nodes_from_left_to_right(vector<Node*> &nodes){
		_enumerate_nodes_from_left_to_right(_root, nodes);
	}
	void _enumerate_nodes_from_left_to_right(Node* node, vector<Node*> &nodes){
		nodes.push_back(node);
		for(const auto child: node->_children){
			_enumerate_nodes_from_left_to_right(child, nodes);
		}
	}
	Node* find_node_with_id(int identifier){
		if(_root->_identifier == identifier){
			return _root;
		}
		return _find_node_with_id(identifier, _root);
	}
	Node* _find_node_with_id(int identifier,  Node* node){
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
	int get_num_nodes(){
		return _get_num_children(_root) + 1;
	}
	int _get_num_children(Node* node){
		int sum = node->_children.size();
		for(const auto child: node->_children){
			sum += _get_num_children(child);
		}
		return sum;
	}
	int get_max_depth(){
		return _get_max_depth(_root);
	}
	int _get_max_depth(Node* node){
		int max_depth = node->_depth_v;
		for(const auto child: node->_children){
			int depth = _get_max_depth(child);
			if(depth > max_depth){
				max_depth = depth;
			}
		}
		return max_depth;
	}
	void dump(){
		_dump(_root);
	}
	void _dump(Node* node){
		string tab = "";
		for(int i = 0;i < node->_depth_v;i++){
			tab += "	";
		}
		cout << tab;
		int pass_count_v = node->get_vertical_pass_count();
		int stop_count_v = node->get_vertical_stop_count();
		int pass_count_h = node->get_horizontal_pass_count();
		int stop_count_h = node->get_horizontal_stop_count();
		cout << (boost::format("%d [vp:%d,vs:%d,hp:%d,hs:%d][len:%f,self:%f,ch:%f,p:%f][ow:%d]") % node->_identifier % pass_count_v % stop_count_v % pass_count_h % stop_count_h % node->_stick_length % (node->_stick_length - node->_children_stick_length) % node->_children_stick_length % node->_probability % node->_htssb_owner_id).str() << endl;
		for(int i = 0;i < node->_children.size();i++){
			Node* child = node->_children[i];
			_dump(child);
		}
	}
};
#endif