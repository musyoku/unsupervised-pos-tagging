#ifndef _node_
#define _node_
#include <unordered_map>
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

class Node{
public:
	static _auto_increment;
	int _identifier;	// ノードID
	Node* _parent;		// 親ノード
	int _depth_v;		// 縦の深さ。 論文中の|s|に相当
	int _depth_h;		// 横の深さ。 論文中のkに相当
	// 各ノードのTSSBは自己同型になっている必要があるため、構造を共有する
	// カウントは各ノードのIDごとに管理
	unordered_map<int, int> _pass_count_v;	// 通過回数。 縦方向のCDP
	unordered_map<int, int> _stop_count_v;	// 停止回数。 縦方向のCDP
	unordered_map<int, int> _pass_count_h;	// 通過回数。 横方向のCDP
	unordered_map<int, int> _stop_count_h;	// 停止回数。 横方向のCDP
	unordered_map<int, Table*> _table_v;	// 客を管理するテーブル。 縦方向のCRP
	unordered_map<int, Table*> _table_h;	// 客を管理するテーブル。 横方向のCRP
	Node(Node* _parent){
		_identifier = _auto_increment;
		_auto_increment++;
		if(_parent != NULL){
		}
	}
	// 縦の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_sbr_param_v(){
		
	}
	// 縦の棒折り過程における、棒を折る比率をサンプリング。論文中のコインの表が出る確率に相当
	double sample_sbr_param_v(){
		
	}
	// 横の棒折り過程における、棒を折る比率を計算。論文中のコインの表が出る確率に相当
	double compute_expectation_of_sbr_param_h(){
		
	}
	// 横の棒折り過程における、棒を折る比率の期待値を計算。論文中のコインの表が出る確率に相当
	double sample_sbr_param_h(){
		
	}
	// TSSBでこのノードに止まる確率。
	double compute_stop_probability(){

	}
};
int Node::_auto_increment = 0;
#endif