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
#include "c_printf.h"
#include "sampler.h"
#include "node.h"
#include "const.h"
#include "vocab.h"

class HPYLM{
private:
	int _hpylm_depth;				// 最大の深さ
	friend class boost::serialization::access;
	template <class Archive>
	// モデルの保存
	void serialize(Archive& archive, unsigned int version)
	{
		static_cast<void>(version); // No use
		archive & _root;
		archive & _hpylm_depth;
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
	int _max_depth;				// VPYLMへ拡張時に使う
	double _g0;					// ゼログラム確率

	// 深さmのノードに関するパラメータ
	vector<double> _d_m;		// Pitman-Yor過程のディスカウント係数
	vector<double> _theta_m;	// Pitman-Yor過程の集中度

	// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
	// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
	vector<double> _a_m;		// ベータ分布のパラメータ	dの推定用
	vector<double> _b_m;		// ベータ分布のパラメータ	dの推定用
	vector<double> _alpha_m;	// ガンマ分布のパラメータ	θの推定用
	vector<double> _beta_m;		// ガンマ分布のパラメータ	θの推定用

	HPYLM(int ngram = 2){
		// 深さは0から始まることに注意
		// 2-gramなら最大深さは1. root(0) -> 2-gram(1)
		// 3-gramなら最大深さは2. root(0) -> 2-gram(1) -> 3-gram(2)
		_hpylm_depth = ngram - 1;
		_max_depth = -1;

		_root = new Node();
		_root->_depth = 0;		// ルートは深さ0

		for(int n = 0;n < ngram;n++){
			_d_m.push_back(PYLM_INITIAL_D);	
			_theta_m.push_back(PYLM_INITIAL_THETA);
			_a_m.push_back(PYLM_INITIAL_A);	
			_b_m.push_back(PYLM_INITIAL_B);	
			_alpha_m.push_back(PYLM_INITIAL_ALPHA);
			_beta_m.push_back(PYLM_INITIAL_BETA);
		}
	}
	int ngram(){
		return _hpylm_depth + 1;
	}
	void set_g0(double g0){
		_g0 = g0;
	}
	// 単語列のindex番目の単語をモデルに追加
	bool add_customer_at_timestep(vector<id> &token_ids, int token_t_index){
		Node* node = find_node_by_tracing_back_context(token_ids, token_t_index, _hpylm_depth, true);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "客を追加できません. ノードが見つかりません.");
			exit(1);
		}
		id token_t = token_ids[token_t_index];
		int added_to_table_k;
		node->add_customer(token_t, _g0, _d_m, _theta_m, true, added_to_table_k);
		return true;
	}
	bool remove_customer_at_timestep(vector<id> &token_ids, int token_t_index){
		Node* node = find_node_by_tracing_back_context(token_ids, token_t_index, _hpylm_depth, false);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "客を除去できません. ノードが見つかりません.");
			exit(1);
		}
		id token_t = token_ids[token_t_index];
		int removed_from_table_k;
		node->remove_customer(token_t, true, removed_from_table_k);
		// 客が一人もいなくなったらノードを削除する
		if(node->need_to_remove_from_parent()){
			node->remove_from_parent();
		}
		return true;
	}
	// token列の位置tからorderだけ遡る
	// token_ids:        [0, 1, 2, 3, 4, 5]
	// token_t_index:4          ^     ^
	// depth_t: 2               |<- <-|
	Node* find_node_by_tracing_back_context(vector<id> &token_ids, int token_t_index, int depth_t, bool generate_node_if_needed = false, bool return_middle_node = false){
		if(token_t_index - depth_t < 0){
			return NULL;
		}
		Node* node = _root;
		for(int depth = 1;depth <= depth_t;depth++){
			id context_token_id = token_ids[token_t_index - depth];
			Node* child = node->find_child_node(context_token_id, generate_node_if_needed);
			if(child == NULL){
				if(return_middle_node){
					return node;
				}
				return NULL;
			}
			node = child;
		}
		return node;
	}
	double Pw_h(vector<id> &token_ids, vector<id> context_token_ids){
		double p = 1;
		for(int n = 0;n < token_ids.size();n++){
			p *= Pw_h(token_ids[n], context_token_ids);
			context_token_ids.push_back(token_ids[n]);
		}
		return p;
	}
	double Pw_h(id token_id, vector<id> &context_token_ids){
		// HPYLMでは深さは固定
		if(context_token_ids.size() < _hpylm_depth){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. context_token_ids.size() < _hpylm_depth");
			exit(1);
		}
		Node* node = find_node_by_tracing_back_context(context_token_ids, context_token_ids.size(), _hpylm_depth, false, true);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. node == NULL");
			exit(1);
		}
		return node->Pw(token_id, _g0, _d_m, _theta_m);
	}
	double Pw(id token_id){
		return _root->Pw(token_id, _g0, _d_m, _theta_m);
	}
	double Pw(vector<id> &token_ids){
		if(token_ids.size() < _hpylm_depth + 1){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. token_ids.size() < _hpylm_depth");
			exit(1);
		}
		double mul_Pw_h = 1;
		vector<id> context_token_ids(token_ids.begin(), token_ids.begin() + _hpylm_depth);
		for(int depth = _hpylm_depth;depth < token_ids.size();depth++){
			id token_id = token_ids[depth];
			mul_Pw_h *= Pw_h(token_id, context_token_ids);;
			context_token_ids.push_back(token_id);
		}
		return mul_Pw_h;
	}
	double log_Pw(vector<id> &token_ids){
		if(token_ids.size() < _hpylm_depth + 1){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. token_ids.size() < _hpylm_depth");
			exit(1);
		}
		double sum_Pw_h = 0;
		vector<id> context_token_ids(token_ids.begin(), token_ids.begin() + _hpylm_depth);
		for(int depth = _hpylm_depth;depth < token_ids.size();depth++){
			id token_id = token_ids[depth];
			sum_Pw_h += log(Pw_h(token_id, context_token_ids) + 1e-10);
			context_token_ids.push_back(token_id);
		}
		return sum_Pw_h;
	}
	double log2_Pw(vector<id> &token_ids){
		if(token_ids.size() < _hpylm_depth + 1){
			c_printf("[r]%s [*]%s\n", "エラー:", "単語確率を計算できません. token_ids.size() < _hpylm_depth");
			exit(1);
		}
		double sum_Pw_h = 0;
		vector<id> context_token_ids(token_ids.begin(), token_ids.begin() + _hpylm_depth);
		for(int depth = _hpylm_depth;depth < token_ids.size();depth++){
			id token_id = token_ids[depth];
			sum_Pw_h += log2(Pw_h(token_id, context_token_ids) + 1e-10);
			context_token_ids.push_back(token_id);
		}
		return sum_Pw_h;
	}
	id sample_next_token(vector<id> &context_token_ids, id eos_id){
		Node* node = find_node_by_tracing_back_context(context_token_ids, context_token_ids.size(), _hpylm_depth, false, true);
		if(node == NULL){
			c_printf("[r]%s [*]%s\n", "エラー:", "トークンを生成できません. ノードが見つかりません.");
			exit(1);
		}
		vector<id> token_ids;
		vector<double> probs;
		double sum_probs = 0;
		for(const auto &elem: node->_arrangement){
			id token_id = elem.first;
			double prob = Pw_h(token_id, context_token_ids);
			if(prob > 0){
				token_ids.push_back(token_id);
				probs.push_back(prob);
				sum_probs += prob;
			}
		}
		if(token_ids.size() == 0){
			return eos_id;
		}
		if(sum_probs == 0){
			return eos_id;
		}
		double ratio = 1.0 / sum_probs;
		double r = Sampler::uniform(0, 1);
		sum_probs = 0;
		id sampled_token_id = token_ids.back();
		for(int i = 0;i < token_ids.size();i++){
			sum_probs += probs[i] * ratio;
			if(sum_probs > r){
				return token_ids[i];
			}
		}
		return sampled_token_id;
	}
	void init_hyperparameters_at_depth_if_needed(int depth){
		if(depth >= _d_m.size()){
			while(_d_m.size() <= depth){
				_d_m.push_back(PYLM_INITIAL_D);
			}
		}
		if(depth >= _theta_m.size()){
			while(_theta_m.size() <= depth){
				_theta_m.push_back(PYLM_INITIAL_THETA);
			}
		}
		if(depth >= _a_m.size()){
			while(_a_m.size() <= depth){
				_a_m.push_back(PYLM_INITIAL_A);
			}
		}
		if(depth >= _b_m.size()){
			while(_b_m.size() <= depth){
				_b_m.push_back(PYLM_INITIAL_B);
			}
		}
		if(depth >= _alpha_m.size()){
			while(_alpha_m.size() <= depth){
				_alpha_m.push_back(PYLM_INITIAL_ALPHA);
			}
		}
		if(depth >= _beta_m.size()){
			while(_beta_m.size() <= depth){
				_beta_m.push_back(PYLM_INITIAL_BETA);
			}
		}
	}
	// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
	// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
	void sum_auxiliary_variables_recursively(Node* node, vector<double> &sum_log_x_u_m, vector<double> &sum_y_ui_m, vector<double> &sum_1_y_ui_m, vector<double> &sum_1_z_uwkj_m, int &bottom){
		for(const auto &elem: node->_children){
			Node* child = elem.second;
			int depth = child->_depth;

			if(depth > bottom){
				bottom = depth;
			}
			init_hyperparameters_at_depth_if_needed(depth);

			double d = _d_m[depth];
			double theta = _theta_m[depth];
			sum_log_x_u_m[depth] += child->auxiliary_log_x_u(theta);	// log(x_u)
			sum_y_ui_m[depth] += child->auxiliary_y_ui(d, theta);		// y_ui
			sum_1_y_ui_m[depth] += child->auxiliary_1_y_ui(d, theta);	// 1 - y_ui
			sum_1_z_uwkj_m[depth] += child->auxiliary_1_z_uwkj(d);		// 1 - z_uwkj

			sum_auxiliary_variables_recursively(child, sum_log_x_u_m, sum_y_ui_m, sum_1_y_ui_m, sum_1_z_uwkj_m, bottom);
		}
	}
	// dとθの推定
	void sample_hyperparams(){
		int max_depth = _d_m.size() - 1;

		// 親ノードの深さが0であることに注意
		vector<double> sum_log_x_u_m(max_depth + 1, 0.0);
		vector<double> sum_y_ui_m(max_depth + 1, 0.0);
		vector<double> sum_1_y_ui_m(max_depth + 1, 0.0);
		vector<double> sum_1_z_uwkj_m(max_depth + 1, 0.0);

		// _root
		sum_log_x_u_m[0] = _root->auxiliary_log_x_u(_theta_m[0]);			// log(x_u)
		sum_y_ui_m[0] = _root->auxiliary_y_ui(_d_m[0], _theta_m[0]);		// y_ui
		sum_1_y_ui_m[0] = _root->auxiliary_1_y_ui(_d_m[0], _theta_m[0]);	// 1 - y_ui
		sum_1_z_uwkj_m[0] = _root->auxiliary_1_z_uwkj(_d_m[0]);				// 1 - z_uwkj

		// それ以外
		_max_depth = 0;
		// _max_depthは以下を実行すると更新される
		// HPYLMでは無意味だがVPYLMで最大深さを求める時に使う
		sum_auxiliary_variables_recursively(_root, sum_log_x_u_m, sum_y_ui_m, sum_1_y_ui_m, sum_1_z_uwkj_m, _max_depth);
		init_hyperparameters_at_depth_if_needed(_max_depth);

		for(int u = 0;u <= _max_depth;u++){
			_d_m[u] = Sampler::beta(_a_m[u] + sum_1_y_ui_m[u], _b_m[u] + sum_1_z_uwkj_m[u]);
			_theta_m[u] = Sampler::gamma(_alpha_m[u] + sum_y_ui_m[u], _beta_m[u] - sum_log_x_u_m[u]);
		}
		// 不要な深さのハイパーパラメータを削除
		int num_remove = _d_m.size() - _max_depth - 1;
		for(int n = 0;n < num_remove;n++){
			_d_m.pop_back();
			_theta_m.pop_back();
			_a_m.pop_back();
			_b_m.pop_back();
			_alpha_m.pop_back();
			_beta_m.pop_back();
		}
	}
	int get_max_depth(bool use_cache = true){
		if(use_cache && _max_depth != -1){
			return _max_depth;
		}
		_max_depth = 0;
		update_max_depth_recursively(_root);
		return _max_depth;
	}
	void update_max_depth_recursively(Node* node){
		for(const auto &elem: node->_children){
			Node* child = elem.second;
			int depth = child->_depth;
			if(depth > _max_depth){
				_max_depth = depth;
			}
			update_max_depth_recursively(child);
		}
	}
	int get_num_nodes(){
		return _root->get_num_nodes() + 1;
	}
	int get_num_customers(){
		return _root->get_num_customers();
	}
	int get_num_tables(){
		return _root->get_num_tables();
	}
	int get_sum_stop_counts(){
		return _root->sum_stop_counts();
	}
	int get_sum_pass_counts(){
		return _root->sum_pass_counts();
	}
	void set_active_tokens(unordered_map<id, bool> &flags){
		_root->set_active_tokens(flags);
	}
	void count_tokens_of_each_depth(unordered_map<int, int> &map){
		_root->count_tokens_of_each_depth(map);
	}
	void enumerate_phrases_at_depth(int depth, vector<vector<id>> &phrases){
		int max_depth = get_max_depth();
		if(depth > max_depth){
			c_printf("[r]%s [*]%s\n", "エラー:", "指定の深さにフレーズが存在しません.");
			return;
		}
		// 指定の深さのノードを探索
		vector<Node*> nodes;
		_root->enumerate_nodes_at_depth(depth, nodes);
		for(int i = 0;i < nodes.size();i++){
			Node* node = nodes[i];
			vector<id> phrase;
			while(node->_parent){
				phrase.push_back(node->_token_id);
				node = node->_parent;
			}
			phrases.push_back(phrase);
		}
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