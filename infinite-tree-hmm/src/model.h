#pragma once
#include "ithmm/ithmm.h"
#include "dictionary.h"

class Model{
public:
	iTHMM* _ithmm;
	Model(){
		// 日本語周り
		setlocale(LC_CTYPE, "ja_JP.UTF-8");
		std::ios_base::sync_with_stdio(false);
		std::locale default_loc("ja_JP.UTF-8");
		std::locale::global(default_loc);
		std::locale ctype_default(std::locale::classic(), default_loc, std::locale::ctype); //※
		std::wcout.imbue(ctype_default);
		std::wcin.imbue(ctype_default);

		_ithmm = new iTHMM();
	}
	~Model(){
		delete _ithmm;
	}
	double get_alpha(){
		return _ithmm->_alpha;
	}
	double get_gamma(){
		return _ithmm->_gamma;
	}
	double get_lambda_alpha(){
		return _ithmm->_lambda_alpha;
	}
	double get_lambda_gamma(){
		return _ithmm->_lambda_gamma;
	}
	double get_strength(){
		return _ithmm->_strength;
	}
	double get_tau0(){
		return _ithmm->_tau0;
	}
	double get_tau1(){
		return _ithmm->_tau1;
	}
	double get_metropolis_hastings_acceptance_rate(){
		return _ithmm->_num_mh_acceptance / (double)(_ithmm->_num_mh_acceptance + _ithmm->_num_mh_rejection);
	}
	boost::python::list get_all_tags(){
		boost::python::list tags;
		std::vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::string indices = "[" + node->_dump_indices() + "]";
			tags.append(indices);
		}
		return tags;
	}
	void set_alpha(double alpha){
		_ithmm->_alpha = alpha;
	}
	void set_gamma(double gamma){
		_ithmm->_gamma = gamma;
	}
	void set_lambda_alpha(double lambda_alpha){
		_ithmm->_lambda_alpha = lambda_alpha;
	}
	void set_lambda_gamma(double lambda_gamma){
		_ithmm->_lambda_gamma = lambda_gamma;
	}
	void set_strength(double strength){
		_ithmm->_strength = strength;
	}
	void set_tau0(double tau0){
		_ithmm->_tau0 = tau0;
	}
	void set_tau1(double tau1){
		_ithmm->_tau1 = tau1;
	}
	void set_depth_limit(int limit){
		_ithmm->set_depth_limit(limit);
	}
	void set_metropolis_hastings_enabled(bool enabled){
		_ithmm->_mh_enabled = enabled;
	}
	bool load(std::string filename){
		return _ithmm->load(filename);
	}
	bool save(std::string filename){
		return _ithmm->save(filename);
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void viterbi_decode_data(std::vector<Word*> &data, std::vector<Node*> &states, std::vector<Node*> &series, double** forward_table, double** decode_table){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < states.size();i++){
			Node* state = states[i];
			Node* state_on_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_on_bos != NULL);
			double Ps = state_on_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			forward_table[0][i] = Pword_given_s * Ps;
			decode_table[0][i] = 0;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int j = 0;j < states.size();j++){
				Node* state = states[j];
				forward_table[t][j] = 0;
				double max_value = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < states.size();i++){
					Node* prev_state = states[i];
					Node* state_on_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_on_prev_htssb != NULL);
					double Ps_given_prev = state_on_prev_htssb->_probability;
					double value = Ps_given_prev * forward_table[t - 1][i];
					if(value > max_value){
						max_value = value;
						forward_table[t][j] = value * Pword_given_s;
						decode_table[t][j] = i;
					}
				}
			}
		}
		// 後ろ向きに系列を復元
		std::vector<int> series_indices;
		int n = data.size() - 1;
		int k = 0;
		double max_value = 0;
		for(int i = 0;i < states.size();i++){
			if(forward_table[n][i] > max_value){
				k = i;
				max_value = forward_table[n][i];
			}
		}
		series_indices.push_back(k);
		for(int t = n - 1;t >= 0;t--){
			k = decode_table[t + 1][series_indices[n - t - 1]];
			series_indices.push_back(k);
		}
		std::reverse(series_indices.begin(), series_indices.end());
		// ノードをセット
		series.clear();
		for(int t = 0;t <= n;t++){
			int k = series_indices[t];
			series.push_back(states[k]);
		}
	}
	// データの対数尤度を計算
	// 前向きアルゴリズム
	double compute_Pdata(std::vector<Word*> &data, std::vector<Node*> &states, double** forward_table){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < states.size();i++){
			Node* state = states[i];
			Node* state_on_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_on_bos != NULL);
			double Ps = state_on_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			forward_table[0][i] = Pword_given_s * Ps;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int j = 0;j < states.size();j++){
				Node* state = states[j];
				forward_table[t][j] = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < states.size();i++){
					Node* prev_state = states[i];
					Node* state_on_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_on_prev_htssb != NULL);
					double Ps_given_prev = state_on_prev_htssb->_probability;
					forward_table[t][j] += Pword_given_s * Ps_given_prev * forward_table[t - 1][i];
				}
			}
		}
		int t = data.size() - 1;
		double Px = 0;
		for(int j = 0;j < states.size();j++){
			Px += forward_table[t][j];
		}
		return Px;
	}
	void show_assigned_words_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag, bool show_probability = true){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			std::string indices = node->_dump_indices();
			// linuxでバグるのでstringとwstring両方作る
			std::string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			std::cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m" << std::endl;
			std::wstring wtab = L"";
			for(int i = 0;i < node->_depth_v;i++){
				wtab += L"	";
			}
			std::wcout << wtab;
			for(const auto &elem: ranking){
				id word_id = elem.first;
				std::wstring &word = dict->_id_to_str[word_id];
				double p = elem.second;
				int count = node->_num_word_assignment[word_id];
				std::wcout << "\x1b[1m" << word << "\x1b[0m" << L" (" << count;
				if(show_probability){
					std::wcout << L";p=" << p;
				} 
				std::wcout << L") ";
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			std::wcout << std::endl;
			ranking.clear();
		}
	}
	void show_assigned_words_and_probability_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			std::string indices = node->_dump_indices();
			std::cout << "\x1b[32;1m" << "[" << indices << "]" << "\x1b[0m" << std::endl;
			for(const auto &elem: ranking){
				id word_id = elem.first;
				std::wstring &word = dict->_id_to_str[word_id];
				double p = elem.second;
				int count = node->_num_word_assignment[word_id];
				std::wcout << "\x1b[1m" << word << "\x1b[0m" << L"	" << count << L"	" << p << std::endl;
				n++;
				if(n > number_to_show_for_each_tag){
					break;
				}
			}
			std::wcout << std::endl;
			ranking.clear();
		}
	}
	void show_hpylm_for_each_tag(Dictionary* dict){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			std::string indices = node->_dump_indices();
			// linuxでバグるのでstringとwstring両方作る
			std::string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			std::cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m" << std::endl;
			std::wstring wtab = L"";
			for(int i = 0;i < node->_depth_v;i++){
				wtab += L"	";
			}
			std::wcout << wtab;
			for(const auto &table: node->_hpylm->_arrangement){
				id word_id = table.first;
				std::wstring &word = dict->_id_to_str[word_id];
				int num_tables = table.second.size();
				int num_customers = std::accumulate(table.second.begin(), table.second.end(), 0);
				std::wcout << "\x1b[1m" << word << "\x1b[0m" << L" (#t=" << num_tables << ";#c=" << num_customers << L") ";
			}
			std::wcout << std::endl;
		}
	}
	void show_sticks(){
		std::vector<Node*> nodes;
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::string indices = node->_dump_indices();
			c_printf("[*]%s\n", ((boost::format("[%s]") % indices.c_str())).str().c_str());
			_show_stick(node);
		}
	}
	void _show_stick(Node* node_on_structure){
		assert(node_on_structure != NULL);
		double p_eos = node_on_structure->compute_transition_probability_to_eos(_ithmm->_tau0, _ithmm->_tau1);
		_ithmm->update_stick_length_of_tssb(node_on_structure->_transition_tssb, 1.0 - p_eos, true);

		std::vector<Node*> nodes;
		node_on_structure->_transition_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::string indices = node->_dump_indices();
			std::string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			std::cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m " << node->_probability << std::endl;
		}
	}
};