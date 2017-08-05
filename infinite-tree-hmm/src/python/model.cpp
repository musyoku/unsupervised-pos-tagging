#include "model.h"

namespace ithmm {
	Model::Model(){
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
	Model::~Model(){
		delete _ithmm;
	}
	double Model::get_alpha(){
		return _ithmm->_alpha;
	}
	double Model::get_gamma(){
		return _ithmm->_gamma;
	}
	double Model::get_lambda_alpha(){
		return _ithmm->_lambda_alpha;
	}
	double Model::get_lambda_gamma(){
		return _ithmm->_lambda_gamma;
	}
	double Model::get_strength(){
		return _ithmm->_strength;
	}
	double Model::get_tau0(){
		return _ithmm->_tau0;
	}
	double Model::get_tau1(){
		return _ithmm->_tau1;
	}
	double Model::get_metropolis_hastings_acceptance_rate(){
		return _ithmm->_num_mh_acceptance / (double)(_ithmm->_num_mh_acceptance + _ithmm->_num_mh_rejection);
	}
	boost::python::list Model::python_get_all_states(){
		boost::python::list tags;
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
		for(const auto &node: nodes){
			std::string indices = "[" + node->_dump_indices() + "]";
			tags.append(indices);
		}
		return tags;
	}
	void Model::set_alpha(double alpha){
		_ithmm->_alpha = alpha;
	}
	void Model::set_gamma(double gamma){
		_ithmm->_gamma = gamma;
	}
	void Model::set_lambda_alpha(double lambda_alpha){
		_ithmm->_lambda_alpha = lambda_alpha;
	}
	void Model::set_lambda_gamma(double lambda_gamma){
		_ithmm->_lambda_gamma = lambda_gamma;
	}
	void Model::set_strength(double strength){
		_ithmm->_strength = strength;
	}
	void Model::set_tau0(double tau0){
		_ithmm->_tau0 = tau0;
	}
	void Model::set_tau1(double tau1){
		_ithmm->_tau1 = tau1;
	}
	void Model::set_depth_limit(int limit){
		_ithmm->set_depth_limit(limit);
	}
	void Model::set_metropolis_hastings_enabled(bool enabled){
		_ithmm->_mh_enabled = enabled;
	}
	bool Model::load(std::string filename){
		return _ithmm->load(filename);
	}
	bool Model::save(std::string filename){
		return _ithmm->save(filename);
	}
	// 存在する全ての状態を集める
	void Model::enumerate_all_states(std::vector<Node*> &nodes){
		assert(nodes.size() == 0);
		_ithmm->_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
	}
	// 全ての棒の長さを計算しておく
	void Model::precompute_all_stick_lengths(std::vector<Node*> &all_states){
		assert(all_states.size() > 0);
		for(auto node: all_states){
			double p_eos_given_s = node->compute_transition_probability_to_eos(_ithmm->_tau0, _ithmm->_tau1);
			double total_stick_length = 1.0 - p_eos_given_s;	// <eos>以外に遷移する確率をTSSBで分配する
			_ithmm->update_stick_length_of_tssb(node->_transition_tssb, total_stick_length, true);
		}
		_ithmm->update_stick_length_of_tssb(_ithmm->_bos_tssb, 1.0, false);
	}
	boost::python::list Model::python_viterbi_decode(boost::python::list py_word_ids){
		// あらかじめ全HTSSBの棒の長さを計算しておく
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
		precompute_all_stick_lengths(nodes);
		// デコード用のテーブルを確保
		int num_words = boost::python::len(py_word_ids);
		double** forward_table = new double*[num_words];
		double** decode_table = new double*[num_words];
		for(int i = 0;i < num_words;i++){
			forward_table[i] = new double[nodes.size()];
			decode_table[i] = new double[nodes.size()];
		}
		// Python側から渡された単語IDリストを変換
		std::vector<Word*> words;
		for(int i = 0;i < num_words;i++){
			Word* word = new Word();
			word->_id = boost::python::extract<id>(py_word_ids[i]);
			word->_state = NULL;
			words.push_back(word);
		}
		// ビタビアルゴリズム
		std::vector<Node*> sampled_state_sequence;
		_viterbi_decode(words, nodes, sampled_state_sequence, forward_table, decode_table);
		// 結果を返す
		boost::python::list result;
		for(int i = 0;i < words.size();i++){
			std::wstring tag = L"[" + sampled_state_sequence[i]->_wdump_indices() + L"]";
			result.append(tag);
		}
		for(int i = 0;i < num_words;i++){
			delete[] forward_table[i];
			delete[] decode_table[i];
		}
		delete[] forward_table;
		delete[] decode_table;
		for(int i = 0;i < words.size();i++){
			delete words[i];
		}
		return result;
	}
	// 状態系列の復号
	// ビタビアルゴリズム
	void Model::_viterbi_decode(std::vector<Word*> &data, std::vector<Node*> &all_states, std::vector<Node*> &sampled_state_sequence, double** forward_table, double** decode_table){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < all_states.size();i++){
			Node* state = all_states[i];
			Node* state_in_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_in_bos != NULL);
			double Ps = state_in_bos->_probability;
			double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
			assert(Ps > 0);
			assert(Pword_given_s > 0);
			forward_table[0][i] = Pword_given_s * Ps;
			decode_table[0][i] = 0;
		}
		for(int t = 1;t < data.size();t++){
			Word* word = data[t];
			for(int j = 0;j < all_states.size();j++){
				Node* state = all_states[j];
				forward_table[t][j] = 0;
				double max_value = 0;
				double Pword_given_s = _ithmm->compute_Pw_given_s(word->_id, state);
				for(int i = 0;i < all_states.size();i++){
					Node* prev_state = all_states[i];
					Node* state_in_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_in_prev_htssb != NULL);
					double Ps_given_prev = state_in_prev_htssb->_probability;
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
		for(int i = 0;i < all_states.size();i++){
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
		sampled_state_sequence.clear();
		for(int t = 0;t <= n;t++){
			int k = series_indices[t];
			sampled_state_sequence.push_back(all_states[k]);
		}
	}
	// データの対数尤度を計算
	// 前向きアルゴリズム
	double Model::compute_Pdata(std::vector<Word*> &data, std::vector<Node*> &states, double** forward_table){
		// 初期化
		Word* word = data[0];
		for(int i = 0;i < states.size();i++){
			Node* state = states[i];
			Node* state_in_bos = _ithmm->_bos_tssb->find_node_by_tracing_horizontal_indices(state);
			assert(state_in_bos != NULL);
			double Ps = state_in_bos->_probability;
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
					Node* state_in_prev_htssb = prev_state->_transition_tssb->find_node_by_tracing_horizontal_indices(state);
					assert(state_in_prev_htssb != NULL);
					double Ps_given_prev = state_in_prev_htssb->_probability;
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
	void Model::show_assigned_words_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag, bool show_probability){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
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
	void Model::show_assigned_words_and_probability_for_each_tag(Dictionary* dict, int number_to_show_for_each_tag){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
		std::wcout << "word      	count	probability" << std::endl;
		for(const auto &node: nodes){
			std::multiset<std::pair<id, double>, multiset_value_comparator> ranking;
			_ithmm->geneerate_word_ranking_of_node(node, ranking);
			int n = 0;
			std::string indices = node->_dump_indices();
			std::cout << "\x1b[32;1m" << "[" << indices << "]" << "\x1b[0m" << std::endl;
			for(const auto &elem: ranking){
				id word_id = elem.first;
				std::wstring &word = dict->_id_to_str[word_id];
				for(int i = 0;i < std::max(0, 15 - (int)word.size());i++){
					word += L" ";
				}
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
	void Model::show_hpylm_for_each_tag(Dictionary* dict){
		auto pair = std::make_pair(0, 0);
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
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
	void Model::show_sticks(){
		std::vector<Node*> nodes;
		enumerate_all_states(nodes);
		for(const auto &node: nodes){
			std::string indices = node->_dump_indices();
			c_printf("[*]%s\n", ((boost::format("[%s]") % indices.c_str())).str().c_str());
			_show_stick(node);
		}
	}
	void Model::_show_stick(Node* node_in_structure){
		assert(node_in_structure != NULL);
		double p_eos = node_in_structure->compute_transition_probability_to_eos(_ithmm->_tau0, _ithmm->_tau1);
		_ithmm->update_stick_length_of_tssb(node_in_structure->_transition_tssb, 1.0 - p_eos, true);

		std::vector<Node*> nodes;
		node_in_structure->_transition_tssb->enumerate_nodes_from_left_to_right(nodes);
		for(const auto &node: nodes){
			std::string indices = node->_dump_indices();
			std::string tab = "";
			for(int i = 0;i < node->_depth_v;i++){
				tab += "	";
			}
			std::cout << "\x1b[32;1m" << tab << "[" << indices << "]" << "\x1b[0m " << node->_probability << std::endl;
		}
	}
	void Model::update_hyperparameters(){
		_ithmm->sample_hpylm_hyperparameters();
	}
}