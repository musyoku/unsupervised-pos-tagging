#include "ithmm.h"

namespace ithmm {
	template <class Archive>
	void iTHMM::serialize(Archive& ar, unsigned int version)
	{
		static_cast<void>(version);
		ar & _structure_tssb;
		ar & _bos_tssb;
		ar & _alpha;
		ar & _gamma;
		ar & _lambda_alpha;
		ar & _lambda_gamma;
		ar & _tau0;
		ar & _tau1;
		ar & _word_g0;
		ar & _current_max_depth;
		ar & _depth_limit;
		ar & _hpylm_d_m;
		ar & _hpylm_theta_m;
		ar & _hpylm_a_m;
		ar & _hpylm_b_m;
		ar & _hpylm_alpha_m;
		ar & _hpylm_beta_m;
	}
	iTHMM::iTHMM(){
		_alpha = sampler::uniform(iTHMM_ALPHA_MIN, iTHMM_ALPHA_MAX);
		_gamma = sampler::uniform(iTHMM_GAMMA_MIN, iTHMM_GAMMA_MAX);
		_lambda_alpha = sampler::uniform(iTHMM_LAMBDA_ALPHA_MIN, iTHMM_LAMBDA_ALPHA_MAX);
		_lambda_gamma = sampler::uniform(iTHMM_LAMBDA_GAMMA_MAX, iTHMM_LAMBDA_GAMMA_MAX);
		_strength = sampler::uniform(iTHMM_STRENGTH_MIN, iTHMM_STRENGTH_MAX);
		_tau0 = iTHMM_TAU_0;
		_tau1 = iTHMM_TAU_1;
		_current_max_depth = 0;
		_word_g0 = -1;

		_structure_tssb = new TSSB();
		_structure_tssb->_root->_owner_id_in_structure = TSSB_STRUCTURE_ID;
		_structure_tssb->_owner_id = TSSB_STRUCTURE_ID;
		Node* root_in_structure = _structure_tssb->_root;
		root_in_structure->init_hpylm();
		Node* root_in_htssb = new Node(NULL, root_in_structure->_identifier);
		root_in_htssb->_owner_id_in_structure = root_in_structure->_identifier;
		root_in_htssb->_owner_in_structure = root_in_structure;
		root_in_htssb->_structure_tssb_myself = root_in_structure;
		root_in_structure->_transition_tssb = new TSSB(root_in_htssb);
		root_in_structure->_transition_tssb->_owner_id = root_in_structure->_identifier;
		root_in_structure->_transition_tssb_myself = root_in_htssb;

		Node* root_in_bos = new Node(NULL, root_in_structure->_identifier);
		root_in_bos->_owner_id_in_structure = TSSB_BOS_ID;		// そもそも木構造上に所有者がいないが気にしない
		root_in_bos->_structure_tssb_myself = root_in_structure;
		_bos_tssb = new TSSB(root_in_bos);
		_bos_tssb->_owner_id = TSSB_BOS_ID;

		_hpylm_d_m.push_back(HPYLM_D);
		_hpylm_theta_m.push_back(HPYLM_THETA);
		_hpylm_a_m.push_back(HPYLM_A);
		_hpylm_b_m.push_back(HPYLM_B);
		_hpylm_alpha_m.push_back(HPYLM_ALPHA);
		_hpylm_beta_m.push_back(HPYLM_BETA);

		_mh_enabled = true;
		_num_mh_rejection = 0;
		_num_mh_acceptance = 0;
	}
	iTHMM::~iTHMM(){
		delete _structure_tssb;
		delete _bos_tssb;
	}
	void iTHMM::initialize_with_training_dataset(std::vector<std::vector<Word*>> &dataset){
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &data = dataset[data_index];
			if(data.size() == 0){
				continue;
			}
			// 状態路ランダムに設定
			for(int i = 0;i < data.size();i++){
				Word* word = data[i];
				Node* state = NULL;
				state = sample_node_in_tssb(_structure_tssb, true);
				assert(state != NULL);
				word->_state = state;
			}
			Node* prev_state = NULL;						// <bos>
			for(int i = 0;i < data.size();i++){
				Word* word = data[i];
				Node* state = word->_state;
				add_initial_parameters(prev_state, state, word->_id);
				prev_state = state;
			}
			add_initial_parameters(prev_state, NULL, 0);	// <eos>
		}
	}
	// デバッグ用
	// これを呼んで全パラメータが消えなかったらバグっている
	void iTHMM::remove_all_data(std::vector<std::vector<Word*>> &dataset){
		for(int data_index = 0;data_index < dataset.size();data_index++){
			std::vector<Word*> &data = dataset[data_index];
			if(data.size() == 0){
				continue;
			}
			Node* prev_state = NULL;
			for(int i = 0;i < data.size();i++){
				Word* word = data[i];
				Node* state = word->_state;
				remove_initial_parameters(prev_state, state, word->_id);
				prev_state = state;
			}
			remove_initial_parameters(prev_state, NULL, 0);
		}
	}
	void iTHMM::set_depth_limit(int limit){
		_depth_limit = limit;
	}
	void iTHMM::set_word_g0(double g0){
		_word_g0 = g0;
	}
	bool iTHMM::is_node_in_bos_tssb(Node* node){
		assert(node != NULL);
		return node->_owner_id_in_structure == TSSB_BOS_ID;
	}
	bool iTHMM::is_node_in_structure_tssb(Node* node){
		assert(node != NULL);
		return node->_owner_id_in_structure == TSSB_STRUCTURE_ID;
	}
	bool iTHMM::is_node_in_htssb(Node* node){
		assert(node != NULL);
		return is_node_in_bos_tssb(node) == false && is_node_in_structure_tssb(node) == false;
	}
	bool iTHMM::is_node_root(Node* node){
		return node->_depth_v == 0;
	}
	bool iTHMM::is_tssb_bos(TSSB* tssb){
		assert(tssb != NULL);
		return tssb->_owner_id == TSSB_BOS_ID;
	}
	bool iTHMM::is_tssb_structure(TSSB* tssb){
		assert(tssb != NULL);
		return tssb->_owner_id == TSSB_STRUCTURE_ID;
	}
	bool iTHMM::is_tssb_htssb(TSSB* tssb){
		return is_tssb_bos(tssb) == false && is_tssb_structure(tssb) == false;
	}
	bool iTHMM::is_node_to_the_left_of_node(Node* left, Node* right){
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
	Node* iTHMM::generate_and_add_new_child_to(Node* parent){
		assert(parent != NULL);
		// まず木構造上で子ノードを作る
		Node* generated_child_in_structure = NULL;
		if(is_node_in_structure_tssb(parent)){	// parentが木構造上のノードの場合
			generated_child_in_structure = parent->generate_child();
		}else{	// parentが別のノードの遷移確率用TSSB上のノードだった場合
			Node* parent_in_structure = _structure_tssb->find_node_by_tracing_horizontal_indices(parent);
			generated_child_in_structure = parent_in_structure->generate_child();
		}
		assert(generated_child_in_structure != NULL);
		// HTSSBをセット
		generated_child_in_structure->_transition_tssb = generate_transition_tssb_belonging_to(generated_child_in_structure);
		Node* myself_in_htssb = generated_child_in_structure->find_same_node_in_transition_tssb();
		assert(myself_in_htssb != NULL);
		generated_child_in_structure->_transition_tssb_myself = myself_in_htssb;
		// HPYLM
		generated_child_in_structure->init_hpylm();

		Node* return_child = generated_child_in_structure;	// 実際に返すノード
		// <bos>TSSB上で子ノードを作成
		Node* generated_child_in_bos = _generate_and_add_new_child_to_bos_tssb(generated_child_in_structure);
		if(is_node_in_bos_tssb(parent)){
			return_child = generated_child_in_bos;
		}
		// 木構造上の全ノードのHTSSBにノードを追加
		_generate_and_add_new_child_to_all_htssb(_structure_tssb->_root, parent, generated_child_in_structure, return_child);
		// ポインタを張る
		Node* generated_child_in_htssb = generated_child_in_structure->_transition_tssb_myself;
		assert(generated_child_in_htssb != NULL);
		// generated_child_in_htssb->_structure_tssb_myself = generated_child_in_structure;
		//// 木構造上の親ノードのHTSSBの自分と同じ位置のノードへのポインタ
		Node* iterator_in_structure = generated_child_in_structure;
		Node* parent_in_structure = iterator_in_structure->_parent;
		Node* iterator_in_htssb = generated_child_in_htssb;
		Node* iterator_in_parent_htssb = NULL;
		while(parent_in_structure != NULL){
			assert(iterator_in_structure->_transition_tssb_myself != NULL);
			// 木構造上での親ノードが持つHTSSBにある対応するノードを取る
			iterator_in_parent_htssb = parent_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(generated_child_in_structure);
			assert(iterator_in_parent_htssb != NULL);
			// ポインタを張る
			iterator_in_htssb->_parent_transition_tssb_myself = iterator_in_parent_htssb;
			iterator_in_htssb->_structure_tssb_myself = generated_child_in_structure;
			assert(iterator_in_htssb->_structure_tssb_myself->_identifier == generated_child_in_structure->_identifier);
			// 木構造上で次の親ノードへ
			iterator_in_structure = parent_in_structure;
			parent_in_structure = iterator_in_structure->_parent;
			iterator_in_htssb = iterator_in_parent_htssb;
		}
		// HPYLM用のハイパーパラメータを追加
		if(return_child->_depth_v > _current_max_depth){
			_current_max_depth = return_child->_depth_v;		
			while(_current_max_depth >= _hpylm_d_m.size()){
				_hpylm_d_m.push_back(HPYLM_D);
			}			
			while(_current_max_depth >= _hpylm_theta_m.size()){
				_hpylm_theta_m.push_back(HPYLM_THETA);
			}		
			while(_current_max_depth >= _hpylm_a_m.size()){
				_hpylm_a_m.push_back(HPYLM_A);
			}			
			while(_current_max_depth >= _hpylm_b_m.size()){
				_hpylm_b_m.push_back(HPYLM_B);
			}		
			while(_current_max_depth >= _hpylm_alpha_m.size()){
				_hpylm_alpha_m.push_back(HPYLM_ALPHA);
			}			
			while(_current_max_depth >= _hpylm_beta_m.size()){
				_hpylm_beta_m.push_back(HPYLM_BETA);
			}	
		}
		return return_child;
	}
	void iTHMM::_generate_and_add_new_child_to_all_htssb(Node* iterator_in_structure, Node* parent, Node* generated_child_in_structure, Node* &return_child){
		// iteratorとgenerated_childが同一の場合はすでに追加されているのでスキップ
		if(iterator_in_structure->_identifier != generated_child_in_structure->_identifier){
			assert(iterator_in_structure->_transition_tssb != NULL);
			int owner_id_of_htssb_parent_belongs = parent->_owner_id_in_structure;
			int child_id_to_generate = generated_child_in_structure->_identifier;
			// 遷移確率用TSSBの同じ位置に子ノードを挿入
			Node* parent_in_htssb = iterator_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(parent);
			assert(parent_in_htssb != NULL);
			assert(parent_in_htssb->_identifier == generated_child_in_structure->_parent->_identifier);
			Node* child_in_htssb = new Node(parent_in_htssb, child_id_to_generate);
			child_in_htssb->_structure_tssb_myself = generated_child_in_structure;
			if(child_in_htssb->_owner_id_in_structure == owner_id_of_htssb_parent_belongs){	// 親と同じTSSB上の子ノードを返す
				return_child = child_in_htssb;
			}
			parent_in_htssb->add_child(child_in_htssb);
		}
		for(const auto &child: iterator_in_structure->_children){
			_generate_and_add_new_child_to_all_htssb(child, parent, generated_child_in_structure, return_child);
		}
	}
	Node* iTHMM::_generate_and_add_new_child_to_bos_tssb(Node* generated_child_in_structure){
		// 木構造上での親ノードが<bos>TSSBのどのノードに対応するかを調べる
		Node* parent = _bos_tssb->find_node_by_tracing_horizontal_indices(generated_child_in_structure->_parent);
		assert(parent != NULL);
		Node* child = new Node(parent, generated_child_in_structure->_identifier);
		child->_owner_id_in_structure = TSSB_BOS_ID;
		parent->add_child(child);
		// ポインタを張る
		generated_child_in_structure->_bos_tssb_myself = child;
		child->_structure_tssb_myself = generated_child_in_structure;
		return child;
	}
	// 木構造上のノードにHTSSBを追加
	TSSB* iTHMM::generate_transition_tssb_belonging_to(Node* owner_in_structure){
		assert(is_node_in_structure_tssb(owner_in_structure));
		Node* root_in_structure = _structure_tssb->_root;
		Node* root_in_htssb = new Node(NULL, root_in_structure->_identifier);
		root_in_htssb->_owner_id_in_structure = owner_in_structure->_identifier;
		root_in_htssb->_owner_in_structure = owner_in_structure;
		root_in_htssb->_parent_transition_tssb_myself = NULL;
		if(owner_in_structure->_parent != NULL){
			root_in_htssb->_parent_transition_tssb_myself = owner_in_structure->_parent->_transition_tssb->_root;
		}
		root_in_htssb->_structure_tssb_myself = root_in_structure;
		copy_children_in_structure_to_transition_tssb(root_in_structure, root_in_htssb, owner_in_structure);
		TSSB* target = new TSSB(root_in_htssb);
		target->_owner_id = owner_in_structure->_identifier;
		target->_owner = owner_in_structure;
		return target;
	}
	// 生成したHTSSBを木構造と同一の形状にするために子ノードを生成・追加
	void iTHMM::copy_children_in_structure_to_transition_tssb(Node* source_in_structure, Node* target_in_htssb, Node* owner_in_structure){
		for(const auto source_child_in_structure: source_in_structure->_children){
			Node* child = new Node(target_in_htssb, source_child_in_structure->_identifier);
			child->_owner_id_in_structure = owner_in_structure->_identifier;
			child->_owner_in_structure = owner_in_structure;
			child->_structure_tssb_myself = source_child_in_structure;
			// child->_owner_id_in_structure = owner_in_structure;
			target_in_htssb->add_child(child);
			copy_children_in_structure_to_transition_tssb(source_child_in_structure, child, owner_in_structure);
		}
	}
	// コインを投げる操作を繰り返して到達したノードを返す
	Node* iTHMM::sample_node_in_tssb(TSSB* tssb, bool ignore_root){
		assert(is_tssb_structure(tssb));
		Node* node = _sample_node_in_tssb_by_iterating_node(tssb->_root, false, ignore_root);
		return node;
	}
	// HTSSB上でノードをサンプリング
	Node* iTHMM::sample_node_in_htssb(TSSB* tssb, bool ignore_root){
		assert(is_tssb_htssb(tssb));
		Node* node = _sample_node_in_tssb_by_iterating_node(tssb->_root, true, ignore_root);
		assert(node->_owner_id_in_structure == tssb->_owner_id);
		return node;
	}
	// 止まるノードを決定する
	// htssb_modeがtrueの場合、停止確率は親のHTSSBから生成する
	// htssb_modeがfalseの場合は普通のTSSBによるクラスタリング
	Node* iTHMM::_sample_node_in_tssb_by_iterating_node(Node* iterator, bool htssb_mode, bool ignore_root){
		assert(iterator != NULL);
		double head = compute_expectation_of_vertical_sbr_ratio(iterator, htssb_mode);
		iterator->_children_stick_length = iterator->_stick_length * (1 - head);
		double bernoulli = sampler::uniform(0, 1);
		if(bernoulli <= head){			// 表が出たらこのノードに降りる
			if(is_node_root(iterator) == false || (is_node_root(iterator) && ignore_root == false)){
				return iterator;
			}
		}
		// 子ノードがある場合
		for(int i = 0;i < iterator->_children.size();i++){
			Node* child = iterator->_children[i];
			assert(child != NULL);
			double head = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double bernoulli = sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _sample_node_in_tssb_by_iterating_node(child, htssb_mode, ignore_root);
			}
		}
		// ない場合生成しながらコインを投げる
		while(true){
			Node* child = generate_and_add_new_child_to(iterator);
			double head = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double bernoulli = sampler::uniform(0, 1);
			if(bernoulli <= head){		// 表が出たら次に止まるかどうかを決める
				return _sample_node_in_tssb_by_iterating_node(child, htssb_mode, ignore_root);
			}
		}
	}
	// [0, 1)の一様分布からノードをサンプリング
	// HTSSBの場合（木構造に対してはそもそも行わない）
	Node* iTHMM::retrospective_sampling(double uniform, TSSB* tssb, double total_stick_length, bool htssb_mode){
		assert(uniform < total_stick_length);
		assert( (htssb_mode && is_tssb_htssb(tssb)) || (htssb_mode == false && is_tssb_htssb(tssb) == false) );
		Node* root = tssb->_root;
		double ratio_v = compute_expectation_of_vertical_sbr_ratio(root, htssb_mode);
		double sum_probability = total_stick_length * ratio_v;
		root->_stick_length = total_stick_length;
		root->_probability = total_stick_length * ratio_v;
		root->_children_stick_length = total_stick_length * (1.0 - ratio_v);
		root->_sum_probability = sum_probability;
		// uniform = uniform * root->_children_stick_length + sum_probability; // ルートを除外する場合
		Node* node =  _retrospective_sampling_by_iterating_node(uniform, sum_probability, root, htssb_mode);
		assert(node != NULL);
		return node;
	}
	Node* iTHMM::_retrospective_sampling_by_iterating_node(double uniform, double &sum_probability, Node* iterator, bool htssb_mode){
		if(uniform < sum_probability){
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
			double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
			child->_stick_length = rest_stick_length * ratio_h;
			assert(child->_stick_length > 0);
			child->_probability = child->_stick_length * ratio_v;
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
			child->_sum_probability = sum_probability + sum_stick_length_over_children + child->_probability;
			if(uniform < sum_probability + sum_stick_length_over_children + child->_stick_length){
				// child->_stick_lengthだけだとこのノードの棒の長さ（つまりこのノード+子ノードに割り当てる棒）なので
				// ratio_vも掛けてこのノードで止まる確率にする必要がある
				sum_probability += sum_stick_length_over_children + child->_stick_length * ratio_v;
				if(uniform < sum_probability){
					return child;
				}
				return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);
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
			child->_sum_probability = sum_probability + sum_stick_length_over_children + child->_probability;
			if(uniform < sum_probability + sum_stick_length_over_children + child->_stick_length){
				sum_probability += sum_stick_length_over_children + child->_probability;
				if(uniform < sum_probability){
					return child;
				}
				return _retrospective_sampling_by_iterating_node(uniform, sum_probability, child, htssb_mode);
			}
			sum_stick_length_over_children += child->_stick_length;
			rest_stick_length *= 1.0 - ratio_h;
		}
	}
	void iTHMM::perform_gibbs_sampling_data(std::vector<Word*> &data){
		assert(data.size() > 0);
		Node* prev_state = NULL;
		Node* next_state = data.size() == 1 ? NULL : data[1]->_state;
		for(int i = 0;i < data.size();i++){
			Word* word = data[i];
			Node* state = word->_state;
			remove_parameters(prev_state, state, next_state, word->_id);
			state = draw_state(prev_state, state, next_state, word->_id);
			add_parameters(prev_state, state, next_state, word->_id);
			prev_state = state;
			next_state = i < data.size() - 2 ? data[i + 2]->_state : NULL;
			word->_state = state;
		}
	}
	// データ読み込み時の状態初期化時にのみ使う
	void iTHMM::add_initial_parameters(Node* prev_state_in_structure, Node* state_in_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(state_in_structure->_transition_tssb != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			add_customer_to_tssb_node(state_in_bos);
			add_customer_to_tssb_node(state_in_structure);			// 参照カウント用
			add_customer_to_hpylm(state_in_structure, word_id);
			return;
		}
		// <eos>への遷移を含む場合
		if(state_in_structure == NULL){
			assert(prev_state_in_structure != NULL);
			assert(prev_state_in_structure->_transition_tssb != NULL);
			assert(is_node_in_structure_tssb(prev_state_in_structure));
			prev_state_in_structure->increment_transition_count_to_eos();
			return;
		}
		assert(prev_state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));

		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
		add_customer_to_htssb_node(state_in_prev_state_htssb);
		add_customer_to_tssb_node(state_in_structure);			// 参照カウント用
		add_customer_to_hpylm(state_in_structure, word_id);
		prev_state_in_structure->increment_transition_count_to_other();
	}
	void iTHMM::add_temporal_parameters(Node* prev_state_in_structure, Node* state_in_structure){
		assert(prev_state_in_structure != NULL);
		assert(state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));

		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
		add_customer_to_htssb_node(state_in_prev_state_htssb);
		add_customer_to_tssb_node(state_in_structure);			// 参照カウント用
		prev_state_in_structure->increment_transition_count_to_other();
	}
	void iTHMM::add_parameters(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(state_in_structure->_transition_tssb != NULL);
			assert(next_state_in_structure != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			assert(is_node_in_structure_tssb(next_state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			add_customer_to_tssb_node(state_in_bos);
			add_customer_to_tssb_node(state_in_structure);			// 参照カウント用
			Node* next_state_in_state_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
			assert(next_state_in_state_htssb != NULL);
			assert(next_state_in_structure->_identifier == next_state_in_state_htssb->_identifier);
			add_customer_to_htssb_node(next_state_in_state_htssb);
			add_customer_to_tssb_node(next_state_in_structure);		// 参照カウント用
			add_customer_to_hpylm(state_in_structure, word_id);
			state_in_structure->increment_transition_count_to_other();
			return;
		}
		// <eos>への遷移を含む場合
		if(next_state_in_structure == NULL){
			assert(prev_state_in_structure != NULL);
			assert(prev_state_in_structure->_transition_tssb != NULL);
			assert(state_in_structure != NULL);
			assert(is_node_in_structure_tssb(prev_state_in_structure));
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_prev_state_htssb != NULL);
			assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
			add_customer_to_htssb_node(state_in_prev_state_htssb);
			add_customer_to_tssb_node(state_in_structure);		// 参照カウント用
			add_customer_to_hpylm(state_in_structure, word_id);
			prev_state_in_structure->increment_transition_count_to_other();
			state_in_structure->increment_transition_count_to_eos();
			return;
		}
		// <bos>と<eos>両方を含む場合
		if(prev_state_in_structure == NULL && next_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			add_customer_to_tssb_node(state_in_bos);
			add_customer_to_tssb_node(state_in_structure);			// 参照カウント用
			add_customer_to_hpylm(state_in_structure, word_id);
			state_in_structure->increment_transition_count_to_eos();
			return;
		}
		assert(prev_state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(next_state_in_structure != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));
		assert(is_node_in_structure_tssb(next_state_in_structure));

		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
		add_customer_to_htssb_node(state_in_prev_state_htssb);
		add_customer_to_tssb_node(state_in_structure);			// 参照カウント用

		Node* next_state_in_state_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
		assert(next_state_in_state_htssb != NULL);
		assert(next_state_in_structure->_identifier == next_state_in_state_htssb->_identifier);
		add_customer_to_htssb_node(next_state_in_state_htssb);
		add_customer_to_tssb_node(next_state_in_structure);		// 参照カウント用

		add_customer_to_hpylm(state_in_structure, word_id);
		prev_state_in_structure->increment_transition_count_to_other();
		state_in_structure->increment_transition_count_to_other();
	}
	// データをモデルから全部消す時はこれを使う
	void iTHMM::remove_initial_parameters(Node* prev_state_in_structure, Node* state_in_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(state_in_structure->_transition_tssb != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			remove_customer_from_tssb_node(state_in_bos);
			remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用
			remove_customer_from_hpylm(state_in_structure, word_id);
			return;
		}
		// <eos>への遷移を含む場合
		if(state_in_structure == NULL){
			assert(prev_state_in_structure != NULL);
			assert(prev_state_in_structure->_transition_tssb != NULL);
			assert(is_node_in_structure_tssb(prev_state_in_structure));
			prev_state_in_structure->decrement_transition_count_to_eos();
			return;
		}
		assert(prev_state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));

		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
		remove_customer_from_htssb_node(state_in_prev_state_htssb);
		remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用
		remove_customer_from_hpylm(state_in_structure, word_id);
		prev_state_in_structure->decrement_transition_count_to_other();
	}
	void iTHMM::remove_temporal_parameters(Node* prev_state_in_structure, Node* state_in_structure){
		assert(prev_state_in_structure != NULL);
		assert(state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));
		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		assert(state_in_structure->_identifier == state_in_prev_state_htssb->_identifier);
		remove_customer_from_htssb_node(state_in_prev_state_htssb, true);
		remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用
		prev_state_in_structure->decrement_transition_count_to_other();
	}
	void iTHMM::remove_parameters(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id){
		// <bos>からの遷移を含む場合
		if(prev_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(state_in_structure->_transition_tssb != NULL);
			assert(next_state_in_structure != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			assert(is_node_in_structure_tssb(next_state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			remove_customer_from_tssb_node(state_in_bos);
			remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用
			Node* next_state_in_state_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
			assert(next_state_in_state_htssb != NULL);
			remove_customer_from_htssb_node(next_state_in_state_htssb);
			remove_customer_from_tssb_node(next_state_in_structure);		// 参照カウント用
			remove_customer_from_hpylm(state_in_structure, word_id);
			state_in_structure->decrement_transition_count_to_other();
			return;
		}
		// <eos>への遷移を含む場合
		if(next_state_in_structure == NULL){
			assert(prev_state_in_structure != NULL);
			assert(prev_state_in_structure->_transition_tssb != NULL);
			assert(state_in_structure != NULL);
			assert(is_node_in_structure_tssb(prev_state_in_structure));
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_prev_state_htssb != NULL);
			remove_customer_from_htssb_node(state_in_prev_state_htssb);
			remove_customer_from_tssb_node(state_in_structure);		// 参照カウント用
			remove_customer_from_hpylm(state_in_structure, word_id);
			prev_state_in_structure->decrement_transition_count_to_other();
			state_in_structure->decrement_transition_count_to_eos();
			return;
		}
		// <bos>と<eos>両方を含む場合
		if(prev_state_in_structure == NULL && next_state_in_structure == NULL){
			assert(state_in_structure != NULL);
			assert(is_node_in_structure_tssb(state_in_structure));
			Node* state_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
			assert(state_in_bos);
			remove_customer_from_tssb_node(state_in_bos);
			remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用
			remove_customer_from_hpylm(state_in_structure, word_id);
			state_in_structure->decrement_transition_count_to_eos();
			return;
		}
		assert(prev_state_in_structure != NULL);
		assert(prev_state_in_structure->_transition_tssb != NULL);
		assert(state_in_structure != NULL);
		assert(state_in_structure->_transition_tssb != NULL);
		assert(next_state_in_structure != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));
		assert(is_node_in_structure_tssb(next_state_in_structure));

		Node* state_in_prev_state_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
		assert(state_in_prev_state_htssb != NULL);
		remove_customer_from_htssb_node(state_in_prev_state_htssb);
		remove_customer_from_tssb_node(state_in_structure);			// 参照カウント用

		Node* next_state_in_state_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
		assert(next_state_in_state_htssb != NULL);
		remove_customer_from_htssb_node(next_state_in_state_htssb);
		remove_customer_from_tssb_node(next_state_in_structure);		// 参照カウント用

		remove_customer_from_hpylm(state_in_structure, word_id);
		prev_state_in_structure->decrement_transition_count_to_other();
		state_in_structure->decrement_transition_count_to_other();
	}
	// 新しい状態のギブスサンプリング
	// なるべく論文の記号を使う
	Node* iTHMM::draw_state(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id){
		if(prev_state_in_structure == NULL){
			return _draw_state_from_bos(state_in_structure, next_state_in_structure, word_id);
		}
		if(next_state_in_structure == NULL){
			return _draw_state_to_eos(prev_state_in_structure, state_in_structure, word_id);
		}
		return _draw_state(prev_state_in_structure, state_in_structure, next_state_in_structure, word_id);
	}
	Node* iTHMM::_draw_state(Node* prev_state_in_structure, Node* state_in_structure, Node* next_state_in_structure, id word_id){
		assert(prev_state_in_structure != NULL);
		assert(state_in_structure != NULL);
		assert(next_state_in_structure != NULL);
		assert(is_node_in_structure_tssb(state_in_structure));
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(next_state_in_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_in_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		// s_tから<eos>へ接続する確率
		assert(state_in_structure->_transition_tssb != NULL);
		Node* next_state_in_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
		assert(next_state_in_htssb != NULL);
		assert(next_state_in_htssb->_identifier == next_state_in_structure->_identifier);
		// <eos>以外に接続する確率を棒全体の長さとする
		// if(state_in_structure->_identifier == next_state_in_structure->_identifier){
			// s_t == s_{t+1}の場合は正しい確率を求めるためにp(s_t|s_{t-1})に客を追加
			// word_idは使わないので何を指定しても良い
			add_temporal_parameters(prev_state_in_structure, state_in_structure);
		// }
		double Peos_given_s = state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		double Pnext_given_s = (1.0 - Peos_given_s) * compute_node_probability_in_tssb(state_in_structure->_transition_tssb, next_state_in_htssb, 1.0);
		assert(0 < Pnext_given_s && Pnext_given_s <= 1);

		// if(state_in_structure->_identifier == next_state_in_structure->_identifier){
			remove_temporal_parameters(prev_state_in_structure, state_in_structure);
		// }

		// スライス
		double slice = Pw_given_s * Pnext_given_s * sampler::uniform(0, 1);
		assert(slice > 0);

		double st = 0;
		double ed = 1;

		while(true){
			double u = sampler::uniform(st, ed);
			if( (st <= u && u < ed) == false){	// 見つからなかったら元の状態を返す
				return state_in_structure;
			}
			// assert(st <= u && u < ed);
			Node* new_state_in_prev_htssb = retrospective_sampling(u, prev_state_in_structure->_transition_tssb, 1.0, true);
			assert(new_state_in_prev_htssb != NULL);
			Node* new_state_in_structure = new_state_in_prev_htssb->_structure_tssb_myself;
			assert(new_state_in_structure != NULL);
			assert(new_state_in_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_in_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			if(new_state_in_structure->_identifier == state_in_structure->_identifier){
				assert(new_Pw_given_s == Pw_given_s);	// 一致しないならバグ
			}
			// 遷移確率
			//// s_{new}からs_{t+1}へ接続する確率
			Node* next_state_in_new_state_htssb = new_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
			//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
			double new_Pnext_given_s;
			double new_Peos_given_s;
			if(new_state_in_structure->_identifier == state_in_structure->_identifier){
				new_Pnext_given_s = Pnext_given_s;
				new_Peos_given_s = Peos_given_s;
			}else{
				add_temporal_parameters(prev_state_in_structure, new_state_in_structure);
				new_Peos_given_s = new_state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
				new_Pnext_given_s = (1.0 - new_Peos_given_s) * compute_node_probability_in_tssb(new_state_in_structure->_transition_tssb, next_state_in_new_state_htssb, 1.0);
				remove_temporal_parameters(prev_state_in_structure, new_state_in_structure);
			}
			assert(0 < new_Pnext_given_s && new_Pnext_given_s <= 1);

			// 尤度を計算
			double likelihoood = new_Pw_given_s * new_Pnext_given_s;

			if(likelihoood > slice){
				if(_mh_enabled == false){
					_num_mh_acceptance += 1;
					return new_state_in_structure;
				}
				// メトロポリス・ヘイスティングス法による補正
				// 未完成なためサンプリングした新しい品詞をそのまま返す
				_num_mh_acceptance += 1;
				return new_state_in_structure;
				Node* state_in_prev_htssb = prev_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
				assert(state_in_prev_htssb != NULL);
				assert(state_in_prev_htssb->_identifier == state_in_structure->_identifier);
				double Peos_given_prev = prev_state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
				double Ps_given_prev = (1.0 - Peos_given_prev) * compute_node_probability_in_tssb(prev_state_in_structure->_transition_tssb, state_in_prev_htssb, 1.0);
				double Pnew_s_given_prev = (1.0 - Peos_given_prev) * compute_node_probability_in_tssb(prev_state_in_structure->_transition_tssb, new_state_in_prev_htssb, 1.0);

				Node* state_in_root_htssb = _structure_tssb->_root->_transition_tssb->find_node_by_tracing_horizontal_indices(state_in_structure);
				assert(state_in_root_htssb != NULL);
				assert(state_in_root_htssb->_identifier == state_in_structure->_identifier);
				Node* new_state_in_root_htssb = _structure_tssb->_root->_transition_tssb->find_node_by_tracing_horizontal_indices(new_state_in_structure);
				assert(new_state_in_root_htssb != NULL);
				assert(new_state_in_root_htssb->_identifier == new_state_in_structure->_identifier);
				compute_node_probability_in_tssb(_structure_tssb->_root->_transition_tssb, state_in_root_htssb, 1.0);
				compute_node_probability_in_tssb(_structure_tssb->_root->_transition_tssb, new_state_in_root_htssb, 1.0);
				double Ps = state_in_root_htssb->_probability;
				double Pnew_s = new_state_in_root_htssb->_probability;
				double Pw_given_new_s = compute_Pw_given_s(word_id, new_state_in_structure);
				assert(Ps_given_prev * Pw_given_new_s > 0);
				// 採択率の計算式が不明
				double adoption = Pnew_s_given_prev / Ps_given_prev;
				adoption = std::min(1.0, adoption);
				double u = sampler::uniform(0, 1);
				if(u <= adoption){
					_num_mh_acceptance += 1;
					return new_state_in_structure;
				}
				_num_mh_rejection += 1;
				return state_in_structure;
			}
			assert(new_state_in_structure->_identifier != state_in_structure->_identifier);	// 同じになる場合バグっている
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_in_structure, state_in_structure)){
				assert(new_state_in_prev_htssb->_sum_probability >= u);
				st = new_state_in_prev_htssb->_sum_probability;
			}else{
				assert(new_state_in_prev_htssb->_sum_probability >= u);
				assert(new_state_in_prev_htssb->_sum_probability - new_state_in_prev_htssb->_probability <= u);
				ed = new_state_in_prev_htssb->_sum_probability - new_state_in_prev_htssb->_probability;
			}
		}
	}
	Node* iTHMM::_draw_state_from_bos(Node* state_in_structure, Node* next_state_in_structure, id word_id){
		assert(state_in_structure != NULL);
		assert(next_state_in_structure != NULL);
		assert(is_node_in_structure_tssb(state_in_structure));
		assert(is_node_in_structure_tssb(next_state_in_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_in_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		//// s_tから<eos>へ接続する確率
		double Peos_given_s = state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		assert(state_in_structure->_transition_tssb != NULL);
		Node* next_state_in_htssb = state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
		assert(next_state_in_htssb != NULL);
		//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double stick_length = 1.0 - Peos_given_s;
		double Pt_given_s = compute_node_probability_in_tssb(state_in_structure->_transition_tssb, next_state_in_htssb, 1.0);
		Pt_given_s *= stick_length;
		assert(0 < Pt_given_s && Pt_given_s <= 1);
		// スライス
		double slice = Pw_given_s * Pt_given_s * sampler::uniform(0, 1);
		assert(slice > 0);

		double st = 0;
		double ed = 1;
		while(true){
			double u = sampler::uniform(st, ed);
			assert(st <= u && u < ed);
			Node* new_state_in_bos = retrospective_sampling(u, _bos_tssb, 1.0, false);
			assert(is_node_in_bos_tssb(new_state_in_bos));
			assert(new_state_in_bos != NULL);
			Node* new_state_in_structure = new_state_in_bos->_structure_tssb_myself;
			assert(new_state_in_structure != NULL);
			assert(new_state_in_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_in_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			// 遷移確率
			//// s_{new}からs_{t+1}へ接続する確率
			Node* next_state_in_new_state_htssb = new_state_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(next_state_in_structure);
			double Peos_given_new_s = new_state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
			double new_Pnext_given_s = compute_node_probability_in_tssb(new_state_in_structure->_transition_tssb, next_state_in_new_state_htssb, 1.0);
			//// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
			double total_stick_length_of_new_tssb = 1.0 - Peos_given_new_s;
			new_Pnext_given_s *= total_stick_length_of_new_tssb;
			assert(0 < new_Pnext_given_s && new_Pnext_given_s <= 1);
			// 尤度を計算
			double likelihoood = new_Pw_given_s * new_Pnext_given_s;
			if(likelihoood > slice){
				return new_state_in_structure;
			}
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_in_structure, state_in_structure)){
				st = u;
			}else{
				ed = u;
			}
		}
	}
	Node* iTHMM::_draw_state_to_eos(Node* prev_state_in_structure, Node* state_in_structure, id word_id){
		assert(prev_state_in_structure != NULL);
		assert(state_in_structure != NULL);
		assert(is_node_in_structure_tssb(prev_state_in_structure));
		assert(is_node_in_structure_tssb(state_in_structure));
		// 出力確率
		double Pw_given_s = compute_Pw_given_s(word_id, state_in_structure);
		assert(0 < Pw_given_s && Pw_given_s <= 1);
		// 遷移確率
		//// s_tから<eos>へ接続する確率
		double Peos_given_s = state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		assert(0 < Peos_given_s && Peos_given_s <= 1);
		// スライス
		double slice = Pw_given_s * Peos_given_s * sampler::uniform(0, 1);
		assert(slice > 0);
		// // s_{t-1}から<eos>へ接続する確率
		// double Peos_given_prev_s = prev_state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
		// <eos>以外に接続する確率を棒全体の長さとし、TSSBで分配
		double st = 0;
		double ed = 1;
		// c_printf("[r]%s\n", "_draw_state_to_eos");
		// state_in_structure->dump();
		while(true){
			double u = sampler::uniform(st, ed);
			assert(st <= u && u < ed);
			Node* new_state_in_htssb = retrospective_sampling(u, prev_state_in_structure->_transition_tssb, 1.0, true);
			// prev_state_in_structure->_transition_tssb->dump();
			// cout << u << endl;
			// new_state_in_htssb->dump();
			assert(new_state_in_htssb != NULL);
			Node* new_state_in_structure = new_state_in_htssb->_structure_tssb_myself;
			assert(new_state_in_structure != NULL);
			assert(new_state_in_structure->_transition_tssb != NULL);

			// 出力確率
			double new_Pw_given_s = compute_Pw_given_s(word_id, new_state_in_structure);
			assert(0 < new_Pw_given_s && new_Pw_given_s <= 1);
			// 遷移確率
			//// s_{new}から<eos>へ接続する確率
			double Peos_given_new_s = new_state_in_structure->compute_transition_probability_to_eos(_tau0, _tau1);
			assert(0 < Peos_given_new_s && Peos_given_new_s <= 1);
			// 尤度を計算
			double likelihoood = new_Pw_given_s * Peos_given_new_s;
			if(likelihoood > slice){
				return new_state_in_structure;
			}
			// 辞書順で前にあるかどうか
			if(is_node_to_the_left_of_node(new_state_in_structure, state_in_structure)){
				st = u;
			}else{
				ed = u;
			}
		}
	}
	void iTHMM::add_customer_to_hpylm(Node* target_in_structure, id token_id){
		assert(target_in_structure != NULL);
		assert(is_node_in_structure_tssb(target_in_structure));	// 木構造上のノードのみ
		assert(target_in_structure->_depth_v == target_in_structure->_hpylm->_depth);
		assert(_hpylm_d_m.size() > target_in_structure->_depth_v);
		assert(_hpylm_theta_m.size() > target_in_structure->_depth_v);
		assert(_word_g0 > 0);
		target_in_structure->_hpylm->add_customer(token_id, _word_g0, _hpylm_d_m, _hpylm_theta_m);
		target_in_structure->increment_word_assignment(token_id);	// 参照カウント用
	}
	void iTHMM::add_customer_to_tssb_node(Node* target_in_tssb){
		assert(target_in_tssb != NULL);
		assert(is_node_in_htssb(target_in_tssb) == false);
		double alpha = _alpha * pow(_lambda_alpha, target_in_tssb->_depth_v);
		double gamma = _gamma * pow(_lambda_gamma, std::max(0, target_in_tssb->_depth_v - 1));
		bool new_table_generated = false;
		// double ratio_v = compute_expectation_of_vertical_tssb_sbr_ratio(target_in_tssb);		// 特に計算しても意味はない
		target_in_tssb->add_customer_to_vertical_crp(alpha, 0, new_table_generated);
		// double ratio_h = compute_expectation_of_horizontal_tssb_sbr_ratio(target_in_tssb);		// 特に計算しても意味はない
		target_in_tssb->add_customer_to_horizontal_crp(gamma, 0, new_table_generated);
		// 総客数のインクリメント
		if(is_node_in_structure_tssb(target_in_tssb)){
			_structure_tssb->increment_num_customers();
		}else if(is_node_in_bos_tssb(target_in_tssb)){
			_bos_tssb->increment_num_customers();
		}
		// 参照カウントのインクリメント
		//// <bos>からの接続のカウント
		if(is_node_in_bos_tssb(target_in_tssb)){
			Node* target_in_structure = target_in_tssb->_structure_tssb_myself;
			assert(target_in_structure != NULL);
			target_in_structure->increment_ref_count();
		}
	}
	void iTHMM::add_customer_to_htssb_node(Node* target_in_htssb){
		assert(target_in_htssb != NULL);
		assert(is_node_in_htssb(target_in_htssb));
		_add_customer_to_htssb_vertical_crp(_strength, target_in_htssb);
		_add_customer_to_htssb_horizontal_crp(_gamma, target_in_htssb);
	}
	void iTHMM::_add_customer_to_htssb_vertical_crp(double alpha, Node* iterator){
		assert(iterator != NULL);
		assert(is_node_in_htssb(iterator));
		Node* iterator_in_parent_htssb = iterator->_parent_transition_tssb_myself;
		double ratio_v = 0;	// 親の場合はテーブルの増加は無視してよい
		if(iterator_in_parent_htssb != NULL){
			ratio_v = compute_expectation_of_vertical_htssb_sbr_ratio(iterator_in_parent_htssb);	// g0は親の停止確率なので注意
		}
		bool new_table_generated = false;
		iterator->add_customer_to_vertical_crp(alpha, ratio_v, new_table_generated);
		// 総客数のインクリメント
		Node* owner_in_structure = iterator->_owner_in_structure;
		assert(owner_in_structure != NULL);
		TSSB* htssb = owner_in_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_in_structure);
		htssb->increment_num_customers();
		// 参照カウントのインクリメント
		Node* iterator_in_structure = iterator->_structure_tssb_myself;
		assert(iterator_in_structure != NULL);
		iterator_in_structure->increment_ref_count();
		// 親TSSBに代理客を追加
		if(new_table_generated && iterator_in_parent_htssb != NULL){
			_add_customer_to_htssb_vertical_crp(alpha, iterator_in_parent_htssb);
		}
	}
	void iTHMM::_add_customer_to_htssb_horizontal_crp(double gamma, Node* iterator){
		assert(iterator != NULL);
		assert(is_node_in_htssb(iterator));
		Node* iterator_in_parent_htssb = iterator->_parent_transition_tssb_myself;
		double ratio_h = 0;	// 親の場合はテーブルの増加は無視してよい
		if(iterator_in_parent_htssb != NULL){
			ratio_h = compute_expectation_of_horizontal_htssb_sbr_ratio(iterator_in_parent_htssb);	// g0は親の停止確率なので注意
		}
		bool new_table_generated = false;
		iterator->add_customer_to_horizontal_crp(gamma, ratio_h, new_table_generated);
		// 総客数のインクリメント
		Node* owner_in_structure = iterator->_owner_in_structure;
		assert(owner_in_structure != NULL);
		TSSB* htssb = owner_in_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_in_structure);
		for(int d = 0;d <= iterator->_depth_v;d++){	// 水平方向には深さの数だけ別のSBRがあり、別の客が追加されることに注意
			htssb->increment_num_customers();
		}
		// 参照カウントのインクリメント
		Node* iterator_in_structure = iterator->_structure_tssb_myself;
		assert(iterator_in_structure != NULL);
		iterator_in_structure->increment_ref_count();
		// 親TSSBに代理客を追加
		if(new_table_generated && iterator_in_parent_htssb != NULL){
			_add_customer_to_htssb_horizontal_crp(gamma, iterator_in_parent_htssb);
		}
	}
	void iTHMM::remove_customer_from_hpylm(Node* target_in_structure, id token_id){
		assert(target_in_structure != NULL);
		assert(is_node_in_structure_tssb(target_in_structure));	// 木構造上のノードのみ
		assert(target_in_structure->_depth_v == target_in_structure->_hpylm->_depth);
		target_in_structure->_hpylm->remove_customer(token_id);
		target_in_structure->decrement_word_assignment(token_id);	// 参照カウント用
	}
	void iTHMM::remove_customer_from_tssb_node(Node* target_in_tssb){
		assert(target_in_tssb != NULL);
		assert(is_node_in_htssb(target_in_tssb) == false);
		bool empty_table_deleted = false;
		target_in_tssb->remove_customer_from_vertical_crp(empty_table_deleted);
		target_in_tssb->remove_customer_from_horizontal_crp(empty_table_deleted);
		// 総客数のインクリメント
		if(is_node_in_structure_tssb(target_in_tssb)){
			_structure_tssb->decrement_num_customers();
		}else if(is_node_in_bos_tssb(target_in_tssb)){
			_bos_tssb->decrement_num_customers();
		}
		// 参照カウントのインクリメント
		//// <bos>からの接続のカウント
		if(is_node_in_bos_tssb(target_in_tssb)){
			Node* target_in_structure = target_in_tssb->_structure_tssb_myself;
			assert(target_in_structure != NULL);
			target_in_structure->decrement_ref_count();
		}
	}
	void iTHMM::remove_customer_from_htssb_node(Node* target_in_htssb, bool remove_last_customer){
		assert(target_in_htssb != NULL);
		assert(is_node_in_htssb(target_in_htssb));
		_remove_customer_from_htssb_vertical_crp(target_in_htssb, remove_last_customer);
		_remove_customer_from_htssb_horizontal_crp(target_in_htssb, remove_last_customer);
	}
	void iTHMM::_remove_customer_from_htssb_vertical_crp(Node* iterator, bool remove_last_customer){
		assert(iterator != NULL);
		assert(is_node_in_htssb(iterator));
		bool empty_table_deleted = false;
		if(remove_last_customer){
			iterator->remove_last_customer_from_vertical_crp(empty_table_deleted);
		}else{
			iterator->remove_customer_from_vertical_crp(empty_table_deleted);
		}
		// 総客数のインクリメント
		Node* owner_in_structure = iterator->_owner_in_structure;
		assert(owner_in_structure != NULL);
		TSSB* htssb = owner_in_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_in_structure);
		htssb->decrement_num_customers();
		// 参照カウントのインクリメント
		Node* iterator_in_structure = iterator->_structure_tssb_myself;
		assert(iterator_in_structure != NULL);
		iterator_in_structure->decrement_ref_count();
		// 親TSSBから代理客を削除
		Node* iterator_in_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(empty_table_deleted && iterator_in_parent_htssb != NULL){
			_remove_customer_from_htssb_vertical_crp(iterator_in_parent_htssb, remove_last_customer);
		}
	}
	void iTHMM::_remove_customer_from_htssb_horizontal_crp(Node* iterator, bool remove_last_customer){
		assert(iterator != NULL);
		bool empty_table_deleted = false;
		if(remove_last_customer){
			iterator->remove_last_customer_from_horizontal_crp(empty_table_deleted);
		}else{
			iterator->remove_customer_from_horizontal_crp(empty_table_deleted);
		}
		// 総客数のインクリメント
		Node* owner_in_structure = iterator->_owner_in_structure;
		assert(owner_in_structure != NULL);
		TSSB* htssb = owner_in_structure->_transition_tssb;
		assert(htssb != NULL);
		assert(htssb->_owner_id == iterator->_owner_id_in_structure);
		for(int d = 0;d <= iterator->_depth_v;d++){	// 水平方向には深さの数だけ別のSBRがあり、別の客が追加されることに注意
			htssb->decrement_num_customers();
		}
		// 参照カウントのインクリメント
		Node* iterator_in_structure = iterator->_structure_tssb_myself;
		assert(iterator_in_structure != NULL);
		iterator_in_structure->decrement_ref_count();
		// 親TSSBから代理客を削除
		Node* iterator_in_parent_htssb = iterator->_parent_transition_tssb_myself;
		if(empty_table_deleted && iterator_in_parent_htssb != NULL){
			_remove_customer_from_htssb_horizontal_crp(iterator_in_parent_htssb, remove_last_customer);
		}
	}
	// update_stick_length_of_tssbは全ノードを更新するのに対しこっちは対象ノードのみ正確に計算する
	double iTHMM::compute_node_probability_in_tssb(TSSB* tssb, Node* node, double total_stick_length){
		assert(tssb->_owner_id == node->_owner_id_in_structure);
		bool htssb_mode = is_tssb_htssb(tssb);
		Node* iterator = tssb->_root;
		iterator->_stick_length = total_stick_length;
		double ratio_v = compute_expectation_of_vertical_sbr_ratio(iterator, htssb_mode);
		iterator->_probability = iterator->_stick_length * ratio_v;
		iterator->_children_stick_length = iterator->_stick_length * (1.0 - ratio_v);
		double min_probability = 1;		// 0が返るのを防ぐ
		for(int n = 0;n < node->_depth_v;n++){
			int depth_h = node->_horizontal_indices_from_root[n];
			double rest_stick_length = iterator->_children_stick_length;
			for(int m = 0;m <= depth_h;m++){
				Node* child = iterator->_children[m];
				double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
				double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
				child->_stick_length = rest_stick_length * ratio_h;
				double probability = child->_stick_length * ratio_v;
				if(probability == 0){		// ものすごく深いノードがサンプリングされた時にdouble型の限界を超える
					child->_probability = min_probability;
					c_printf("[r]%s\n", "stick length == 0");
					std::cout << "fixed to " << min_probability << std::endl;
				}else{
					child->_probability = probability;
					if(probability < min_probability){
						min_probability = probability;
					}
				}
				child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);
				rest_stick_length *= (1.0 - ratio_h);
			}
			iterator = iterator->_children[depth_h];
		}
		return iterator->_probability;
	}
	double iTHMM::compute_expectation_of_vertical_sbr_ratio(Node* iterator, bool htssb_mode){
		if(_depth_limit > 0){
			if(iterator->_depth_v >= _depth_limit){
				return 1;
			}
		}
		if(htssb_mode){
			assert(is_node_in_htssb(iterator));
			return compute_expectation_of_vertical_htssb_sbr_ratio(iterator);
		}
		assert(is_node_in_htssb(iterator) == false);
		return compute_expectation_of_vertical_tssb_sbr_ratio(iterator);
	}
	double iTHMM::compute_expectation_of_horizontal_sbr_ratio(Node* iterator, bool htssb_mode){
		if(htssb_mode){
			assert(is_node_in_htssb(iterator));
			return compute_expectation_of_horizontal_htssb_sbr_ratio(iterator);
		}
		assert(is_node_in_htssb(iterator) == false);
		return compute_expectation_of_horizontal_tssb_sbr_ratio(iterator);
	}
	double iTHMM::compute_expectation_of_vertical_tssb_sbr_ratio(Node* target_in_tssb){
		int pass_count = target_in_tssb->_pass_count_v;
		int stop_count = target_in_tssb->_stop_count_v;
		double alpha = _alpha * pow(_lambda_alpha, target_in_tssb->_depth_v);
		return (1.0 + stop_count) / (1.0 + alpha + stop_count + target_in_tssb->_pass_count_v);
	}
	double iTHMM::compute_expectation_of_horizontal_tssb_sbr_ratio(Node* target_in_tssb){
		int pass_count = target_in_tssb->_pass_count_h;
		int stop_count = target_in_tssb->_stop_count_h;
		double gamma = _gamma * pow(_lambda_gamma, std::max(0, target_in_tssb->_depth_v - 1));
		return (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count);
	}
	double iTHMM::compute_expectation_of_vertical_htssb_sbr_ratio(Node* target_in_htssb){
		// c_printf("[*]%s\n", "compute_expectation_of_vertical_htssb_sbr_ratio");
		assert(target_in_htssb != NULL);
		assert(is_node_in_htssb(target_in_htssb));	// 木構造上のノードだった場合は計算できない
		Node* owner_in_structure = target_in_htssb->_owner_in_structure;
		assert(owner_in_structure != NULL);
		double sbr_ratio = -1;
		// target_in_htssb->dump();

		// 木構造上での基準となるノードを選ぶ
		assert(owner_in_structure != NULL);
		// 階層TSSBなので親TSSBのSBPを全て睿珊しないと次ノードのSBPを計算できない
		int num_itr_in_structure = owner_in_structure->_depth_v + 1;
		int num_itr_in_htssb = target_in_htssb->_depth_v + 1;
		// キャッシュ用配列
		double* stop_ratio_over_parent = target_in_htssb->_stop_ratio_v_over_parent;
		double* stop_probability_over_parent = target_in_htssb->_stop_probability_v_over_parent;
		double parent_ratio_v = 0;
		// 計算
		for(int n = 0;n < num_itr_in_structure;n++){
			// cout << "n = " << n << endl;
			// 木構造上での基準となるノードを選ぶ
			// nが増えるごとに木構造を下に降りていく
			Node* iterator_in_structure = owner_in_structure->_nodes_from_root_to_myself[n];
			assert(iterator_in_structure != NULL);
			assert(iterator_in_structure->_transition_tssb != NULL);
			// iterator_in_structure->dump();
			// トップレベルのノードから順に停止確率を計算
			double sum_parent_stop_probability = 0;
			Node* iterator_in_htssb = iterator_in_structure->_transition_tssb->_root;
			assert(iterator_in_htssb != NULL);
			for(int m = 0;m < num_itr_in_htssb;m++){
				assert(iterator_in_htssb != NULL);
				// cout << "m = " << m << endl;
				// iterator_in_htssb->dump();
				if(n == 0){	// 木構造の親ノードの場合
					int pass_count = iterator_in_htssb->_pass_count_v;
					int stop_count = iterator_in_htssb->_stop_count_v;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					// cout << "depth = " << iterator_in_htssb->_depth_v << endl;
					double alpha = _alpha * pow(_lambda_alpha, iterator_in_htssb->_depth_v);
					// cout << "alpha = " << alpha << endl;
					double ratio_v = (1.0 + stop_count) / (1.0 + alpha + stop_count + pass_count + EPS);
					// cout << "ratio_v = " << ratio_v << endl;
					// assert(ratio_v < 1);
					stop_ratio_over_parent[m] = ratio_v;
					sbr_ratio = ratio_v;
				}else{	// 親の遷移確率用HTSSBから生成
					int pass_count = iterator_in_htssb->_pass_count_v;
					int stop_count = iterator_in_htssb->_stop_count_v;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					// cout << "depth = " << iterator_in_htssb->_depth_v << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					parent_ratio_v = (stop_ratio_over_parent[m] > 0) ?  stop_ratio_over_parent[m] : parent_ratio_v;
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					// ここでのαは集中度であることに注意
					double ratio_v = (_strength * parent_stop_probability + stop_count) / (_strength * (1.0 - sum_parent_stop_probability) + stop_count + pass_count + EPS);
					assert(ratio_v >= 0);
					// cout << "ratio_v = " << ratio_v << endl;
					// assert(ratio_v < 1);
					stop_ratio_over_parent[m] = ratio_v;
					sum_parent_stop_probability += parent_stop_probability;
					sbr_ratio = ratio_v;
				}
				// HTSSB上で親から子へ降りていく
				if(m < num_itr_in_htssb - 1){
					int index_h = target_in_htssb->_horizontal_indices_from_root[m];
					assert(index_h < iterator_in_htssb->_children.size());
					iterator_in_htssb = iterator_in_htssb->_children[index_h];
				}
			}
			// 計算した棒を折る比率から確率を計算
			if(n < num_itr_in_structure - 1){
				double rest_stick_length = 1;
				for(int m = 0;m < target_in_htssb->_depth_v + 1;m++){
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
		if(sbr_ratio <= 0){		// 仕方ない
			c_printf("[r]%s\n", "sbr_ratio <= 0");
			assert(parent_ratio_v > 0);
			return parent_ratio_v;
		}
		return sbr_ratio;
	}
	// 横の棒折り過程における、棒を折る比率を計算。親のTSSBから階層的に生成
	double iTHMM::compute_expectation_of_horizontal_htssb_sbr_ratio(Node* target_in_htssb){
		// c_printf("[*]%s\n", "compute_expectation_of_horizontal_htssb_sbr_ratio");
		assert(target_in_htssb != NULL);
		if(target_in_htssb->_depth_v == 0){	// ルートノードなら必ず止まる
			return 1;
		}
		Node* owner_in_structure = target_in_htssb->_owner_in_structure;
		assert(owner_in_structure != NULL);
		int num_itr_in_structure = owner_in_structure->_depth_v + 1;
		int num_itr_horizontal = target_in_htssb->_depth_h + 1;
		double parent_ratio_h = 0;
		double sbr_ratio = 0;
		double* stop_ratio_over_parent = target_in_htssb->_stop_ratio_h_over_parent;
		double* stop_probability_over_parent = target_in_htssb->_stop_probability_h_over_parent;
		for(int n = 0;n < num_itr_in_structure;n++){
			// cout << "n = " << n << endl;
			// トップレベルのノードから順に停止確率を計算
			Node* iterator_in_structure = owner_in_structure->_nodes_from_root_to_myself[n];
			Node* iterator_in_htssb = iterator_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(target_in_htssb);
			assert(iterator_in_htssb != NULL);
			double sum_parent_stop_probability = 0;
			Node* parent_contains_target_in_htssb = iterator_in_htssb->_parent;
			assert(parent_contains_target_in_htssb != NULL);
			// cout << "parent contains" << endl << "	";
			// parent_contains_target_in_htssb->dump();
			for(int m = 0;m < num_itr_horizontal;m++){		// 自分自身も含めるので+1
				Node* child_in_htssb = parent_contains_target_in_htssb->_children[m];
				assert(child_in_htssb);
				// cout << "m = " << m << endl << "	";
				if(n == 0){		// 親ノードの場合
					int pass_count = child_in_htssb->_pass_count_h;
					int stop_count = child_in_htssb->_stop_count_h;
					double gamma = _gamma * pow(_lambda_gamma, std::max(0, iterator_in_htssb->_depth_v - 1));
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double ratio_h = (1.0 + stop_count) / (1.0 + gamma + stop_count + pass_count + EPS);
					// assert(ratio_h < 1);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
					sbr_ratio = ratio_h;
				}else{
					int pass_count = child_in_htssb->_pass_count_h;
					int stop_count = child_in_htssb->_stop_count_h;
					// cout << "pass_count = " << pass_count << ", stop_count = " << stop_count << endl;
					double parent_stop_probability = stop_probability_over_parent[m];
					parent_ratio_h = (stop_ratio_over_parent[m] > 0) ?  stop_ratio_over_parent[m] : parent_ratio_h;
					// cout << "parent_stop_probability = " << parent_stop_probability << endl;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					// ルートノードではガンマを使うがそれ以外は集中度を使う
					double ratio_h = (_strength * parent_stop_probability + stop_count) / (_strength * (1.0 - sum_parent_stop_probability) + stop_count + pass_count + EPS);
					assert(ratio_h >= 0);
					// assert(ratio_h < 1);
					// cout << "ratio_h = " << ratio_h << endl;
					stop_ratio_over_parent[m] = ratio_h;
					sum_parent_stop_probability += parent_stop_probability;
					// cout << "sum_parent_stop_probability = " << sum_parent_stop_probability << endl;
					sbr_ratio = ratio_h;
				}

			}
			if(n < num_itr_in_structure - 1){
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
		if(sbr_ratio <= 0){		// 仕方ない
			c_printf("[r]%s\n", "sbr_ratio_h <= 0");
			assert(parent_ratio_h > 0);
			return parent_ratio_h;
		}
		return sbr_ratio;
	}
	double iTHMM::compute_Pw_given_s(id token_id, Node* node_in_structure){
		assert(node_in_structure != NULL);
		assert(node_in_structure->_hpylm != NULL);
		assert(is_node_in_structure_tssb(node_in_structure));
		assert(_hpylm_d_m.size() > node_in_structure->_depth_v);
		assert(_hpylm_theta_m.size() > node_in_structure->_depth_v);
		assert(_word_g0 > 0);
		return node_in_structure->_hpylm->compute_Pw(token_id, _word_g0, _hpylm_d_m, _hpylm_theta_m);
	}
	// TSSBの全ての棒の長さを計算
	void iTHMM::update_stick_length_of_tssb(TSSB* tssb, double total_stick_length, bool htssb_mode){
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
	void iTHMM::_update_stick_length_of_parent_node(double &sum_probability, Node* parent, bool htssb_mode){
		assert(parent->_children_stick_length > 0);
		double rest_stick_length = parent->_children_stick_length;	// 親ノードが持っている子ノードに割り当てる棒の長さ
		double sum_stick_length_from_left_to_current_node = sum_probability;
		for(int i = 0;i < parent->_children.size();i++){
			Node* child = parent->_children[i];
			double ratio_h = compute_expectation_of_horizontal_sbr_ratio(child, htssb_mode);
			double ratio_v = compute_expectation_of_vertical_sbr_ratio(child, htssb_mode);
			child->_stick_length = rest_stick_length * ratio_h;		// このノードかこのノードの子ノードに止まる確率
			child->_probability = child->_stick_length * ratio_v;	// このノードに止まる確率
			child->_children_stick_length = child->_stick_length * (1.0 - ratio_v);	// 子ノードに割り当てる長さはこのノードに降りない確率
			sum_probability += child->_probability;
			sum_stick_length_from_left_to_current_node += child->_stick_length;
			child->_sum_probability = sum_stick_length_from_left_to_current_node - child->_children_stick_length;				// このノードより左側の全ての棒の長さの総和
			rest_stick_length *= 1.0 - ratio_h;
			if(child->has_child()){
				_update_stick_length_of_parent_node(child->_sum_probability, child, htssb_mode);
			}
		}
	}
	// 不要なノードの削除
	void iTHMM::delete_invalid_children(){
		_delete_invalid_children_in_structure_tssb(_structure_tssb);
	}
	void iTHMM::_delete_invalid_children_in_structure_tssb(TSSB* tssb){
		assert(is_tssb_structure(tssb));
		_delete_invalid_children_of_node_in_structure(tssb->_root);
	}
	void iTHMM::_delete_invalid_children_of_node_in_structure(Node* parent){
		assert(is_node_in_structure_tssb(parent));
		std::vector<Node*> &children = parent->_children;
		for(int i = children.size() - 1;i >= 0;i--){
			Node* child = children[i];
			_delete_invalid_children_of_node_in_structure(child);
			bool success = delete_node_in_structure_if_needed(child);
			if(success == false){	// 失敗したらそれ以上は消さない
				// break;
			}
		}
	}
	bool iTHMM::delete_node_in_structure_if_needed(Node* target_in_structure){
		assert(target_in_structure != NULL);
		assert(is_node_in_structure_tssb(target_in_structure));
		if(target_in_structure->_depth_v == 0){
			return false;
		}
		assert(target_in_structure->_parent != NULL);
		int delete_id = target_in_structure->_identifier;
		if(target_in_structure->_pass_count_v != 0){
			return false;
		}
		if(target_in_structure->_stop_count_v != 0){
			return false;
		}
		if(target_in_structure->_pass_count_h != 0){
			return false;
		}
		if(target_in_structure->_stop_count_h != 0){
			return false;
		}
		if(target_in_structure->_ref_count != 0){
			return false;
		}
		if(target_in_structure->_num_transitions_to_eos != 0){
			return false;
		}
		if(target_in_structure->_num_transitions_to_other != 0){
			return false;
		}
		Node* parent_in_structure = target_in_structure->_parent;
		Node* delete_node = target_in_structure->_parent->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		// 全てのHTSSBから削除
		_delete_node_in_all_htssb(delete_id, _structure_tssb->_root, parent_in_structure);
		// <bos>から削除
		Node* parent_in_bos = _bos_tssb->find_node_by_tracing_horizontal_indices(parent_in_structure);
		assert(parent_in_bos != NULL);
		assert(is_node_in_bos_tssb(parent_in_bos));
		delete_node = parent_in_bos->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		return true;
	}
	void iTHMM::_delete_node_in_all_htssb(int delete_id, Node* iterator_in_structure, Node* target_parent_in_structure){
		assert(target_parent_in_structure != NULL);
		assert(iterator_in_structure->_transition_tssb != NULL);
		// 遷移確率用TSSBでの同じ位置の子ノードを削除
		Node* parent_in_htssb = iterator_in_structure->_transition_tssb->find_node_by_tracing_horizontal_indices(target_parent_in_structure);
		assert(parent_in_htssb != NULL);
		assert(is_node_in_htssb(parent_in_htssb));
		Node* delete_node = parent_in_htssb->delete_child_node(delete_id);
		if(delete_node != NULL){
			TSSB* delete_tssb = delete_node->_transition_tssb;
			delete delete_node;
			delete delete_tssb;
		}
		for(const auto &child: iterator_in_structure->_children){
			_delete_node_in_all_htssb(delete_id, child, target_parent_in_structure);
		}
	}
	// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
	// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
	void iTHMM::sum_auxiliary_variables_recursively_for_hpylm(Node* parent, std::vector<double> &sum_log_x_u_m, std::vector<double> &sum_y_ui_m, std::vector<double> &sum_1_y_ui_m, std::vector<double> &sum_1_z_uwkj_m){
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
	// HPYLMのdとθの推定
	void iTHMM::sample_hpylm_hyperparameters(){
		assert(_current_max_depth < _hpylm_d_m.size());
		assert(_current_max_depth < _hpylm_theta_m.size());
		assert(_current_max_depth < _hpylm_a_m.size());
		assert(_current_max_depth < _hpylm_b_m.size());
		assert(_current_max_depth < _hpylm_alpha_m.size());
		assert(_current_max_depth < _hpylm_beta_m.size());

		// ルートノードの深さが0であることに注意
		std::vector<double> sum_log_x_u_m(_current_max_depth + 1, 0.0);
		std::vector<double> sum_y_ui_m(_current_max_depth + 1, 0.0);
		std::vector<double> sum_1_y_ui_m(_current_max_depth + 1, 0.0);
		std::vector<double> sum_1_z_uwkj_m(_current_max_depth + 1, 0.0);

		// ルートノード
		HPYLM* root = _structure_tssb->_root->_hpylm;
		sum_log_x_u_m[0] = root->auxiliary_log_x_u(_hpylm_theta_m[0]);					// log(x_u)
		sum_y_ui_m[0] = root->auxiliary_y_ui(_hpylm_d_m[0], _hpylm_theta_m[0]);			// y_ui
		sum_1_y_ui_m[0] = root->auxiliary_1_y_ui(_hpylm_d_m[0], _hpylm_theta_m[0]);		// 1 - y_ui
		sum_1_z_uwkj_m[0] = root->auxiliary_1_z_uwkj(_hpylm_d_m[0]);					// 1 - z_uwkj

		// それ以外
		sum_auxiliary_variables_recursively_for_hpylm(_structure_tssb->_root, sum_log_x_u_m, sum_y_ui_m, sum_1_y_ui_m, sum_1_z_uwkj_m);

		// サンプリング
		for(int u = 0;u <= _current_max_depth;u++){
			_hpylm_d_m[u] = sampler::beta(_hpylm_a_m[u] + sum_1_y_ui_m[u], _hpylm_b_m[u] + sum_1_z_uwkj_m[u]);
			_hpylm_theta_m[u] = sampler::gamma(_hpylm_alpha_m[u] + sum_y_ui_m[u], _hpylm_beta_m[u] - sum_log_x_u_m[u]);
		}

		assert(_hpylm_d_m.size() == _hpylm_theta_m.size());
		assert(_hpylm_theta_m.size() == _hpylm_a_m.size());
		assert(_hpylm_a_m.size() == _hpylm_b_m.size());
		assert(_hpylm_b_m.size() == _hpylm_alpha_m.size());
		assert(_hpylm_alpha_m.size() == _hpylm_beta_m.size());

		// 不要な深さのハイパーパラメータを削除
		int num_remove = _hpylm_d_m.size() - _current_max_depth - 1;
		for(int n = 0;n < num_remove;n++){
			_hpylm_d_m.pop_back();
			_hpylm_theta_m.pop_back();
			_hpylm_a_m.pop_back();
			_hpylm_b_m.pop_back();
			_hpylm_alpha_m.pop_back();
			_hpylm_beta_m.pop_back();
		}
	}
	void iTHMM::geneerate_word_ranking_of_node(Node* node_in_structure, std::multiset<std::pair<id, double>, multiset_value_comparator> &ranking){
		assert(node_in_structure != NULL);
		assert(is_node_in_structure_tssb(node_in_structure));
		HPYLM* hpylm = node_in_structure->_hpylm;
		assert(hpylm != NULL);
		assert(_word_g0 > 0);
		std::pair<id, double> pair = std::make_pair(0, 0);
		for(const auto &elem: node_in_structure->_num_word_assignment){
			id word_id = elem.first;
			double Pw = hpylm->compute_Pw(word_id, _word_g0, _hpylm_d_m, _hpylm_theta_m);
			pair.first = word_id;
			pair.second = Pw;
			ranking.insert(pair);
		}
	}
	bool iTHMM::save(std::string filename){
		bool success = false;
		std::ofstream ofs(filename);
		if(ofs.good()){
			boost::archive::binary_oarchive oarchive(ofs);
			oarchive << static_cast<const iTHMM&>(*this);
			success = true;
		}
		ofs.close();
		return success;
	}
	bool iTHMM::load(std::string filename){
		bool success = false;
		std::ifstream ifs(filename);
		if(ifs.good()){
			boost::archive::binary_iarchive iarchive(ifs);
			iarchive >> *this;
			assert(_structure_tssb != NULL);
			assert(_structure_tssb->_root != NULL);
			std::vector<Node*> nodes;
			_structure_tssb->enumerate_nodes_from_left_to_right(nodes);
			for(auto node: nodes){
				// 配列を確保
				node->init_arrays();
				node->init_horizontal_indices();
				node->init_pointers_from_root_to_myself();
				std::vector<Node*> nodes_in_htssb;
				node->_transition_tssb->enumerate_nodes_from_left_to_right(nodes_in_htssb);
				for(auto node_in_htssb: nodes_in_htssb){
					// 配列を確保
					node_in_htssb->init_arrays();
					node_in_htssb->init_horizontal_indices();
					node_in_htssb->init_pointers_from_root_to_myself();
				}
				std::vector<Node*>().swap(nodes_in_htssb);	// 解放
			}
			nodes.clear();
			std::vector<Node*>().swap(nodes);				// 解放
			_bos_tssb->enumerate_nodes_from_left_to_right(nodes);
			for(auto node: nodes){
				// 配列を確保
				node->init_arrays();
				node->init_horizontal_indices();
				node->init_pointers_from_root_to_myself();
			}
			std::vector<Node*>().swap(nodes);				// 解放
			success = true;
		}
		ifs.close();

		return success;
	}
}