#ifndef _hpylm_
#define _hpylm_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/map.hpp>
#include <vector>
#include <unordered_map>
#include <map>
using namespace std;
using id = int;

class Node;
class HPYLM{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & archive, unsigned int version)
	{
		static_cast<void>(version);
		archive & _arrangement;
		archive & _num_tables;
		archive & _num_customers;
		archive & _depth;
		archive & _parent;
		archive & _state_node;
	}
	bool add_customer_to_table(id token_id, int table_k, double parent_Pw, vector<double> &d_m, vector<double> &theta_m);
	bool add_customer_to_new_table(id token_id, double parent_Pw, vector<double> &d_m, vector<double> &theta_m);
	bool remove_customer_from_table(id token_id, int table_k, vector<int> &num_customers_at_table);
public:
	map<id, vector<int>> _arrangement;	// 客の配置 vector<int>のk番目の要素がテーブルkの客数を表す
	int _num_tables;					// 総テーブル数
	int _num_customers;					// 客の総数
	int _depth;
	HPYLM* _parent;
	Node* _state_node;
	HPYLM();
	HPYLM(Node* node);
	bool child_exists(id token_id);
	bool need_to_remove_from_parent();
	int get_num_tables_serving_word(id token_id);
	int get_num_customers_eating_word(id token_id);
	HPYLM* find_child_node(id token_id, bool generate_if_not_exist = false);
	bool add_customer(id token_id, double g0, vector<double> &d_m, vector<double> &theta_m);
	bool remove_customer(id token_id);
	double compute_Pw(id token_id, double g0, vector<double> &d_m, vector<double> &theta_m);
	bool remove_from_parent();
	void delete_child_node(id token_id);
	id sample_token(double g0, vector<double> &d_m, vector<double> &theta_m);
	int get_max_depth(int base);
	int get_num_tables();
	int get_num_customers();
	// dとθの推定用
	// "A Bayesian Interpretation of Interpolated Kneser-Ney" Appendix C参照
	// http://www.gatsby.ucl.ac.uk/~ywteh/research/compling/hpylm.pdf
	double auxiliary_log_x_u(double theta_u);
	double auxiliary_y_ui(double d_u, double theta_u);
	double auxiliary_1_y_ui(double d_u, double theta_u);
	double auxiliary_1_z_uwkj(double d_u);
	void dump();
};
#endif