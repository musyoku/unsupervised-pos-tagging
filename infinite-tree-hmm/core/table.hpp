#ifndef _table_
#define _table_
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>
using namespace std;

// 中華料理店過程のテーブル
// 通常CRPではテーブルが各クラスタを表すが、TSSBでは全テーブルが同じクラスタに属する

class Table{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive & archive, unsigned int version)
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
	Table();
	Table(int token_id);
	bool is_empty();
	void add_customer(double concentration_parameter, double g0, int num_total_customers, bool &new_table_generated);
	void remove_customer(bool &empty_table_deleted);
};
#endif