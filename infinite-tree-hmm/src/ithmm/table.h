#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>

// 中華料理店過程のテーブル
// 通常CRPではテーブルが各クラスタを表すが、TSSBでは全テーブルが同じクラスタに属する

namespace ithmm {
	class Table {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, unsigned int version);
	public:
		std::vector<int> _arrangement;
		int _num_customers;
		int _token_id;
		int _last_added_index;
		Table();
		Table(int token_id);
		bool is_empty();
		int get_num_customers();
		int get_num_tables();
		void add_customer(double concentration_parameter, double g0, int num_total_customers, bool &new_table_generated);
		void remove_customer(bool &empty_table_deleted);
		void remove_last_customer(bool &empty_table_deleted);
	};
}