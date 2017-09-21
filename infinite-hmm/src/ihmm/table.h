#include <boost/serialization/serialization.hpp>

// 中華料理店過程のテーブル
namespace ihmm {
	class Table {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive & ar, unsigned int version);
	public:
		std::vector<int> _arrangement;
		int _num_customers;
		int _identifier;
		Table();
		Table(int identifier);
		bool is_empty();
		void add_customer(double concentration_parameter, bool &new_table_generated);
		void remove_customer(bool &empty_table_deleted);
		int get_num_customers();
	};
}