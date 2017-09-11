
namespace ihmm {
	class Table {
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version);
	public:
		vector<int> _arrangement;
		int _num_customers;
		int _token_id;
		Table();
		Table(int token_id);
		bool is_empty();
		void add_customer(double concentration_parameter, bool &new_table_generated);
		void remove_customer(bool &empty_table_deleted);
	};
}