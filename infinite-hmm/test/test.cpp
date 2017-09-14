#include  <iostream>
#include  <string>
#include "../src/ihmm/utils.h"
#include "../src/ihmm/ihmm.h"
using namespace ihmm;
using namespace std;

void test1(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	int num_tags = ihmm->get_num_tags();
	for(int tag = num_tags;tag >= 1;tag--){
		ihmm->_delete_tag(tag);
	}
	for(int i = 0;i < 100;i++){
		ihmm->_add_new_tag();
		ihmm->_delete_tag(1);
	}
	cout << ihmm->get_num_tags() << endl;
	for(int i = 0;i < 10;i++){
		ihmm->_add_new_tag();
	}
	num_tags = ihmm->get_num_tags();
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_increment_tag_bigram_count(context_tag, tag);
			}
		}
	}
	for(int context_tag = 1;context_tag <= num_tags;context_tag++){
		for(int tag = 1;tag <= num_tags + 1;tag++){
			double prob = ihmm->compute_p_tag_given_context(tag, context_tag);
			cout << "p(" << tag << "|" << context_tag << ") = " << prob << endl;
		}
	}
	num_tags = ihmm->get_num_tags();
	for(int tag = 0;tag <= num_tags;tag++){
		cout << "tag: " << tag << endl;
		cout << ihmm->_sum_n_i_over_j[tag] << endl;
		cout << ihmm->_oracle_n_j_counts[tag] << endl;
	}
	cout << "oracle: " << ihmm->_oracle_sum_n_over_j << endl;
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_decrement_tag_bigram_count(context_tag, tag);
			}
		}
	}
	num_tags = ihmm->get_num_tags();
	for(int tag = 0;tag <= num_tags;tag++){
		cout << "tag: " << tag << endl;
		cout << ihmm->_sum_n_i_over_j[tag] << endl;
		cout << ihmm->_oracle_n_j_counts[tag] << endl;
	}
	for(int tag = num_tags;tag >= 1;tag--){
		ihmm->_delete_tag(tag);
	}
	cout << ihmm->_sum_n_i_over_j.size() << endl;
	cout << ihmm->_oracle_n_j_counts.size() << endl;
	cout << ihmm->_sum_word_count_of_tag.size() << endl;
	cout << "oracle: " << ihmm->_oracle_sum_n_over_j << endl;
	delete ihmm;
}

int main(){
	test1();
	return 0;
}