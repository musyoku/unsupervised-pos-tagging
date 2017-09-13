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
		for(int tag = 1;tag <= num_tags;tag++){
			ihmm->_increment_tag_unigram_count(tag);
		}
	}
	for(int i = 0;i < 100;i++){
		for(int tag = 1;tag <= num_tags;tag++){
			ihmm->_decrement_tag_unigram_count(tag);
		}
	}
	cout << ihmm->get_num_tags() << endl;
	for(int i = 0;i < 10;i++){
		ihmm->_add_new_tag();
	}
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_increment_tag_bigram_count(context_tag, tag);
			}
		}
	}
	for(int i = 0;i < 100;i++){
		for(int context_tag = 1;context_tag <= num_tags;context_tag++){
			for(int tag = 1;tag <= num_tags;tag++){
				ihmm->_decrement_tag_bigram_count(context_tag, tag);
			}
		}
	}
	num_tags = ihmm->get_num_tags();
	for(int i = 0;i < 100;i++){
		for(int tag = 1;tag <= num_tags;tag++){
			ihmm->_increment_tag_unigram_count(tag);
		}
	}
	for(int i = 0;i < 100;i++){
		for(int tag = 1;tag <= num_tags;tag++){
			ihmm->_decrement_tag_unigram_count(tag);
		}
	}
	delete ihmm;
}

int main(){
	test1();
	return 0;
}