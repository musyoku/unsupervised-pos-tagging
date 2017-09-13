#include  <iostream>
#include  <string>
#include "../src/ihmm/utils.h"
#include "../src/ihmm/ihmm.h"
using namespace ihmm;

void test1(){
	InfiniteHMM* ihmm = new InfiniteHMM(10, 100);
	ihmm->_increment_tag_unigram_count(1);
}

int main(){
	test1();
	return 0;
}