#include "model.cpp"

int main(){
	PyHpylmHMM* hmm = new PyHpylmHMM(30);
	hmm->load_textfile("../test.txt");
	hmm->prepare_for_training();
	for(int epoch = 0;epoch < 1000;epoch++){
		hmm->perform_gibbs_sampling();
	}
	return 0;
}