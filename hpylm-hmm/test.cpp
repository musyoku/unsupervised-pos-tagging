#include "model.cpp"

int main(){
	PyHpylmHMM* hmm = new PyHpylmHMM(30);
	hmm->load_textfile("../wiki.txt");
	hmm->prepare_for_training();
	hmm->perform_gibbs_sampling();
	return 0;
}