#include "model.cpp"

int main(){
	PyHpylmHMM* hmm = new PyHpylmHMM(30);
	hmm->load_textfile("../wiki.txt");
	hmm->prepare_for_training();
	double ppl = hmm->compute_perplexity();
	cout << ppl << endl;
	for(int epoch = 0;epoch < 1000;epoch++){
		hmm->perform_gibbs_sampling();
	}
	return 0;
}