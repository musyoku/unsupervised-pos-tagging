# -*- coding: utf-8 -*-
import argparse, time
import model

def main(args):
	if args.filename is None:
		raise Exception()

	hmm = model.bayesian_hmm()
	hmm.set_max_epoch(args.epoch);
	hmm.set_num_tags(args.num_tags);
	hmm.load_textfile(args.filename)

	if hmm.load("hmm.model") == False:
		hmm.initialize()
	hmm.perform_gibbs_sampling()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-n", "--num_tags", type=int, default=30, help="品詞の数.")
	args = parser.parse_args()
	main(args)