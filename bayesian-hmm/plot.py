# -*- coding: utf-8 -*-
import argparse, sys
import model

def main(args):
	if args.filename is None:
		raise Exception()
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")
		
	hmm.show_typical_words_for_each_tag(20);	# それぞれのタグにつき上位n個の単語を表示

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="hmm.model", help="モデルファイル.")
	args = parser.parse_args()
	main(args)