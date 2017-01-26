# -*- coding: utf-8 -*-
import argparse, sys, os, time, re, codecs
import model

def main(args):
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	hmm.show_typical_words_for_each_tag(args.num_words_to_show);	# それぞれのタグにつき上位n個の単語を表示

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--num-words-to-show", type=int, default=20, help="各予測タグに属する単語をいくつ表示するか.")
	main(parser.parse_args())