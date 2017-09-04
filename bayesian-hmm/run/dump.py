# coding: utf-8
import argparse, os
import bhmm

def main(args):
	dictionary = bhmm.dictionary()
	dictionary.load(os.path.join(args.model, "bhmm.dict"))
	model = bhmm.model(os.path.join(args.model, "bhmm.model"))
	model.print_typical_words_of_each_tag(args.num_words_to_show, dictionary);
	model.print_alpha_and_beta()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--num-words-to-show", type=int, default=20, help="各予測タグに属する単語をいくつ表示するか.")
	main(parser.parse_args())