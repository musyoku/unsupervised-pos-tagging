import argparse, os
import bhmm

def main(args):
	dictionary = bhmm.dictionary()
	dictionary.load(os.path.join(args.working_directory, "bhmm.dict"))
	model = bhmm.model(os.path.join(args.working_directory, "bhmm.model"))
	model.print_typical_words_assigned_to_each_tag(args.num_words_to_show, dictionary);
	model.print_alpha_and_beta()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-cwd", "--working-directory", type=str, default="out", help="ワーキングディレクトリ.")
	parser.add_argument("-n", "--num-words-to-show", type=int, default=20, help="各予測タグに属する単語をいくつ表示するか.")
	main(parser.parse_args())