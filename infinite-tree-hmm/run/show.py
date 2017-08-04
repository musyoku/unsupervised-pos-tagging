# coding: utf-8
import sys, argparse, os
import ithmm

def main(args):
	model = ithmm.model()
	if model.load(os.path.join(args.model, "ithmm.model")) == False:
		raise Exception("モデルが見つかりません.")

	dictionary = ithmm.dictionary()
	dictionary.load(os.path.join(args.model, "ithmm.dict"))

	if args.words:
		model.show_assigned_words_and_probability_for_each_tag(dictionary, 30)

	if args.tssb:
		model.show_assigned_words_for_each_tag(dictionary, 30, False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("--words", action="store_true", default=False, help="各状態に割り当てられた単語を確率の高い順に表示.")
	parser.add_argument("--tssb", action="store_true", default=False, help="TSSBと各状態に割り当てられた単語を表示.")
	main(parser.parse_args())