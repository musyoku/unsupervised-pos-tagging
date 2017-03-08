# coding: utf-8
import sys, argparse
import model

def main(args):
	ithmm = model.ithmm()
	if ithmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	if args.words:
		show_words(ithmm)

def show_words(ithmm):
	ithmm.show_assigned_words_and_probability_for_each_tag(30)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("--words", action="store_true", default=False, help="各状態に割り当てられた単語を確率の高い順に表示.")
	main(parser.parse_args())