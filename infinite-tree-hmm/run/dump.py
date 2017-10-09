import argparse, os
import ithmm

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-cwd", "--working-directory", type=str, default="out", help="ワーキングディレクトリ.")
	parser.add_argument("-n", "--num-words-to-show", type=int, default=20, help="各予測タグに属する単語をいくつ表示するか.")
	args = parser.parse_args()

	dictionary = ithmm.dictionary()
	dictionary.load(os.path.join(args.working_directory, "ithmm.dict"))
	model = ithmm.model(os.path.join(args.working_directory, "ithmm.model"))

	model.show_sticks()
	model.show_assigned_words_for_each_tag(dictionary, 30, False)
	model.show_assigned_words_and_probability_for_each_tag(dictionary, 30)
	# model.show_hpylm_for_each_tag(dictionary)

if __name__ == "__main__":
	main()