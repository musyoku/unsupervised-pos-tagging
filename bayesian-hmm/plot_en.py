# -*- coding: utf-8 -*-
import argparse, sys, re
import treetaggerwrapper
import model

def main(args):
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	all_types_of_pos = set()	# TreeTaggerが出力する品詞の総数
	num_occurrence_of_pos_for_tag = []
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	threashold = 1	# 出現回数がこの値以下の単語は除外
	tags = hmm.get_all_words_for_each_tag(threashold)
	for tag, words in enumerate(tags):
		sys.stdout.write("\r集計しています ... {} / {}".format(tag, len(tags)))
		sys.stdout.flush()
		num_occurrence_of_pos = {}
		for word, count in words:
			poses = tagger.tag_text(word)
			if len(poses) == 0:
				continue
			poses = poses[0].split("\t")
			if len(poses) != 3:
				continue
			pos = poses[1]
			all_types_of_pos.add(pos)
			if pos not in num_occurrence_of_pos:
				num_occurrence_of_pos[pos] = 1
			else:
				num_occurrence_of_pos[pos] += 1
		num_occurrence_of_pos_for_tag.append(num_occurrence_of_pos)

	print "\n英語の品詞数: {}".format(len(all_types_of_pos))
	print all_types_of_pos

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルファイル.")
	args = parser.parse_args()
	main(args)