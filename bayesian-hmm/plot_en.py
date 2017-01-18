# -*- coding: utf-8 -*-
import argparse, sys, re, pylab
import treetaggerwrapper
import pandas as pd
import seaborn as sns
sns.set()
import model
from train_en import colapse_pos, posset

def main(args):
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	all_types_of_pos = set()	# TreeTaggerが出力する品詞の総数
	num_occurrence_of_pos_for_tag = []
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	threashold = 0	# 出現回数がこの値以下の単語は除外
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
			pos = colapse_pos(pos)
			all_types_of_pos.add(pos)
			if pos not in num_occurrence_of_pos:
				num_occurrence_of_pos[pos] = 1
			else:
				num_occurrence_of_pos[pos] += 1
		num_occurrence_of_pos_for_tag.append(num_occurrence_of_pos)

	# 存在しない部分を0埋め
	for occurrence in num_occurrence_of_pos_for_tag:
		for pos in all_types_of_pos:
			if pos not in occurrence:
				occurrence[pos] = 0

	# タグごとに正規化
	for pos in all_types_of_pos:
		z = 0
		for occurrence in num_occurrence_of_pos_for_tag:
			z += occurrence[pos]
		for occurrence in num_occurrence_of_pos_for_tag:
			occurrence[pos] = float(occurrence[pos]) / float(z)

	fig = pylab.gcf()
	fig.set_size_inches(len(all_types_of_pos), hmm.get_num_tags())
	pylab.clf()
	dataframe = pd.DataFrame(num_occurrence_of_pos_for_tag)
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
	heatmap = ax.get_figure()
	heatmap.savefig("pos.png")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルファイル.")
	args = parser.parse_args()
	main(args)