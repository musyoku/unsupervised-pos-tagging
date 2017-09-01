# coding: utf-8
from __future__ import print_function
import argparse, sys, re, pylab, codecs, os
import treetaggerwrapper
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import bhmm
from train_en import collapse_pos, posset, printr, printb

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=3)

def main(args):
	# 辞書
	dictionary = bhmm.dictionary()
	dictionary.save(os.path.join(args.model, "bhmm.dict"))

	# モデル
	model = bhmm.model(os.path.join(args.model, "bhmm.model"))

	# 訓練データを形態素解析して集計
	printb("データを集計しています ...")
	num_true_tags_for_found_tag = {}
	all_types_of_pos = set()
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	with codecs.open(args.filename, "r", "utf-8") as f:
		for i, line in enumerate(f):
			if i % 500 == 0:
				printr("{}行目を処理中です ...".format(i))
			word_id_seq = []
			true_tag_seq = []
			line = re.sub(ur"\n", "", line)	# 開業を消す
			poses = tagger.tag_text(line)	# 形態素解析
			for i, word_pos_lowercase in enumerate(poses):
				pos = collapse_pos(word_pos_lowercase.split("\t")[1])
				lowercase = collapse_pos(word_pos_lowercase.split("\t")[2])
				all_types_of_pos.add(pos)
				true_tag_seq.append(pos)
				word_id_seq.append(dictionary.string_to_word_id(lowercase))

			found_tag_seq = model.viterbi_decode(word_id_seq)
			for true_tag, found_tag in zip(true_tag_seq, found_tag_seq):
				if found_tag not in num_true_tags_for_found_tag:
					num_true_tags_for_found_tag[found_tag] = {}
				if true_tag not in num_true_tags_for_found_tag[found_tag]:
					num_true_tags_for_found_tag[found_tag][true_tag] = 0
				num_true_tags_for_found_tag[found_tag][true_tag] += 1

	# 存在しない部分を0埋め
	for tag, occurrence in num_true_tags_for_found_tag.items():
		for pos in all_types_of_pos:
			if pos not in occurrence:
				occurrence[pos] = 0
	for tag in xrange(hmm.get_num_tags()):
		if tag not in num_true_tags_for_found_tag:
			num_true_tags_for_found_tag[tag] = {}
			for pos in all_types_of_pos:
				num_true_tags_for_found_tag[tag][pos] = 0

	# 正解品詞ごとに正規化
	for pos in all_types_of_pos:
		z = 0
		for tag, occurrence in num_true_tags_for_found_tag.items():
			z += occurrence[pos]
		if z > 0:
			for tag, occurrence in num_true_tags_for_found_tag.items():
				occurrence[pos] = float(occurrence[pos]) / float(z)

	fig = pylab.gcf()
	fig.set_size_inches(hmm.get_num_tags() + 3, len(all_types_of_pos))
	pylab.clf()
	dataframe = pd.DataFrame(num_true_tags_for_found_tag)
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
	ax.tick_params(labelsize=20) 
	plt.yticks(rotation=0)
	plt.xlabel(u"予測タグ")
	plt.ylabel(u"正解品詞")
	heatmap = ax.get_figure()
	heatmap.savefig("pos.png")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="学習に使ったテキストファイルのパス.")
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルファイル.")
	args = parser.parse_args()
	main(args)