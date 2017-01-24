# -*- coding: utf-8 -*-
import argparse, sys, re, pylab
import MeCab
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()
import model

def main(args):
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	# 訓練データを分かち書きする
	print stdout.BOLD + "データを集計しています ..." + stdout.END
	num_occurrence_of_pos_for_tag = []
	word_count = set()	# 単語の種類の総数
	pos_count = set()	# 品詞数
	major_pos_count = set()	# 品詞数（大分類）
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = MeCab.Tagger()
		for i, line in enumerate(f):
			if i % 100 == 0:
				sys.stdout.write("\r{}行目を処理中です ...".format(i))
				sys.stdout.flush()
			segmentation = ""
			line = re.sub(ur"\n", "", line)	# 開業を消す
			string = line.encode("utf-8")
			m = tagger.parseToNode(string)
			while m:
				word = m.surface
				word_count.add(word)
				features = m.feature.split(",")
				major_pos_count.add(features[0].decode("utf-8"))
				pos = (features[0] + "," + features[1]).decode("utf-8")
				pos_count.add(pos)
				if pos == u"名詞,数":
					word = "##"		# 数字は全て置き換える
				segmentation += word + " "
				m = m.next
			segmentation = re.sub(ur" +$", "",  segmentation)	# 行末の空白を除去
			segmentation = re.sub(ur"^ +", "",  segmentation)	# 行頭の空白を除去
			hmm.add_line(segmentation.decode("utf-8"))	# 学習用データに追加











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
	fig.set_size_inches(hmm.get_num_tags(), len(all_types_of_pos))
	pylab.clf()
	dataframe = pd.DataFrame(num_occurrence_of_pos_for_tag)
	ax = sns.heatmap(dataframe.T, annot=False, fmt="f", linewidths=0)
	plt.yticks(rotation=0)
	heatmap = ax.get_figure()
	heatmap.savefig("pos.png")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="学習に使ったテキストファイルのパス. 分かち書きされていない必要がある.")
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルファイル.")
	args = parser.parse_args()
	main(args)