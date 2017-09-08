import argparse, sys, re, pylab, codecs, os
import treetaggerwrapper
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import bhmm
from train_en import collapse_true_tag, printr, printb

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=3)

def main():
	# 辞書
	dictionary = bhmm.dictionary()
	dictionary.load(os.path.join(args.working_directory, "bhmm.dict"))

	# モデル
	model = bhmm.model(os.path.join(args.working_directory, "bhmm.model"))

	# 訓練データを形態素解析して集計
	num_true_tags_of_found_tag = {}
	words_of_true_tag = {}
	true_tag_set = set()
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	with codecs.open(args.filename, "r", "utf-8") as f:
		for i, sentence_str in enumerate(f):
			printr("データを集計しています ... {}".format(i + 1))
			word_id_seq = []
			true_tag_seq = []
			sentence_str = sentence_str.strip()	# 開業を消す
			results = tagger.tag_text(sentence_str)	# 形態素解析
			for metadata in results:
				metadata = metadata.split("\t")
				if len(metadata) == 3:
					word, true_tag, lowercase = metadata
					true_tag = collapse_true_tag(true_tag, lowercase)
				else:
					lowercase = metadata[0]
					true_tag = "X"
				word_id_seq.append(dictionary.string_to_word_id(lowercase))
				true_tag_seq.append(true_tag)

				if dictionary.is_unk(lowercase):	# <unk>は無視
					pass
				else:
					true_tag_set.add(true_tag)
					if true_tag not in words_of_true_tag:
						words_of_true_tag[true_tag] = set()
					words_of_true_tag[true_tag].add(lowercase)

			found_tag_seq = model.viterbi_decode(word_id_seq)
			for word_id, true_tag, found_tag in zip(word_id_seq, true_tag_seq, found_tag_seq):
				if word_id == 0:	# <unk>は無視
					pass
				else:
					if found_tag not in num_true_tags_of_found_tag:
						num_true_tags_of_found_tag[found_tag] = {}
					if true_tag not in num_true_tags_of_found_tag[found_tag]:
						num_true_tags_of_found_tag[found_tag][true_tag] = 0
					num_true_tags_of_found_tag[found_tag][true_tag] += 1

	# 存在しない部分を0埋め
	for tag, occurrence in num_true_tags_of_found_tag.items():
		for true_tag in true_tag_set:
			if true_tag not in occurrence:
				occurrence[true_tag] = 0
	for tag in range(1, model.get_num_tags() + 1):
		if tag not in num_true_tags_of_found_tag:
			num_true_tags_of_found_tag[tag] = {}
			for true_tag in true_tag_set:
				num_true_tags_of_found_tag[tag][true_tag] = 0

	# 予測品詞ごとに正規化
	for tag, occurrence in num_true_tags_of_found_tag.items():
		z = 0
		for true_tag in true_tag_set:
			z += occurrence[true_tag]
		if z > 0:
			for true_tag in true_tag_set:
				occurrence[true_tag] = occurrence[true_tag] / z

	fig = pylab.gcf()
	fig.set_size_inches(model.get_num_tags() + 3, len(true_tag_set))
	pylab.clf()
	dataframe = pd.DataFrame(num_true_tags_of_found_tag)
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0, cmap="gray_r")
	ax.tick_params(labelsize=20) 
	plt.yticks(rotation=0)
	plt.xlabel(u"予測タグ")
	plt.ylabel(u"正解品詞")
	heatmap = ax.get_figure()
	heatmap.savefig("pos.png")

	printr("")
	for true_tag in words_of_true_tag:
		printb("[{}]".format(true_tag))
		words = []
		sys.stdout.write("	")
		for i, word in enumerate(words_of_true_tag[true_tag]):
			sys.stdout.write(word)
			sys.stdout.write(" ")
			if i >= 100:
				break
		sys.stdout.write("\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="学習に使ったテキストファイルのパス.")
	parser.add_argument("-cwd", "--working-directory", type=str, default="out", help="ワーキングディレクトリ.")
	args = parser.parse_args()
	main()