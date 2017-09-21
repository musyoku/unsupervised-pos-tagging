import argparse, sys, re, pylab, codecs, os
import treetaggerwrapper
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import bhmm
from train_en import collapse_true_tag, printr, printb

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=2, style="ticks")

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

	true_tag_to_id = {}
	for true_tag in true_tag_set:
		true_tag_to_id[true_tag] = len(true_tag_to_id)

	num_true_tags = len(true_tag_to_id)
	confusion_mat = np.zeros((num_true_tags, model.get_num_tags()), dtype=int)
	for tag, occurrence in num_true_tags_of_found_tag.items():
		for true_tag in true_tag_set:
			num = 0
			if true_tag in occurrence:
				num = occurrence[true_tag]
			confusion_mat[true_tag_to_id[true_tag]][tag - 1] = num
	normalized_confusion_mat = confusion_mat / np.sum(confusion_mat, axis=0)

	fig = pylab.gcf()
	fig.set_size_inches(model.get_num_tags() + 1, len(true_tag_set))
	pylab.clf()
	dataframe = pd.DataFrame(normalized_confusion_mat)
	yticks = np.arange(0, num_true_tags)
	yticks_labeles = [None] * len(true_tag_to_id)
	for true_tag, index in true_tag_to_id.items():
		yticks_labeles[index] = true_tag
	ax = sns.heatmap(dataframe, annot=confusion_mat, annot_kws={"size": 14}, fmt="d", linewidths=0, cmap="gray_r", 
		yticklabels=yticks_labeles, xticklabels=np.arange(1, model.get_num_tags() + 1), cbar=False, square=True)
	for _, spine in ax.spines.items():
		spine.set_visible(True)	# 枠線を表示
	ax.tick_params(labelsize=20)
	plt.yticks(rotation=0)
	plt.xlabel("Found Tags", fontname="Arial", fontsize=28, labelpad=20)
	plt.ylabel("True Tags", fontname="Arial", fontsize=28, labelpad=20)
	heatmap = ax.get_figure()
	heatmap.savefig("{}/confusion_mat.png".format(args.working_directory))

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