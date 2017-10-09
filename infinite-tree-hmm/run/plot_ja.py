import argparse, sys, re, pylab, codecs, os, copy
import MeCab
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ithmm
from train_ja import printr, printb

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=2, style="ticks")

def main():
	# 辞書
	dictionary = ithmm.dictionary()
	dictionary.load(os.path.join(args.working_directory, "ithmm.dict"))

	# モデル
	model = ithmm.model(os.path.join(args.working_directory, "ithmm.model"))

	# 訓練データを形態素解析して集計
	found_tags = model.get_tags()
	found_tags_to_id = {}
	num_true_tags_of_found_tag = {}
	for tag in found_tags:
		found_tags_to_id[tag] = len(found_tags_to_id)
		num_true_tags_of_found_tag[tag] = {}
	words_of_true_tag = {}
	true_tag_set = set()
	tagger = MeCab.Tagger()
	tagger.parse("")  # 追加
	with codecs.open(args.filename, "r", "utf-8") as f:
		for i, sentence_str in enumerate(f):
			if i % 10 == 0:
				printr("データを集計しています ... {}".format(i + 1))
			word_id_seq = []
			true_tag_seq = []
			sentence_str = sentence_str.strip()	# 改行を消す
			m = tagger.parseToNode(sentence_str)	# 形態素解析
			words = []
			while m:
				word = m.surface
				features = m.feature.split(",")
				pos_major = features[0]
				pos = (pos_major + "," + features[1])
				if pos == "名詞,数":
					word = "##"		# 数字は全て置き換える

				word_id_seq.append(dictionary.string_to_word_id(word))
				true_tag_seq.append(pos_major)

				if dictionary.is_string_unk(word):	# <unk>は無視
					pass
				else:
					true_tag_set.add(pos_major)
					if pos_major not in words_of_true_tag:
						words_of_true_tag[pos_major] = set()
					words_of_true_tag[pos_major].add(word)

				m = m.next

			found_tag_seq = model.viterbi_decode(word_id_seq)
			for word_id, true_tag, found_tag in zip(word_id_seq, true_tag_seq, found_tag_seq):
				if word_id == 0:	# <unk>は無視
					pass
				else:
					assert found_tag in num_true_tags_of_found_tag
					if true_tag not in num_true_tags_of_found_tag[found_tag]:
						num_true_tags_of_found_tag[found_tag][true_tag] = 0
					num_true_tags_of_found_tag[found_tag][true_tag] += 1

	true_tag_to_id = {}
	for true_tag in true_tag_set:
		true_tag_to_id[true_tag] = len(true_tag_to_id)

	num_true_tags = len(true_tag_to_id)
	confusion_mat = np.zeros((num_true_tags, len(found_tags)), dtype=int)
	for found_tag, occurrence in num_true_tags_of_found_tag.items():
		for true_tag in true_tag_set:
			num = 0
			if true_tag in occurrence:
				num = occurrence[true_tag]
			confusion_mat[true_tag_to_id[true_tag]][found_tags_to_id[found_tag]] = num
	sum_counts = np.sum(confusion_mat, axis=0)
	delete_indices = []
	for index, count in enumerate(sum_counts):
		if count == 0:
			delete_indices.append(index)
	num_nonzero_found_tags = len(found_tags) - len(delete_indices)
	if num_nonzero_found_tags > 0:
		new_confusion_mat = np.zeros((num_true_tags, num_nonzero_found_tags), dtype=int)
		new_found_tags = []
		pointer = 0
		for tag_index, found_tag in enumerate(found_tags):
			if tag_index in delete_indices:
				continue
			new_confusion_mat[:, pointer] = confusion_mat[:, tag_index]
			new_found_tags.append(found_tag)
			pointer += 1
		confusion_mat = new_confusion_mat
		found_tags = new_found_tags
	normalized_confusion_mat = confusion_mat / np.sum(confusion_mat, axis=0)

	fig = pylab.gcf()
	fig.set_size_inches(len(found_tags) + 1, len(true_tag_set))
	pylab.clf()
	dataframe = pd.DataFrame(normalized_confusion_mat)
	yticks = np.arange(0, num_true_tags)
	yticks_labeles = [None] * len(true_tag_to_id)
	for true_tag, index in true_tag_to_id.items():
		yticks_labeles[index] = true_tag
	ax = sns.heatmap(dataframe, annot=confusion_mat, annot_kws={"size": 14}, fmt="d", linewidths=0, cmap="gray_r", 
		yticklabels=yticks_labeles, xticklabels=found_tags, cbar=False, square=True)
	for _, spine in ax.spines.items():
		spine.set_visible(True)	# 枠線を表示
	ax.tick_params(labelsize=20)
	plt.yticks(rotation=0)
	plt.xticks(rotation=90)
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