# Python 3のみ対応
import argparse, sys, re, pylab, codecs, os
import MeCab
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import bhmm
from train_ja import printr, printb

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
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = MeCab.Tagger()
		for i, sentence_str in enumerate(f):
			if i % 10 == 0:
				printr("データを集計しています ... {}".format(i + 1))
			word_id_seq = []
			true_tag_seq = []
			sentence_str = sentence_str.strip().encode("utf-8")	# 改行を消す
			m = tagger.parseToNode(sentence_str)	# 形態素解析
			words = []
			while m:
				word = m.surface.decode("utf-8")
				features = m.feature.split(",")
				pos_major = features[0]
				pos = (pos_major + "," + features[1]).decode("utf-8")
				if pos == u"名詞,数":
					word = u"##"		# 数字は全て置き換える

				print(word.encode("utf-8"))
				print(type(word))
				raise Exception()
				word_id_seq.append(dictionary.string_to_word_id(word))
				true_tag_seq.append(pos)

				if dictionary.is_unk(word):	# <unk>は無視
					print(word, "is <unk>")
					pass
				else:
					true_tag_set.add(true_tag)
					if true_tag not in words_of_true_tag:
						words_of_true_tag[true_tag] = set()
					words_of_true_tag[true_tag].add(lowercase)


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
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
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

# def main(args):
# 	hmm = model.bayesian_hmm()
# 	if hmm.load(args.model) == False:
# 		raise Exception("モデルが見つかりません.")

# 	# 訓練データを分かち書きする
# 	print stdout.BOLD + "データを集計しています ..." + stdout.END
# 	num_occurrence_of_pos_for_tag = {}
# 	all_types_of_pos = set()
# 	major_pos_count = set()	# 品詞数（大分類）
# 	with codecs.open(args.filename, "r", "utf-8") as f:
# 		tagger = MeCab.Tagger()
# 		for i, line in enumerate(f):
# 			if i % 500 == 0:
# 				sys.stdout.write("\r{}行目を処理中です ...".format(i))
# 				sys.stdout.flush()
# 			tag_ids = [0, 0]	# <bos>の品詞IDは0. 3-gramなので文脈は2つ
# 			line = re.sub(ur"\n", "", line)	# 開業を消す
# 			string = line.encode("utf-8")
# 			m = tagger.parseToNode(string)
# 			while m:
# 				word = m.surface
# 				features = m.feature.split(",")
# 				major_pos_count.add(features[0].decode("utf-8"))
# 				pos = (features[0] + "," + features[1]).decode("utf-8")
# 				major_pos = features[0].decode("utf-8")
# 				if pos == u"名詞,数":
# 					word = "##"		# 数字は全て置き換える
# 				word_id = hmm.string_to_word_id(word.decode("utf-8"))
# 				if major_pos == "BOS/EOS":
# 					word_id = 1
# 				if args.major:	# 品詞の大分類を使う場合
# 					pos = major_pos
# 				all_types_of_pos.add(pos)
# 				tag_id = hmm.argmax_tag_from_Pt_w(tag_ids[-2], tag_ids[-1], word_id)
# 				tag_ids.append(tag_id)
# 				if tag_id not in num_occurrence_of_pos_for_tag:
# 					num_occurrence_of_pos_for_tag[tag_id] = {}
# 				if pos not in num_occurrence_of_pos_for_tag[tag_id]:
# 					num_occurrence_of_pos_for_tag[tag_id][pos] = 0
	# 			num_occurrence_of_pos_for_tag[tag_id][pos] += 1
	# 			m = m.next

	# # 存在しない部分を0埋め
	# for tag, occurrence in num_occurrence_of_pos_for_tag.items():
	# 	for pos in all_types_of_pos:
	# 		if pos not in occurrence:
	# 			occurrence[pos] = 0

	# # 正解品詞ごとに正規化
	# for pos in all_types_of_pos:
	# 	z = 0
	# 	for tag, occurrence in num_occurrence_of_pos_for_tag.items():
	# 		z += occurrence[pos]
	# 	if z > 0:
	# 		for tag, occurrence in num_occurrence_of_pos_for_tag.items():
	# 			occurrence[pos] = float(occurrence[pos]) / float(z)

	# fig = pylab.gcf()
	# fig.set_size_inches(hmm.get_num_tags(), len(all_types_of_pos))
	# pylab.clf()
	# dataframe = pd.DataFrame(num_occurrence_of_pos_for_tag)
	# ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
	# ax.tick_params(labelsize=20) 
	# plt.yticks(rotation=0)
	# plt.xlabel(u"予測タグ")
	# plt.ylabel(u"正解品詞")
	# heatmap = ax.get_figure()
	# heatmap.subplots_adjust(left=0.25, right=0.9)
	# heatmap.savefig("pos.png")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="学習に使ったテキストファイルのパス.")
	parser.add_argument("-cwd", "--working-directory", type=str, default="out", help="ワーキングディレクトリ.")
	args = parser.parse_args()
	main()