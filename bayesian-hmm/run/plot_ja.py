# -*- coding: utf-8 -*-
import argparse, sys, re, pylab, codecs, matplotlib
import MeCab
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import model
from train_ja import stdout

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=3)

def main(args):
	hmm = model.bayesian_hmm()
	if hmm.load(args.model) == False:
		raise Exception("モデルが見つかりません.")

	# 訓練データを分かち書きする
	print stdout.BOLD + "データを集計しています ..." + stdout.END
	num_occurrence_of_pos_for_tag = {}
	all_types_of_pos = set()
	major_pos_count = set()	# 品詞数（大分類）
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = MeCab.Tagger()
		for i, line in enumerate(f):
			if i % 500 == 0:
				sys.stdout.write("\r{}行目を処理中です ...".format(i))
				sys.stdout.flush()
			tag_ids = [0, 0]	# <bos>の品詞IDは0. 3-gramなので文脈は2つ
			line = re.sub(ur"\n", "", line)	# 開業を消す
			string = line.encode("utf-8")
			m = tagger.parseToNode(string)
			while m:
				word = m.surface
				features = m.feature.split(",")
				major_pos_count.add(features[0].decode("utf-8"))
				pos = (features[0] + "," + features[1]).decode("utf-8")
				major_pos = features[0].decode("utf-8")
				if pos == u"名詞,数":
					word = "##"		# 数字は全て置き換える
				word_id = hmm.string_to_word_id(word.decode("utf-8"))
				if major_pos == "BOS/EOS":
					word_id = 1
				if args.major:	# 品詞の大分類を使う場合
					pos = major_pos
				all_types_of_pos.add(pos)
				tag_id = hmm.argmax_tag_from_Pt_w(tag_ids[-2], tag_ids[-1], word_id)
				tag_ids.append(tag_id)
				if tag_id not in num_occurrence_of_pos_for_tag:
					num_occurrence_of_pos_for_tag[tag_id] = {}
				if pos not in num_occurrence_of_pos_for_tag[tag_id]:
					num_occurrence_of_pos_for_tag[tag_id][pos] = 0
				num_occurrence_of_pos_for_tag[tag_id][pos] += 1
				m = m.next

	# 存在しない部分を0埋め
	for tag, occurrence in num_occurrence_of_pos_for_tag.items():
		for pos in all_types_of_pos:
			if pos not in occurrence:
				occurrence[pos] = 0

	# 正解品詞ごとに正規化
	for pos in all_types_of_pos:
		z = 0
		for tag, occurrence in num_occurrence_of_pos_for_tag.items():
			z += occurrence[pos]
		if z > 0:
			for tag, occurrence in num_occurrence_of_pos_for_tag.items():
				occurrence[pos] = float(occurrence[pos]) / float(z)

	fig = pylab.gcf()
	fig.set_size_inches(hmm.get_num_tags(), len(all_types_of_pos))
	pylab.clf()
	dataframe = pd.DataFrame(num_occurrence_of_pos_for_tag)
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
	ax.tick_params(labelsize=20) 
	plt.yticks(rotation=0)
	plt.xlabel(u"予測タグ")
	plt.ylabel(u"正解品詞")
	heatmap = ax.get_figure()
	heatmap.subplots_adjust(left=0.25, right=0.9)
	heatmap.savefig("pos.png")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="学習に使ったテキストファイルのパス. 分かち書きされていない必要がある.")
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルファイル.")
	parser.add_argument("--major", dest="major", default=False, action="store_true", help="品詞の大分類を使うかどうか.")
	args = parser.parse_args()
	main(args)