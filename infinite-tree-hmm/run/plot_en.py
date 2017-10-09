import argparse, sys, pylab, codecs, os
import treetaggerwrapper
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ithmm
from train_en import collapse_pos, posset, parse_tagger_result_str

# フォントをセット
# UbuntuならTakaoGothicなどが標準で入っている
sns.set(font=["MS Gothic"], font_scale=3)

def main():
	assert args.train_filename is not None

	# モデル
	model = ithmm.model()
	if model.load(os.path.join(args.model, "ithmm.model")) == False:
		raise Exception("モデルが見つかりません.")

	# 単語辞書
	dictionary = ithmm.dictionary()
	dictionary.load(os.path.join(args.model, "ithmm.dict"))

	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	all_types_of_pos = set()
	state_sequence_array_true = []		# 真の品詞列
	state_sequence_array_viterbi = []	# ビタビアルゴリズムによる推定
	
	def append(f):
		for i, line in enumerate(f):
			line = line.strip()				# 改行を消す
			results = tagger.tag_text(line)	# 形態素解析
			word_ids = []
			state_sequence_true = []
			for result_str in results:
				pos, word = parse_tagger_result_str(result_str)
				all_types_of_pos.add(pos)
				word_ids.append(dictionary.string_to_word_id(word))
				state_sequence_true.append(pos)
			state_sequence_array_true.append(state_sequence_true)
			state_sequence_viterbi = model.viterbi_decode(word_ids)
			assert len(state_sequence_viterbi) == len(state_sequence_true)
			state_sequence_array_viterbi.append(state_sequence_viterbi)

	# 読み込み
	if args.train_filename.endswith(".txt"):
		with codecs.open(args.train_filename, "r", "utf-8") as f:
			append(f)
	else:
		train_dir = args.train_filename
		files = os.listdir(train_dir)
		for filename in files:
			with codecs.open(os.path.join(train_dir, filename), "r", "utf-8") as f:
				append(f)

	# モデルの予測品詞と正解品詞の対応関係
	num_occurrence_of_pos_for_tag = {}
	for state_sequence_true, state_sequence_viterbi in zip(state_sequence_array_true, state_sequence_array_viterbi):
		for pos_true, tag_viterbi in zip(state_sequence_true, state_sequence_viterbi):
			if tag_viterbi not in num_occurrence_of_pos_for_tag:
				num_occurrence_of_pos_for_tag[tag_viterbi] = {}
			if pos_true not in num_occurrence_of_pos_for_tag[tag_viterbi]:
				num_occurrence_of_pos_for_tag[tag_viterbi][pos_true] = 0
			num_occurrence_of_pos_for_tag[tag_viterbi][pos_true] += 1

	tags = model.get_all_states()

	# 存在しない部分を0埋め
	for tag, occurrence in num_occurrence_of_pos_for_tag.items():
		for pos in all_types_of_pos:
			if pos not in occurrence:
				occurrence[pos] = 0
	
	for tag in tags:
		if tag not in num_occurrence_of_pos_for_tag:
			num_occurrence_of_pos_for_tag[tag] = {}
			for pos in all_types_of_pos:
				num_occurrence_of_pos_for_tag[tag][pos] = 0

	# 正解品詞ごとに正規化
	for pos in all_types_of_pos:
		z = 0
		for tag, occurrence in num_occurrence_of_pos_for_tag.items():
			z += occurrence[pos]
		if z > 0:
			for tag, occurrence in num_occurrence_of_pos_for_tag.items():
				occurrence[pos] = float(occurrence[pos]) / float(z)

	fig = pylab.gcf()
	fig.set_size_inches(len(tags) + 3, len(all_types_of_pos))
	pylab.clf()
	dataframe = pd.DataFrame(num_occurrence_of_pos_for_tag)
	ax = sns.heatmap(dataframe, annot=False, fmt="f", linewidths=0)
	ax.tick_params(labelsize=20) 
	plt.yticks(rotation=0)
	plt.xticks(rotation=90)
	plt.xlabel(u"予測タグ")
	plt.ylabel(u"正解品詞")
	heatmap = ax.get_figure()
	heatmap.savefig("{}/pos.png".format(args.model))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-file", "--train-filename", type=str, default=None, help="訓練用のテキストファイルのパス.ディレクトリも可.")
	parser.add_argument("-m", "--model", type=str, default="out", help="モデルの保存ディレクトリ.")
	args = parser.parse_args()
	main()