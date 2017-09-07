# coding: utf-8
from __future__ import print_function
from __future__ import division
import argparse, sys, os, time, codecs, random
import treetaggerwrapper
import bhmm

class stdout:
	BOLD = "\033[1m"
	END = "\033[0m"
	CLEAR = "\033[2K"

def printb(string):
	print(stdout.BOLD + string + stdout.END)

def printr(string):
	sys.stdout.write("\r" + stdout.CLEAR)
	sys.stdout.write(string)
	sys.stdout.flush()

# https://spacy.io/docs/usage/pos-tagging
# https://courses.washington.edu/hypertxt/csar-v02/penntable.html
class POS:
	PUNCT = {":", ",", "'", "\"", "HYPH", "LS", "NFP", "(", ")"}	# punctuation mark
	SYM = {"SYM", "SENT", "#", "$"}	# symbol
	NUM = {"CD"}	# number
	X = {"ADD", "FW", "GW", "XX"}	
	ADJ = {"AFX", "JJ", "JJR", "JJS", "PDT", "PRP$", "WDT", "WP$"}	# adjective
	VERB = {"BES", "HVS", "MD", "VB", "VBD", "VBG", "VBN", "VBZ", "VBP", "VH", "VH", "VHG", "VHN", 
		"VHZ", "VHP", "VHD", "VD", "VDD", "VDG", "VDN", "VDZ", "VDP", "VV", "VV", "VVG", "VVN", "VVZ", "VVP", "VVD"}
	CONJ = {"CC"}	# conjunction
	DET = {"DT"}	# determiner
	ADV = {"EX", "RB", "RBR", "RBS", "WRB"}	# adverb
	ADP = {"IN"}	# conjunction, subordinating or preposition
	NOUN = {"NN", "NNS", "WP", "NP", "NPS"}	# noun
	PROPN = {"NNP", "NNPS"}	# pronoun
	PART = {"POS", "RP", "TO"}
	INTJ = {"UH"}	# interjection

# 品詞をまとめる
def collapse_true_tag(tag, word):
	if word == "##":
		return "NUM"
	if tag in POS.PUNCT:
		return "PUNCT"
	if tag in POS.SYM:
		return "SYM"
	if tag in POS.X:
		return "X"
	if tag in POS.ADJ:
		return "ADJ"
	if tag in POS.VERB:
		return "VERB"
	if tag in POS.CONJ:
		return "CONJ"
	if tag in POS.DET:
		return "DET"
	if tag in POS.ADV:
		return "ADV"
	if tag in POS.ADP:
		return "ADP"
	if tag in POS.NOUN:
		return "NOUN"
	if tag in POS.PROPN:
		return "PROPN"
	if tag in POS.PART:
		return "PART"
	if tag in POS.INTJ:
		return "INTJ"
	if tag in POS.NUM:
		return "NUM"
	return tag

def build_corpus(filename):
	dataset = bhmm.dataset()
	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
	sentence_list = []
	with codecs.open(filename, "r", "utf-8") as f:
		for sentence_str in f:
			sentence_list.append(sentence_str)
	random.shuffle(sentence_list)	# データをシャッフル
	train_split = int(len(sentence_list) * args.train_split)

	word_count = set()	# 単語の種類の総数
	Wt_count = {}
	with codecs.open(filename, "r", "utf-8") as f:
		tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
		for i, sentence_str in enumerate(f):
			sentence_str = sentence_str.strip()
			if i % 10 == 0:
				printr("データを準備しています ... {}".format(i + 1))
			result = tagger.tag_text(sentence_str)
			if len(result) == 0:
				continue
			# 形態素解析を行いながら訓練データも作る
			# 英語は通常スペース区切りなので不要と思うかもしれないが、TreeTaggerを使うと$600が$ 600に分割されたりする
			# そのためplot_en.pyで評価の際に文の単語数が[スペース区切り]と[TreeTagger]で異なる場合があり正しく評価を行えなくなる
			# よって単語分割は全てTreeTaggerによるものに統一しておく
			words = []
			for metadata in result:
				metadata = metadata.split("\t")
				if len(metadata) == 3:
					word, true_tag, lowercase = metadata
					true_tag = collapse_true_tag(true_tag, lowercase)
					if true_tag not in Wt_count:
						Wt_count[true_tag] = {}
					if lowercase not in Wt_count[true_tag]:
						Wt_count[true_tag][lowercase] = 1
					else:
						Wt_count[true_tag][lowercase] += 1
				else:
					lowercase = metadata[0]
				word_count.add(lowercase)
				words.append(lowercase)
			# データを追加
			if i > train_split:
				dataset.add_words_dev(words)		# 評価用データに追加
			else:
				dataset.add_words_train(words)		# 学習用データに追加

	if args.supervised:
		# Wtは各品詞について、その品詞になりうる単語の数が入っている
		Wt = [len(words) for tag, words in Wt_count.items()]
	else:
		# Wtに制限をかけない場合
		Wt = [int(len(word_count) / args.num_tags)] * args.num_tags

	return dataset, Wt

def main():
	assert args.train_filename is not None
	try:
		os.mkdir(args.working_directory)
	except:
		pass

	# 訓練データを追加
	dataset, Wt = build_corpus(args.train_filename)
	dataset.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える

	# 単語辞書を保存
	dictionary = dataset.get_dict()
	dictionary.save(os.path.join(args.working_directory, "bhmm.dict"))

	# モデル
	num_tags = len(Wt) if args.supervised else args.num_tags
	model = bhmm.model(num_tags, dataset, Wt)

	# ハイパーパラメータの設定
	model.set_temperature(args.start_temperature)		# 温度の初期設定
	model.set_minimum_temperature(args.min_temperature)	# 温度の下限
	model.set_initial_alpha(args.initial_alpha)
	model.set_initial_beta(args.initial_beta)

	# 学習の準備
	trainer = bhmm.trainer(dataset, model)

	# 学習ループ
	decay = (args.start_temperature - args.min_temperature) / args.epochs 
	for epoch in range(1, args.epochs + 1):
		start = time.time()
		trainer.perform_gibbs_sampling()	# 新しい状態系列をギブスサンプリング
		trainer.anneal_temperature(decay)	# 温度を下げる

		# ログ
		elapsed_time = time.time() - start
		printr("Iteration {} / {} - temp {:.3f} - {:.3f} sec".format(epoch, args.epochs, model.get_temperature(), elapsed_time))
		if epoch % 1000 == 0:
			printr("")
			model.print_typical_words_of_each_tag(20, dictionary)
		if epoch % 100 == 0:
			printr("ハイパーパラメータのサンプリング ...")
			trainer.update_hyperparameters()	# ハイパーパラメータをサンプリング
			printr("")
			print("log_likelihood: train {} - dev {}".format(trainer.compute_log_p_dataset_train(), trainer.compute_log_p_dataset_dev()))
			model.save(os.path.join(args.working_directory, "bhmm.model"))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-file", "--train-filename", type=str, default=None, help="訓練用のテキストファイルのパス.ディレクトリも可.")
	parser.add_argument("-epochs", "--epochs", type=int, default=100000, help="総epoch.")
	parser.add_argument("-cwd", "--working-directory", type=str, default="out", help="ワーキングディレクトリ.")
	parser.add_argument("--supervised", dest="supervised", default=False, action="store_true", help="各タグのWtを訓練データで制限するかどうか.指定した場合num_tagsは無視される.")
	parser.add_argument("--unsupervised", dest="supervised", action="store_false", help="各タグのWtを訓練データで制限するかどうか.")
	parser.add_argument("-tags", "--num-tags", type=int, default=20, help="タグの種類（semi_supervisedがFalseの時のみ有効）.")
	parser.add_argument("-unk", "--unknown-threshold", type=int, default=1, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("-split", "--train-split", type=float, default=0.9, help="テキストデータの何割を訓練データにするか.")
	parser.add_argument("--start-temperature", type=float, default=1.5, help="開始温度.")
	parser.add_argument("--min-temperature", type=float, default=0.08, help="最小温度.")
	parser.add_argument("--initial-alpha", "-alpha", type=float, default=0.003, help="alphaの初期値.")
	parser.add_argument("--initial-beta", "-beta", type=float, default=1.0, help="betaの初期値.")
	args = parser.parse_args()
	main()