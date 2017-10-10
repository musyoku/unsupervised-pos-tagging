import argparse, sys, os, time, codecs, random
import treetaggerwrapper
import ihmm

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
	corpus = ihmm.corpus()
	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
	sentence_list = []
	with codecs.open(filename, "r", "utf-8") as f:
		for sentence_str in f:
			sentence_list.append(sentence_str)
	with codecs.open(filename, "r", "utf-8") as f:
		tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
		for i, sentence_str in enumerate(f):
			sentence_str = sentence_str.strip()
			if (i + 1) % 10 == 0:
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
				else:
					lowercase = metadata[0]
				words.append(lowercase)
			# データを追加
			corpus.add_words(words)

	return corpus

def main():
	assert args.train_filename is not None
	try:
		os.mkdir(args.working_directory)
	except:
		pass

	# 訓練データを追加
	corpus = build_corpus(args.train_filename)
	dataset = ihmm.dataset(corpus, args.train_split, args.unknown_threshold, args.seed)	# 低頻度語を全て<unk>に置き換える

	# 単語辞書を保存
	dictionary = dataset.get_dict()
	dictionary.save(os.path.join(args.working_directory, "ihmm.dict"))

	# モデル
	model = ihmm.model(args.initial_num_tags, dataset)

	# ハイパーパラメータの設定
	model.set_initial_alpha(args.initial_alpha)
	model.set_initial_beta(args.initial_beta)
	model.set_initial_gamma(args.initial_gamma)
	model.set_initial_gamma_emission(args.initial_gamma_emission)
	model.set_initial_beta_emission(args.initial_beta_emission)

	# 学習の準備
	trainer = ihmm.trainer(dataset, model)

	# 学習ループ
	for epoch in range(1, args.epochs + 1):
		start = time.time()
		trainer.gibbs()	# 新しい状態系列をギブスサンプリング

		# ログ
		elapsed_time = time.time() - start
		printr("Iteration {} / {} - {:.3f} sec".format(epoch, args.epochs, elapsed_time))
		if epoch % 1000 == 0:
			printr("")
			model.print_typical_words_assigned_to_each_tag(20, dictionary)
		if epoch % 100 == 0:
			printr("")
			print("log_likelihood: train {} - dev {}".format(trainer.compute_log_p_dataset_train(), trainer.compute_log_p_dataset_dev()))
			model.save(os.path.join(args.working_directory, "ihmm.model"))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--train-filename", "-file", type=str, default=None, help="訓練用のテキストファイルのパス..")
	parser.add_argument("--seed", type=int, default=1)
	parser.add_argument("--epochs", "-e", type=int, default=100000, help="総epoch.")
	parser.add_argument("--working-directory", "-cwd", type=str, default="out", help="ワーキングディレクトリ.")
	parser.add_argument("--initial-num-tags", "-tags", type=int, default=20, help="タグの種類（semi_supervisedがFalseの時のみ有効）.")
	parser.add_argument("--unknown-threshold", "-unk",  type=int, default=1, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("--train-split", "-split", type=float, default=0.9, help="テキストデータの何割を訓練データにするか.")
	parser.add_argument("--initial-alpha", "-alpha", type=int, default=1, help="alphaの初期値.")
	parser.add_argument("--initial-beta", "-beta", type=int, default=1, help="betaの初期値.")
	parser.add_argument("--initial-gamma", "-gamma", type=int, default=1, help="gammaの初期値.")
	parser.add_argument("--initial-gamma-emission", "-egamma", type=int, default=1, help="gamma_emissionの初期値.")
	parser.add_argument("--initial-beta-emission", "-ebeta", type=int, default=1, help="beta_emissionの初期値.")
	args = parser.parse_args()
	main()