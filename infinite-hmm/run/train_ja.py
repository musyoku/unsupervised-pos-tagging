import argparse, sys, os, time, codecs, random
import MeCab
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

# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
def build_corpus(filename):
	corpus = ihmm.corpus()
	sentence_list = []
	with codecs.open(filename, "r", "utf-8") as f:
		for sentence_str in f:
			sentence_list.append(sentence_str)
	tagger = MeCab.Tagger()
	tagger.parse("")
	for i, sentence_str in enumerate(sentence_list):
		sentence_str = sentence_str.strip()
		if (i + 1) % 10 == 0:
			printr("データを準備しています ... {}".format(i + 1))
		m = tagger.parseToNode(sentence_str)	# 形態素解析
		words = []
		while m:
			word = m.surface
			features = m.feature.split(",")
			pos_major = features[0]
			pos = (pos_major + "," + features[1])
			if pos == "名詞,数":
				word = "##"		# 数字は全て置き換える
			words.append(word)
			m = m.next

		if len(words) == 0:
			continue

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
		trainer.perform_gibbs_sampling()	# 新しい状態系列をギブスサンプリング

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