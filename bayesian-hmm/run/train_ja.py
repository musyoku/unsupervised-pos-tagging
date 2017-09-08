import argparse, sys, os, time, codecs, random
import MeCab
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

# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
def build_corpus(filename):
	dataset = bhmm.dataset()
	sentence_list = []
	with codecs.open(filename, "r", "utf-8") as f:
		for sentence_str in f:
			sentence_list.append(sentence_str)
	random.shuffle(sentence_list)	# データをシャッフル
	train_split = int(len(sentence_list) * args.train_split)

	word_count = set()	# 単語の種類の総数
	pos_count = set()	# 品詞数
	major_pos_count = set()	# 品詞数（大分類）
	Wt_count = {}

	tagger = MeCab.Tagger()
	for i, sentence_str in enumerate(sentence_list):
		sentence_str = sentence_str.strip()
		if i % 10 == 0:
			printr("データを準備しています ... {}".format(i + 1))
		m = tagger.parseToNode(sentence_str)	# 形態素解析
		words = []
		while m:
			word = m.surface
			features = m.feature.split(",")
			pos_major = features[0]
			pos = (pos_major + "," + features[1])
			major_pos_count.add(pos_major)
			pos_count.add(pos)
			if pos == "名詞,数":
				word = "##"		# 数字は全て置き換える
			words.append(word)
			word_count.add(word)
			if pos_major not in Wt_count:
				Wt_count[pos_major] = {}
			if word not in Wt_count[pos_major]:
				Wt_count[pos_major][word] = 1
			else:
				Wt_count[pos_major][word] += 1
			m = m.next

		if len(words) == 0:
			continue

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
	parser.add_argument("-file", "--train-filename", type=str, default=None, help="訓練用のテキストファイルのパス..")
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