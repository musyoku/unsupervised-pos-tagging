# coding: utf-8
from __future__ import print_function
import argparse, sys, os, time, re, codecs, random
import pandas as pd
import treetaggerwrapper
import ithmm

class posset:
	sym = {"SYM", "SENT", ":", ",", "$", "(", ")", "'", "\""}
	nn = {"NN", "NNS", "NP", "NPS"}
	jj = {"JJ", "JJR", "JJS"}
	rb = {"RB", "RBR", "RBS"}
	vb = {"VB", "VBD", "VBG", "VBN", "VBZ", "VBP"}
	vd = {"VD", "VDD", "VDG", "VDN", "VDZ", "VDP"}
	vh = {"VH", "VH", "VHG", "VHN", "VHZ", "VHP", "VHD"}
	vv = {"VV", "VV", "VVG", "VVN", "VVZ", "VVP", "VVD"}

# 品詞をまとめる
def collapse_pos(pos):
	if pos in posset.sym:
		return "SYM"
	if pos in posset.nn:
		return "NN"
	if pos in posset.rb:
		return "RB"
	if pos in posset.vb:
		return "VB"
	if pos in posset.vd:
		return "VD"
	if pos in posset.vh:
		return "VH"
	if pos in posset.vv:
		return "VV"
	if pos in posset.jj:
		return "JJ"
	return pos

def parse_tagger_result_str(result_str):
	result = result_str.split("\t")
	if len(result) == 1:		# URLなど
		word = result[0]
		if word == "<eos>":
			pos = "EOS"
		else:
			match = re.search(r"<([^ ]+)", word)
			word = "<" + match.group(1) + ">"
			pos = "SYM"
	else:
		word, pos, orig = result
		word = word.lower()
		if orig == "@card@":
			word = "##"
		if orig == "@ord@":
			word = "##"
	pos = collapse_pos(pos)
	return pos, word

class stdout:
	BOLD = "\033[1m"
	END = "\033[0m"
	CLEAR = "\033[2K"

def add_file_to_dataset(filename, dataset):
	num_words_in_train_dataset = 0
	print(stdout.BOLD + "データを準備しています ..." + stdout.END)
	word_count = set()	# 単語の種類の総数
	sentences = []
	with codecs.open(filename, "r", "utf-8") as f:
		for line in f:
			sentences.append(line)
	# データをシャッフル
	random.shuffle(sentences)
	train_split = int(len(sentences) * args.train_split)

	# 形態素解析
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	for i, data in enumerate(sentences):
		data = re.sub(ur"\n", "", data)
		data = re.sub(ur" +$", "",  data)	# 行末の空白を除去
		data = re.sub(ur"^ +", "",  data)	# 行頭の空白を除去
		if i % 10 == 0:
			sys.stdout.write("\r{}行目を処理中です ...".format(i + 1))
			sys.stdout.flush()
		results = tagger.tag_text(data)
		if len(results) == 0:
			continue
		# 形態素解析を行いながら訓練データも作る
		# 英語は通常スペース区切りなので不要と思うかもしれないが、TreeTaggerを使うと$600が$ 600に分割されたりする
		# そのためplot_en.pyで評価の際に文の単語数が[スペース区切り]と[TreeTagger]で異なる場合があり正しく評価を行えなくなる
		# よって単語分割は全てTreeTaggerによるものに統一しておく
		segmentation = ""
		for result_str in results:
			pos, lowercase = parse_tagger_result_str(result_str)
			word_count.add(lowercase)
			segmentation += lowercase + " "
		segmentation = re.sub(r" +$", "",  segmentation)	# 行末の空白を除去
		if i > train_split:
			dataset.add_test_data(segmentation)	# テストデータに追加
		else:
			dataset.add_train_data(segmentation)	# 学習用データに追加
			num_words_in_train_dataset += len(results)

	sys.stdout.write("\r")
	sys.stdout.flush()

	return dataset, num_words_in_train_dataset

def main():
	assert args.train_filename is not None
	try:
		os.mkdir(args.model)
	except:
		pass

	# 単語辞書
	dictionary = ithmm.dictionary()

	# データセット
	dataset = ithmm.dataset(dictionary)

	# 訓練データを追加
	dataset, num_words_in_train_dataset = add_file_to_dataset(args.train_filename, dataset)
	dataset.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える

	# 単語辞書を保存
	dictionary.save(os.path.join(args.model, "ithmm.dict"))

	# モデル
	model = ithmm.model()

	# ハイパーパラメータの設定
	model.set_alpha(random.uniform(10, 20))
	model.set_gamma(random.uniform(0.5, 1))
	model.set_lambda_alpha(random.uniform(0.1, 0.5))
	model.set_lambda_gamma(random.uniform(0.001, 0.05))	# 1にすればオリジナルのiTHMMと同等
	model.set_strength(random.uniform(1, 10))			# HTSSBの集中度
	model.set_tau0(1)
	model.set_tau1(100)
	# 深さを制限する場合
	# コメントアウトか-1指定で無限大
	model.set_depth_limit(args.depth_limit)

	print("alpha:", model.get_alpha(), "gamma:", model.get_gamma(), "lambda_alpha:", model.get_lambda_alpha(), "lambda_gamma:", model.get_lambda_gamma(), "strength:", model.get_strength(), "tau0:", model.get_tau0(), "tau1:", model.get_tau1())

	# 学習の準備
	trainer = ithmm.trainer(dataset, model, dictionary)

	# 初期の割り当てをチェックする場合
	trainer.show_assigned_words_for_each_tag(dictionary, 20, False);

	# グラフプロット用
	csv_likelihood = []
	csv_perplexity = []

	# 学習ループ
	for epoch in range(1, args.epoch + 1):
		start = time.time()

		# 新しい状態をギブスサンプリング
		trainer.perform_gibbs_sampling()

		# ハイパーパラメータをサンプリング
		trainer.update_hyperparameters()

		# ログ
		elapsed_time = time.time() - start
		sys.stdout.write("\rEpoch {} / {} - {:.3f} sec - {:.1f} gibbs/s".format(epoch, args.epoch, elapsed_time, num_words_in_train_dataset / elapsed_time))		
		sys.stdout.flush()

		if epoch % 100 == 0:
			print("\n")
			trainer.show_assigned_words_for_each_tag(dictionary, 20, False);
			log_likelihood = trainer.compute_log_Pdataset_test() 
			perplexity = trainer.compute_perplexity_test()
			print("alpha:", model.get_alpha(), "gamma:", model.get_gamma(), "lambda_alpha:", model.get_lambda_alpha(), "lambda_gamma:", model.get_lambda_gamma(), "strength:", model.get_strength(), "tau0:", model.get_tau0(), "tau1:", model.get_tau1())
			print("log_likelihood:", int(log_likelihood))
			print("perplexity:", int(perplexity))
			# print "MH:", trainer.get_metropolis_hastings_acceptance_rate() 

			# モデルの保存
			model.save(os.path.join(args.model, "ithmm.model"));

			# CSV出力
			csv_likelihood.append([epoch, log_likelihood])
			data = pd.DataFrame(csv_likelihood)
			data.columns = ["epoch", "log_likelihood"]
			data.to_csv(os.path.join(args.model, "likelihood.csv"))
			
			csv_perplexity.append([epoch, perplexity])
			data = pd.DataFrame(csv_perplexity)
			data.columns = ["epoch", "perplexity"]
			data.to_csv(os.path.join(args.model, "perplexity.csv"))

def _main():
	if args.filename is None:
		raise Exception()
	try:
		os.mkdir(args.model)
	except:
		pass

	trainer = ithmm.trainer()
	num_words_in_train_dataset = 0
	print(stdout.BOLD + "データを準備しています ..." + stdout.END)
	word_count = set()	# 単語の種類の総数
	dataset = []
	with codecs.open(args.filename, "r", "utf-8") as f:
		for data in f:
			dataset.append(data)
	# データをシャッフル
	random.shuffle(dataset)
	if args.train_split is None:
		args.train_split = int(len(dataset) * 0.85)
	# 形態素解析
	tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
	for i, data in enumerate(dataset):
		data = re.sub(ur"\n", "", data)
		data = re.sub(ur" +$", "",  data)	# 行末の空白を除去
		data = re.sub(ur"^ +", "",  data)	# 行頭の空白を除去
		if i % 10 == 0:
			sys.stdout.write("\r{}行目を処理中です ...".format(i))
			sys.stdout.flush()
		results = tagger.tag_text(data)
		if len(results) == 0:
			continue
		# 形態素解析を行いながら訓練データも作る
		# 英語は通常スペース区切りなので不要と思うかもしれないが、TreeTaggerを使うと$600が$ 600に分割されたりする
		# そのためplot_en.pyで評価の際に文の単語数が[スペース区切り]と[TreeTagger]で異なる場合があり正しく評価を行えなくなる
		# よって単語分割は全てTreeTaggerによるものに統一しておく
		segmentation = ""
		for result_str in results:
			pos, lowercase = parse_tagger_result_str(result_str)
			word_count.add(lowercase)
			segmentation += lowercase + " "
		segmentation = re.sub(r" +$", "",  segmentation)	# 行末の空白を除去
		if i > args.train_split:
			trainer.add_test_data(segmentation)	# テストデータに追加
		else:
			trainer.add_train_data(segmentation)	# 学習用データに追加
			num_words_in_train_dataset += len(results)

	print("\n", stdout.BOLD, "訓練データ数:", args.train_split )
	print("テストデータ数:", len(dataset) - args.train_split, stdout.END)

	# ハイパーパラメータの設定
	trainer.set_alpha(random.uniform(10, 20))
	trainer.set_gamma(random.uniform(0.5, 1))
	trainer.set_lambda_alpha(random.uniform(0.1, 0.5))
	trainer.set_lambda_gamma(random.uniform(0.001, 0.05))	# 1にすればオリジナルのiTHMMと同等
	trainer.set_strength(random.uniform(1, 10))			# HTSSBの集中度
	trainer.set_tau0(1)
	trainer.set_tau1(100)

	# 深さを制限する場合
	# コメントアウトか-1指定で無限大
	trainer.set_depth_limit(args.depth_limit)

	trainer.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える
	trainer.compile()	# 品詞をランダムに割り当てる初期化
	trainer.show_assigned_words_for_each_tag(20, False);
	print("alpha:", trainer.get_alpha(), "gamma:", trainer.get_gamma(), "lambda_alpha:", trainer.get_lambda_alpha(), "lambda_gamma:", trainer.get_lambda_gamma(), "strength:", trainer.get_strength(), "tau0:", trainer.get_tau0(), "tau1:", trainer.get_tau1())

	# グラフプロット用
	csv_likelihood = []
	csv_perplexity = []

	for epoch in range(1, args.epoch + 1):
		start = time.time()
		trainer.perform_gibbs_sampling()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec - {:.1f} gibbs/s\r".format(epoch, args.epoch, elapsed_time, num_words_in_train_dataset / elapsed_time))		
		sys.stdout.flush()
		if epoch % 100 == 0:
			print("\n")
			trainer.show_assigned_words_for_each_tag(20, False);
			log_likelihood = trainer.compute_log_Pdataset_test() 
			perplexity = trainer.compute_perplexity_test()
			print("alpha:", trainer.get_alpha(), "gamma:", trainer.get_gamma(), "lambda_alpha:", trainer.get_lambda_alpha(), "lambda_gamma:", trainer.get_lambda_gamma(), "strength:", trainer.get_strength(), "tau0:", trainer.get_tau0(), "tau1:", trainer.get_tau1())
			print("log_likelihood:", int(log_likelihood))
			print("perplexity:", int(perplexity))
			# print "MH:", trainer.get_metropolis_hastings_acceptance_rate() 
			trainer.save(args.model);
			# CSV出力
			csv_likelihood.append([epoch, log_likelihood])
			data = pd.DataFrame(csv_likelihood)
			data.columns = ["epoch", "log_likelihood"]
			data.to_csv("{}/likelihood.csv".format(args.model))
			csv_perplexity.append([epoch, perplexity])
			data = pd.DataFrame(csv_perplexity)
			data.columns = ["epoch", "perplexity"]
			data.to_csv("{}/perplexity.csv".format(args.model))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-file", "--train-filename", type=str, default=None, help="訓練用のテキストファイルのパス.ディレクトリも可.")
	parser.add_argument("-epoch", "--epoch", type=int, default=1000000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="モデル保存フォルダ名.")
	parser.add_argument("-unk", "--unknown-threshold", type=int, default=0, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("-depth", "--depth-limit", type=int, default=-1, help="最大の深さ.")
	parser.add_argument("-split", "--train-split", type=float, default=0.9, help="テキストデータの何割を訓練データにするか.")
	args = parser.parse_args()
	main()