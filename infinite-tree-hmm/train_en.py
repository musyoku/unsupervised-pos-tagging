# coding: utf-8
import argparse, sys, os, time, re, codecs, random
import treetaggerwrapper
import model

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
		lowercase = result[0]
		match = re.search(r"<([^ ]+) ", lowercase)
		lowercase = "<" + match.group(1) + ">"
		pos = "SYM"
	else:
		word, pos, lowercase = result
	if lowercase == "@card@":
		lowercase = "##"
	if lowercase == "@ord@":
		lowercase = "##"
	pos = collapse_pos(pos)
	return pos, lowercase

class stdout:
	BOLD = "\033[1m"
	END = "\033[0m"
	CLEAR = "\033[2K"

def main(args):
	if args.filename is None:
		raise Exception()
	try:
		os.mkdir(args.model)
	except:
		pass

	ithmm = model.ithmm()
	num_words_in_train_dataset = 0
	print stdout.BOLD + "データを準備しています ..." + stdout.END
	word_count = set()	# 単語の種類の総数
	# 似たような品詞をまとめる
	# https://courses.washington.edu/hypertxt/csar-v02/penntable.html
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
			ithmm.add_test_data(segmentation)	# 学習用データに追加
		else:
			ithmm.add_train_data(segmentation)	# 学習用データに追加
			num_words_in_train_dataset += len(results)

	print "\n", stdout.BOLD, "訓練データ数:", args.train_split 
	print "テストデータ数:", len(dataset) - args.train_split, stdout.END

	# ハイパーパラメータの設定
	ithmm.set_alpha(6)
	ithmm.set_gamma(0.5)
	ithmm.set_lambda(0.02)
	ithmm.set_strength(3)
	ithmm.set_tau0(1)
	ithmm.set_tau1(100)

	# 深さを制限する場合
	# コメントアウトか-1指定で無限大
	ithmm.set_depth_limit(args.depth_limit)

	ithmm.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える
	ithmm.compile()	# 品詞をランダムに割り当てる初期化
	ithmm.set_metropolis_hastings_enabled(False)

	for epoch in xrange(1, args.epoch + 1):
		start = time.time()
		ithmm.perform_gibbs_sampling()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec - {:.1f} gibbs/s\r".format(epoch, args.epoch, elapsed_time, num_words_in_train_dataset / elapsed_time))		
		sys.stdout.flush()
		if epoch % 100 == 0:
			print "\n"
			ithmm.show_typical_words_for_each_tag(20, False);
			print "alpha:", ithmm.get_alpha(), "gamma:", ithmm.get_gamma(), "lambda:", ithmm.get_lambda(), "strength:", ithmm.get_strength(), "tau0:", ithmm.get_tau0(), "tau1:", ithmm.get_tau1()
			print "logP(x):", ithmm.compute_log_Pdataset_test() 
			print "PPL:", ithmm.compute_perplexity_test() 
			print "MH:", ithmm.get_metropolis_hastings_acceptance_rate() 
			ithmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=1000000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-u", "--unknown-threshold", type=int, default=0, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("-d", "--depth-limit", type=int, default=-1, help="最大の深さ.")
	parser.add_argument("-l", "--train-split", type=int, default=None, help="テキストデータの最初の何行を訓練データにするか.")
	main(parser.parse_args())