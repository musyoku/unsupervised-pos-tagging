# -*- coding: utf-8 -*-
import argparse, sys, os, time, re, codecs
import treetaggerwrapper
import model

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

	print stdout.BOLD + "単語数を計算しています ..." + stdout.END
	word_count = set()	# 単語の種類の総数
	with codecs.open(args.filename, "r", "utf-8") as f:
		for i, line in enumerate(f):
			line = re.sub(ur"\n", "", line)
			sys.stdout.write("\r{}行目を処理中です ...".format(i))
			sys.stdout.flush()
			words = line.split(" ")
			for word in words:
				word_count.add(word)
		print stdout.END
		print stdout.BOLD + "単語数:", len(word_count), stdout.END

	# Wtに制限をかけない場合
	Wt = [len(word_count)] * args.num_tags
	print "Wt:", Wt

	hmm = model.bayesian_hmm()
	hmm.set_num_tags(len(Wt));	# 品詞数を設定

	# テキストファイルの読み込み
	# 複数のファイルを読んでもOK
	hmm.load_textfile(args.filename)

	# 全てのテキストファイルを読み込み終わってから初期化
	hmm.initialize()

	# Wtをセット
	hmm.set_Wt(Wt)

	hmm.set_temperature(2)	# 温度の初期設定
	hmm.set_minimum_temperature(0.08)	# 温度の下限
	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		hmm.perform_gibbs_sampling()
		hmm.sample_new_beta()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec\r".format(epoch, args.epoch, elapsed_time))		
		sys.stdout.flush()
		hmm.anneal_temperature(0.998)	# 温度を下げる
		if epoch % 1 == 0:
			print "\n"
			hmm.show_beta()
			hmm.show_random_line(20, True);	# ランダムなn個の文と推定結果のタグを表示
			hmm.show_typical_words_for_each_tag(20);	# それぞれのタグにつき上位n個の単語を表示
			print "temperature: ", hmm.get_temperature()
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--num-tags", type=int, default=20, help="タグの種類.")
	main(parser.parse_args())