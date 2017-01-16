# -*- coding: utf-8 -*-
import argparse, sys
import model

def main(args):
	if args.filename is None:
		raise Exception()
	hmm = model.bayesian_hmm()
	hmm.set_num_tags(args.num_tags);	# 品詞数を設定

	# テキストファイルの読み込み
	# 複数のファイルを読んでもOK
	hmm.load_textfile(args.filename)

	# 全てのテキストファイルを読み込み終わってから初期化
	hmm.initialize()

	# 温度の調整が面倒なので1で固定
	# hmm.set_temperature(2)	# 温度の初期設定
	# hmm.set_minimum_temperature(0.08)	# 温度の下限
	for epoch in xrange(1, args.epoch + 1):
		sys.stdout.write(" Epoch {} / {}\r".format(epoch, args.epoch))		
		sys.stdout.flush()
		hmm.perform_gibbs_sampling()
		# hmm.anneal_temperature(0.998)	# 温度を下げる
		if epoch % 10 == 0:
			hmm.show_random_line(20, True);	# ランダムなn個の文と推定結果のタグを表示
			hmm.show_typical_words_for_each_tag(20);	# それぞれのタグにつき上位n個の単語を表示
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-n", "--num_tags", type=int, default=30, help="品詞の数.")
	parser.add_argument("-m", "--model", type=str, default="hmm.model", help="モデルファイル.")
	args = parser.parse_args()
	main(args)