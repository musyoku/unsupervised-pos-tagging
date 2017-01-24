# -*- coding: utf-8 -*-
import argparse, sys, os, time, re, codecs
import MeCab
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

	hmm = model.bayesian_hmm()
	# 訓練データを分かち書きする
	print stdout.BOLD + "データを準備しています ..." + stdout.END
	word_count = set()	# 単語の種類の総数
	pos_count = set()	# 品詞数
	major_pos_count = set()	# 品詞数（大分類）
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = MeCab.Tagger()
		for i, line in enumerate(f):
			if i % 100 == 0:
				sys.stdout.write("\r{}行目を処理中です ...".format(i))
				sys.stdout.flush()
			segmentation = ""
			line = re.sub(ur"\n", "", line)	# 開業を消す
			string = line.encode("utf-8")
			m = tagger.parseToNode(string)
			while m:
				word = m.surface
				word_count.add(word)
				features = m.feature.split(",")
				major_pos_count.add(features[0].decode("utf-8"))
				pos = (features[0] + "," + features[1]).decode("utf-8")
				pos_count.add(pos)
				if pos == u"名詞,数":
					word = "##"		# 数字は全て置き換える
				segmentation += word + " "
				m = m.next
			segmentation = re.sub(ur" +$", "",  segmentation)	# 行末の空白を除去
			segmentation = re.sub(ur"^ +", "",  segmentation)	# 行頭の空白を除去
			hmm.add_line(segmentation.decode("utf-8"))	# 学習用データに追加

	print stdout.END
	print stdout.BOLD + "単語数:", len(word_count), stdout.END
	print stdout.BOLD + "品詞数:", len(pos_count), stdout.END
	print repr(pos_count).decode("unicode-escape").encode("utf-8")
	print stdout.BOLD + "品詞数（大分類）:", len(major_pos_count), stdout.END
	print repr(major_pos_count).decode("unicode-escape").encode("utf-8")
	print stdout.BOLD + "1文あたりの単語数:　{}（最大）- {} 最小".format(hmm.get_max_num_words_in_line(), hmm.get_min_num_words_in_line()), stdout.END

	# Wtに制限をかけない場合
	Wt = [len(word_count)] * args.num_tags
	print "Wt:", Wt

	hmm.set_num_tags(len(Wt));	# 品詞数を設定
	hmm.initialize()	# 学習の準備

	# Wtをセット
	hmm.set_Wt(Wt)

	hmm.set_temperature(2)	# 温度の初期設定
	hmm.set_minimum_temperature(0.08)	# 温度の下限
	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		hmm.perform_gibbs_sampling()
		hmm.sample_new_beta()

		elapsed_time = time.time() - start
		sys.stdout.write("\rEpoch {} / {} - {:.3f} sec".format(epoch, args.epoch, elapsed_time))		
		sys.stdout.flush()
		hmm.anneal_temperature(0.9989)	# 温度を下げる
		if epoch % 10 == 0:
			print "\n"
			hmm.show_beta()
			hmm.show_random_line(20, True);	# ランダムなn個の文と推定結果のタグを表示
			hmm.show_typical_words_for_each_tag(20);	# それぞれのタグにつき上位n個の単語を表示
			print "temperature: ", hmm.get_temperature()
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス. 分かち書きされていない必要がある.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--num-tags", type=int, default=20, help="タグの種類.")
	main(parser.parse_args())