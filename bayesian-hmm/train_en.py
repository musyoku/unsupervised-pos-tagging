# -*- coding: utf-8 -*-
import argparse, sys, os, time, re, codecs
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
def colapse_pos(pos):
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
	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
	print stdout.BOLD + "データを準備しています ..." + stdout.END
	Wt_count = {}
	word_count = set()	# 単語の種類の総数
	# 似たような品詞をまとめる
	# https://courses.washington.edu/hypertxt/csar-v02/penntable.html
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
		for i, line in enumerate(f):
			line = re.sub(ur"\n", "", line)
			line = re.sub(ur" +$", "",  line)	# 行末の空白を除去
			line = re.sub(ur"^ +", "",  line)	# 行頭の空白を除去
			sys.stdout.write("\r{}行目を処理中です ...".format(i))
			sys.stdout.flush()
			result = tagger.tag_text(line)
			if len(result) == 0:
				continue
			# 形態素解析を行いながら訓練データも作る
			# 英語は通常スペース区切りなので不要と思うかもしれないが、TreeTaggerを使うと$600が$ 600に分割されたりする
			# そのためplot_en.pyで評価の際に文の単語数が[スペース区切り]と[TreeTagger]で異なる場合があり正しく評価を行えなくなる
			# よって単語分割は全てTreeTaggerによるものに統一しておく
			segmentation = ""
			for poses in result:
				word, pos, lowercase = poses.split("\t")
				word_count.add(lowercase)
				segmentation += lowercase + " "
				pos = colapse_pos(pos)
				if pos not in Wt_count:
					Wt_count[pos] = {}
				if lowercase not in Wt_count[pos]:
					Wt_count[pos][lowercase] = 1
				else:
					Wt_count[pos][lowercase] += 1
			segmentation = re.sub(ur" +$", "",  segmentation)	# 行末の空白を除去
			hmm.add_line(segmentation.decode("utf-8"))	# 学習用データに追加
	if args.supervised:
		# Wtは各タグについて、そのタグになりうる単語の数が入っている
		# タグ0には<bos>と<eos>だけ含まれることにする
		Wt = [2]
		for tag, words in Wt_count.items():
			print tag, ":", len(words)
			if len(words) < 10:
				print words
			Wt.append(len(words))
	else:
		# Wtに制限をかけない場合
		Wt = [len(word_count)] * args.num_tags

	print "Wt:", Wt

	hmm.set_num_tags(len(Wt));	# 品詞数を設定
	hmm.initialize()	# 品詞数をセットしてから初期化

	# Wtをセット
	hmm.set_Wt(Wt)

	hmm.set_temperature(2)	# 温度の初期設定
	hmm.set_minimum_temperature(0.08)	# 温度の下限
	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		hmm.perform_gibbs_sampling()
		hmm.sample_new_alpha()
		hmm.sample_new_beta()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec\r".format(epoch, args.epoch, elapsed_time))		
		sys.stdout.flush()
		hmm.anneal_temperature(0.9989)	# 温度を下げる
		if epoch % 10 == 0:
			print "\n"
			hmm.show_alpha()
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
	parser.add_argument("-s", "--supervised", dest="supervised", default=True, action="store_true", help="各タグのWtを訓練データで制限するかどうか.")
	parser.add_argument("-u", "--unsupervised", dest="supervised", action="store_false", help="各タグのWtを訓練データで制限するかどうか.")
	parser.add_argument("-n", "--num-tags", type=int, default=20, help="タグの種類（semi_supervisedがFalseの時のみ有効）.")
	parser.add_argument("--start-temperature", type=float, default=2, help="開始温度.")
	parser.add_argument("--min-temperature", type=float, default=0.08, help="最小温度.")
	parser.add_argument("--anneal", type=float, default=0.9989, help="温度の減少に使う係数.")
	main(parser.parse_args())