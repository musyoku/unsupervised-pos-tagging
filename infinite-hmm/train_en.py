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

	hmm = model.ihmm(args.initial_num_tags)

	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
	print stdout.BOLD + "データを準備しています ..." + stdout.END
	Wt_count = {}
	word_count = set()	# 単語の種類の総数
	# 似たような品詞をまとめる
	# https://courses.washington.edu/hypertxt/csar-v02/penntable.html
	with codecs.open(args.filename, "r", "utf-8") as f:
		tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
		for i, line in enumerate(f):
			if args.train_split is not None and i > args.train_split:
				break
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
				poses = poses.split("\t")
				if len(poses) == 1:
					lowercase = poses[0]
				else:
					word, pos, lowercase = poses
				if lowercase == "@card@":
					lowercase = "##"
				if lowercase == "@ord@":
					lowercase = "##"
				word_count.add(lowercase)
				segmentation += lowercase + " "
				pos = collapse_pos(pos)
				if pos not in Wt_count:
					Wt_count[pos] = {}
				if lowercase not in Wt_count[pos]:
					Wt_count[pos][lowercase] = 1
				else:
					Wt_count[pos][lowercase] += 1
			segmentation = re.sub(r" +$", "",  segmentation)	# 行末の空白を除去
			hmm.add_line(segmentation)	# 学習用データに追加

	hmm.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える
	hmm.initialize()	# 品詞数をセットしてから初期化

	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		if args.beam:
			hmm.perform_beam_sampling()
		else:
			hmm.perform_gibbs_sampling()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec\r".format(epoch, args.epoch, elapsed_time))		
		sys.stdout.flush()
		if epoch % 10 == 0:
			print "\n"
			hmm.show_typical_words_for_each_tag(20);
			hmm.show_log_Pdata();
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--initial-num-tags", type=int, default=20, help="品詞の個数.")
	parser.add_argument("-u", "--unknown-threshold", type=int, default=1, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("-l", "--train-split", type=int, default=None, help="テキストデータの最初の何行を訓練データにするか.")
	parser.add_argument("--beam", default=False, action="store_true", help="品詞の個数.")
	main(parser.parse_args())