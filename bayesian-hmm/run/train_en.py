# coding: utf-8
from __future__ import print_function
import argparse, sys, os, time, codecs, random
import treetaggerwrapper
import bhmm

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
# https://courses.washington.edu/hypertxt/csar-v02/penntable.html
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

def build_corpus(filename, dataset):
	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
	print("データを準備しています ...")
	sentence_list = []
	with codecs.open(filename, "r", "utf-8") as f:
		for sentence in f:
			sentence_list.append(sentence)
	random.shuffle(sentence_list)	# データをシャッフル
	train_split = int(len(sentence_list) * args.train_split)

	word_count = set()	# 単語の種類の総数
	Wt_count = {}
	with codecs.open(filename, "r", "utf-8") as f:
		tagger = treetaggerwrapper.TreeTagger(TAGLANG="en")
		for i, sentence in enumerate(f):
			sentence = sentence.strip()
			if i % 10 == 0:
				sys.stdout.write("\r{}行目を処理中です ...".format(i))
				sys.stdout.flush()
			result = tagger.tag_text(sentence)
			if len(result) == 0:
				continue
			# 形態素解析を行いながら訓練データも作る
			# 英語は通常スペース区切りなので不要と思うかもしれないが、TreeTaggerを使うと$600が$ 600に分割されたりする
			# そのためplot_en.pyで評価の際に文の単語数が[スペース区切り]と[TreeTagger]で異なる場合があり正しく評価を行えなくなる
			# よって単語分割は全てTreeTaggerによるものに統一しておく
			words = []
			for poses in result:
				word, pos, lowercase = poses.split("\t")
				word_count.add(lowercase)
				words.append(lowercase)
				pos = collapse_pos(pos)
				if pos not in Wt_count:
					Wt_count[pos] = {}
				if lowercase not in Wt_count[pos]:
					Wt_count[pos][lowercase] = 1
				else:
					Wt_count[pos][lowercase] += 1
			# データを追加
			if i > train_split:
				dataset.add_words_dev(words)		# 評価用データに追加
			else:
				dataset.add_words_train(words)		# 学習用データに追加

	if args.supervised:
		# Wtは各タグについて、そのタグになりうる単語の数が入っている
		# タグ0には<bos>と<eos>だけ含まれることにする
		Wt = [2]
		for tag, words in Wt_count.items():
			Wt.append(len(words))
	else:
		# Wtに制限をかけない場合
		Wt = [len(word_count)] * args.num_tags

	return dataset, Wt

def main():
	assert args.train_filename is not None
	try:
		os.mkdir(args.model)
	except:
		pass

	# 単語辞書
	dictionary = bhmm.dictionary()

	# データセット
	dataset = bhmm.dataset(dictionary)

	# 訓練データを追加
	dataset, Wt = build_corpus(args.train_filename, dataset)
	dataset.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える

	# 単語辞書を保存
	dictionary.save(os.path.join(args.model, "bhmm.dict"))

	# モデル
	model = bhmm.model()
	# model.load(os.path.join(args.model, "bhmm.model"))

	# ハイパーパラメータの設定
	model.set_temperature(args.start_temperature)		# 温度の初期設定
	model.set_minimum_temperature(args.min_temperature)	# 温度の下限

	# 学習の準備
	trainer = bhmm.trainer(dataset, model, dictionary, Wt)

	# 学習ループ
	for epoch in range(1, args.epoch + 1):
		start = time.time()

		# 新しい状態をギブスサンプリング
		trainer.perform_gibbs_sampling()

		# ハイパーパラメータをサンプリング
		trainer.update_hyperparameters()

		# ログ
		elapsed_time = time.time() - start
		sys.stdout.write("\rEpoch {} / {} - {:.3f} sec - {:.1f} gibbs/s".format(epoch, args.epoch, elapsed_time, 0))		
		sys.stdout.flush()

def _main(args):
	if args.filename is None:
		raise Exception()
	try:
		os.mkdir(args.model)
	except:
		pass

	hmm = model.bayesian_hmm()
	# 訓練データを形態素解析して各品詞ごとにその品詞になりうる単語の総数を求めておく
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
			words = []
			for poses in result:
				word, pos, lowercase = poses.split("\t")
				word_count.add(lowercase)
				words.append(lowercase)
				pos = collapse_pos(pos)
				if pos not in Wt_count:
					Wt_count[pos] = {}
				if lowercase not in Wt_count[pos]:
					Wt_count[pos][lowercase] = 1
				else:
					Wt_count[pos][lowercase] += 1
			hmm.add_line(segmentation)	# 学習用データに追加
	if args.supervised:
		# Wtは各タグについて、そのタグになりうる単語の数が入っている
		# タグ0には<bos>と<eos>だけ含まれることにする
		Wt = [2]
		for tag, words in Wt_count.items():
			Wt.append(len(words))
	else:
		# Wtに制限をかけない場合
		Wt = [len(word_count)] * args.num_tags

	hmm.set_num_tags(len(Wt));	# 品詞数を設定
	hmm.mark_low_frequency_words_as_unknown(args.unknown_threshold)	# 低頻度語を全て<unk>に置き換える
	hmm.initialize()	# 品詞数をセットしてから初期化

	# Wtをセット
	hmm.set_Wt(Wt)

	hmm.set_temperature(args.start_temperature)	# 温度の初期設定
	hmm.set_minimum_temperature(args.min_temperature)	# 温度の下限
	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		hmm.perform_gibbs_sampling()
		hmm.sample_new_alpha()
		hmm.sample_new_beta()

		elapsed_time = time.time() - start
		sys.stdout.write(" Epoch {} / {} - {:.3f} sec\r".format(epoch, args.epoch, elapsed_time))		
		sys.stdout.flush()
		hmm.anneal_temperature(args.anneal)	# 温度を下げる
		if epoch % 10 == 0:
			hmm.show_alpha()
			hmm.show_beta()
			# hmm.show_random_line(20, True);	# ランダムなn個の文と推定結果のタグを表示
			hmm.show_typical_words_for_each_tag(20);	# それぞれのタグにつき上位n個の単語を表示
			print("temperature: ", hmm.get_temperature())
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-file", "--train-filename", type=str, default=None, help="訓練用のテキストファイルのパス.ディレクトリも可.")
	parser.add_argument("-epoch", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("--supervised", dest="supervised", default=True, action="store_true", help="各タグのWtを訓練データで制限するかどうか.")
	parser.add_argument("--unsupervised", dest="supervised", action="store_false", help="各タグのWtを訓練データで制限するかどうか.")
	parser.add_argument("-tags", "--num-tags", type=int, default=20, help="タグの種類（semi_supervisedがFalseの時のみ有効）.")
	parser.add_argument("-unk", "--unknown-threshold", type=int, default=0, help="出現回数がこの値以下の単語は<unk>に置き換える.")
	parser.add_argument("-split", "--train-split", type=float, default=0.9, help="テキストデータの何割を訓練データにするか.")
	parser.add_argument("--start-temperature", type=float, default=1.5, help="開始温度.")
	parser.add_argument("--min-temperature", type=float, default=0.08, help="最小温度.")
	parser.add_argument("--anneal", type=float, default=0.99989, help="温度の減少に使う係数.")
	args = parser.parse_args()
	main()