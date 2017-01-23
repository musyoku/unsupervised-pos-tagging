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
	vh = {"VH", "VH", "VHG", "VHN", "VHZ", "VHP"}
	vv = {"VV", "VV", "VVG", "VVN", "VVZ", "VVP"}

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
	return pos

def main(args):
	if args.filename is None:
		raise Exception()
	try:
		os.mkdir(args.model)
	except:
		pass

	hmm = model.hpylm_hmm(args.num_tags)

	# テキストファイルの読み込み
	# 複数のファイルを読んでもOK
	hmm.load_textfile(args.filename)

	# 学習
	hmm.prepare_for_training()
	for epoch in xrange(1, args.epoch + 1):
		start = time.time()

		hmm.perform_gibbs_sampling()
		hmm.sample_hyperparams()

		elapsed_time = time.time() - start
		ppl = hmm.compute_perplexity()
		sys.stdout.write("\r\33[2K");	# 行をクリア
		sys.stdout.write("Epoch {} / {} - {:.3f} sec - {:.3f} ppl\n".format(epoch, args.epoch, elapsed_time, ppl))		
		sys.stdout.flush()
		if epoch % 10 == 0:
			hmm.dump_hpylm()
			hmm.show_typical_words_for_each_tag(20)
			hmm.save(args.model);

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename", type=str, default=None, help="訓練用のテキストファイルのパス.")
	parser.add_argument("-e", "--epoch", type=int, default=20000, help="総epoch.")
	parser.add_argument("-m", "--model", type=str, default="out", help="保存フォルダ名.")
	parser.add_argument("-n", "--num_tags", type=int, default=30, help="品詞数.")
	main(parser.parse_args())