# -*- coding: utf-8 -*-
import argparse, codecs, sys, re
import MeCab

def main(args):
	if args.input is None:
		raise Exception()
	if args.output is None:
		raise Exception()

	dataset = []
	with codecs.open(args.input, "r", "utf-8") as f:
		tagger = MeCab.Tagger()
		for i, line in enumerate(f):
			segmentation = ""
			line = re.sub(ur"\n", "", line)
			line = re.sub(ur"[一二三四五六七八九十百千万億]+", "##", line)
			line = re.sub(ur"[0-9]+", "##", line)
			string = line.encode("utf-8")
			m = tagger.parseToNode(string)
			while m:
				segmentation += m.surface + " "
				m = m.next
			segmentation = re.sub(ur" $", "",  segmentation)
			segmentation = re.sub(ur"^ ", "",  segmentation)
			dataset.append(segmentation.decode("utf-8"))
			if i % 100 == 0:
				sys.stdout.write("\r{}行目を処理しています ...".format(i))
				sys.stdout.flush()
		sys.stdout.write("\n")

	with codecs.open(args.output, "w", "utf-8") as f:
		for line in dataset:
			f.write(line)
			f.write("\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type=str, default=None, help="分かち書き前のテキストファイルのパス.")
	parser.add_argument("-o", "--output", type=str, default=None, help="分かち書き後のテキストを保存するファイルのパス.")
	args = parser.parse_args()
	main(args)