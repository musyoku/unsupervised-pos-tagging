# -*- coding: utf-8 -*-
import argparse, codecs, re
import numpy as np

def main(args):
	num_state = 3
	num_symbol = 6

	# 状態遷移確率
	A = np.zeros((num_state, num_state), dtype=np.float64)
	A[0] = 0.1, 0.7, 0.2
	A[1] = 0.2, 0.1, 0.7
	A[2] = 0.7, 0.2, 0.1

	# 出力確率
	B = np.zeros((num_state, num_symbol), dtype=np.float64)
	B[0] = 0.9, 0.1, 0, 0, 0, 0
	B[1] = 0, 0, 0.6, 0.4, 0, 0
	B[2] = 0, 0, 0, 0, 0.1, 0.9

	with codecs.open("../test.txt", "w", "utf-8") as f:
		for n in xrange(args.num_seq):
			state = 0
			sequence = ""
			for l in xrange(args.seq_length):
				state = int(np.argwhere(np.random.multinomial(1, A[state]) == 1))
				emission = int(np.argwhere(np.random.multinomial(1, B[state]) == 1))
				sequence += str(emission) + " "
			sequence = re.sub(r" $", "", sequence)
			print sequence
			f.write(sequence + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-l", "--seq-length", type=int, default=20, help="品詞の個数.")
	parser.add_argument("-n", "--num-seq", type=int, default=20, help="品詞の個数.")
	main(parser.parse_args())