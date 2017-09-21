import argparse, codecs, re
import numpy as np

def main(args):
	num_state = 4
	num_symbol = 6

	# 状態遷移確率
	A = np.zeros((num_state, num_state), dtype=np.float64)
	A[0] = 0.4, 0.1, 0.4, 0.1,
	A[1] = 0.2, 0.3, 0.2, 0.3,
	A[2] = 0.3, 0.2, 0.1, 0.4,
	A[3] = 0.1, 0.4, 0.3, 0.2,

	# 出力確率
	B = np.zeros((num_state, num_symbol), dtype=np.float64)
	B[0] = 0.3, 0.7, 0.0, 0.0, 0.0, 0.0
	B[1] = 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
	B[2] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
	B[3] = 0.0, 0.0, 0.0, 0.0, 0.7, 0.3

	with codecs.open("text/test.txt", "w", "utf-8") as f:
		for n in range(args.num_seq):
			state = 0
			sequence = ""
			for l in range(args.seq_length):
				state = int(np.argwhere(np.random.multinomial(1, A[state]) == 1))
				emission = int(np.argwhere(np.random.multinomial(1, B[state]) == 1))
				sequence += str(emission) + " "
			sequence = re.sub(r" $", "", sequence)
			print(sequence)
			f.write(sequence + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-l", "--seq-length", type=int, default=20, help="1つの文の長さ.")
	parser.add_argument("-n", "--num-seq", type=int, default=20, help="生成する文の個数.")
	main(parser.parse_args())