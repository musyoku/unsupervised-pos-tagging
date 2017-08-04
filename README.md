## Python tools for unsupervised POS tagging

教師なし品詞推定の論文4本の実装・解説を目標にしています。

実装状況

- [x] [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)

- [x] [The Infinite Hidden Markov Model](http://mlg.eng.cam.ac.uk/zoubin/papers/ihmm.pdf)

- [x] [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)

- [ ] [Embedded HMM](https://papers.nips.cc/paper/2391-inferring-state-sequences-for-non-linear-systems-with-embedded-hidden-markov-models.pdf)に基づくiTHMMのforward-backwardによる学習

## データセット

### 英語

#### Alice's Adventures in Wonderland

[http://www.gutenberg.org/files/11/11-0.txt](http://www.gutenberg.org/files/11/11-0.txt)からダウンロードできます。

`text/alice.txt`は上記データに前処理を施したものになります。

#### Penn TreeBank

[https://github.com/wojzaremba/lstm/tree/master/data](https://github.com/wojzaremba/lstm/tree/master/data)からPenn TreeBankのテキストデータをダウンロードできます。

`text/ptb.txt`は上記データの`ptb.train.txt`と`ptb.valid.txt`を結合したものになります。

### 日本語

[http://ch.nicovideo.jp/saikai/blomaga/ar904435](http://ch.nicovideo.jp/saikai/blomaga/ar904435)から青空文庫の全書籍のテキストデータをダウンロードできます。

`aozora.txt`は上記データから10万行を取り出し前処理を施したものになります。

## 形態素解析

学習結果の可視化をするために形態素解析器が必要になります。

### 英語

[TreeTagger](http://www.cis.uni-muenchen.de/~schmid/tools/TreeTagger/)を使用します。

```
pip install treetaggerwrapper
```

### 日本語

MeCabを使用します。

Pythonから利用できるようにインストールしておいてください。