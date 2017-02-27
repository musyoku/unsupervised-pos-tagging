## Python tools for unsupervised POS tagging

教師なし品詞推定の論文3本の実装・解説を目標にしています。

実装状況

:smile: [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)

:confused: [The Infinite Hidden Markov Model](http://mlg.eng.cam.ac.uk/zoubin/papers/ihmm.pdf)

:tired_face: [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)

## データセット

### 英語

#### Alice's Adventures in Wonderland

[http://www.gutenberg.org/files/11/11-0.txt](http://www.gutenberg.org/files/11/11-0.txt)からダウンロードできます。

`alice.txt`は上記データに前処理を施したものになります。

#### Wikicorpus

[http://www.cs.upc.edu/~nlp/wikicorpus/](http://www.cs.upc.edu/~nlp/wikicorpus/)から英語版Wikipediaのテキストデータをダウンロードできます。

`wiki.txt`は上記データから10万行を取り出し前処理を施したものになります。

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