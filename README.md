## Python Tools for Unsupervised POS Tagging

教師なし品詞推定の論文4本の実装・解説を目標にしています。

実装状況

- [x] [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)

- [x] [The Infinite Hidden Markov Model](http://mlg.eng.cam.ac.uk/zoubin/papers/ihmm.pdf)

- [x] [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)

- [ ] [Embedded HMM](https://papers.nips.cc/paper/2391-inferring-state-sequences-for-non-linear-systems-with-embedded-hidden-markov-models.pdf)に基づくiTHMMのforward-backwardによる学習

## データセット

### 英語

#### Penn TreeBank

[https://github.com/wojzaremba/lstm/tree/master/data](https://github.com/wojzaremba/lstm/tree/master/data)からPenn TreeBankのテキストデータをダウンロードできます。

`text/ptb.txt`は上記データの`ptb.train.txt`と`ptb.valid.txt`を結合したものになります。

### 日本語

#### こころ

[http://www.aozora.gr.jp/cards/000148/card773.html](http://www.aozora.gr.jp/cards/000148/card773.html)からダウンロードできます。

`text/kokoro.txt`は上記データに前処理を施したものになります。

#### 吾輩は猫である

[http://www.aozora.gr.jp/cards/000148/card789.html](http://www.aozora.gr.jp/cards/000148/card789.html)からダウンロードできます。

`text/neko.txt`は上記データに前処理を施したものになります。