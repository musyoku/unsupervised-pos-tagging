## A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging

- [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)
- [実装について](http://musyoku.github.io/2017/01/28/A-Fully-Bayesian-Approach-to-Unsupervised-Part-of-Speech-Tagging/)

## ビルド

```
make install
```

## 学習

### 英語

```
python train_en.py -f ../alice.txt -n 7
```

### 日本語

```
python train_ja.py -f ../aozora.txt -n 20
```

## 結果の可視化

```
python tags.py -n 200
```

獲得した品詞とそれに属する単語を一覧表示します。