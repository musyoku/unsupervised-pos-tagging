## 無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習

- [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)
- [実装について](http://musyoku.github.io/2017/03/09/%E7%84%A1%E9%99%90%E6%9C%A8%E6%A7%8B%E9%80%A0%E9%9A%A0%E3%82%8CMarkov%E3%83%A2%E3%83%87%E3%83%AB%E3%81%AB%E3%82%88%E3%82%8B%E9%9A%8E%E5%B1%A4%E7%9A%84%E5%93%81%E8%A9%9E%E3%81%AE%E6%95%99%E5%B8%AB%E3%81%AA%E3%81%97%E5%AD%A6%E7%BF%92/)

#### Todo:

- [ ] メトロポリス・ヘイスティングス法による補正
- [ ] ハイパーパラメータのサンプリング

## 準備

### macOS

macOSの場合、PythonとBoostはともにbrewでインストールする必要があります。

#### Python 3のインストール

```
brew install python3
```

`PYTHONPATH`を変更する必要があるかもしれません。

#### Boostのインストール

```
brew install boost-python --with-python3
```

### Ubuntu

#### Boostのインストール

```
./bootstrap.sh --with-python=python3 --with-python-version=3.5
./b2 python=3.5 -d2 -j4 --prefix BOOST_DIR install
```

Pythonのバージョンを自身のものと置き換えてください。

### ビルド

以下のコマンドで`ithmm.so`が生成され、Pythonから利用できるようになります。

```
make install
```

`makefile`内のBoostのパスを自身の環境に合わせて書き換えてください。

Ubuntuでエラーが出る場合は代わりに以下を実行します。

```
make install_ubuntu
```

### MeCabのインストール

```
pip install mecab-python3
```

## 学習

英語のテキストファイルの場合は以下のコマンドで学習できます。

```
python3 train_en.py  -f ../../text/ptb.txt -split 1 -tags 1 -e 40000
```

日本語のテキストファイルの場合は以下のコマンドで学習できます。

```
python3 train_ja.py  -f ../../text/neko.txt -split 1 -tags 1 -e 40000
```

## 結果の可視化

各予測タグとそれに割り当てられた単語を表示するには以下のコマンドを実行します。

```
python3 tags.py -n 200
```

混同行列をプロットするには以下のコマンドを実行します。

```
python3 plot_ja.py -f ../../text/kokoro.txt
```

学習時に使用したテキストファイルを指定します。

