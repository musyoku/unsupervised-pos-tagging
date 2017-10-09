CC = g++
BOOST = /usr/local/Cellar/boost/1.65.0
INCLUDE = `python3-config --includes` -std=c++11 -I$(BOOST)/include
LDFLAGS = `python3-config --ldflags` -lboost_serialization -lboost_python3 -L$(BOOST)/lib
SOFLAGS = -shared -fPIC

install: ## Python用ライブラリをビルドします.
	$(CC) $(INCLUDE) $(LDFLAGS) $(SOFLAGS) src/python.cpp src/bhmm/*.cpp src/python/*.cpp -o run/bhmm.so -O3

install_ubuntu: ## Python用ライブラリをビルドします.
	$(CC) -Wl,--no-as-needed -Wno-deprecated $(INCLUDE) $(LDFLAGS) $(SOFLAGS) src/python.cpp src/bhmm/*.cpp src/python/*.cpp -o run/bhmm.so -O3

check_includes:	## Python.hの場所を確認
	python3-config --includes

check_ldflags:	## libpython3の場所を確認
	python3-config --ldflags

.PHONY: test
test: ## LLDB用.
	$(CC) test/train.cpp src/bhmm/*.cpp src/python/*.cpp -o test/train $(INCLUDE) $(LDFLAGS) -O0 -g
	$(CC) test/viterbi.cpp src/bhmm/*.cpp src/python/*.cpp -o test/viterbi $(INCLUDE) -O0 -g

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help