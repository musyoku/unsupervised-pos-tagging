CC = g++
BOOST = /usr/local/Cellar/boost/1.65.0/include
INCLUDE = `python3-config --includes` -std=c++11 -I$(BOOST)
LDFLAGS = `python3-config --ldflags` -lboost_serialization -lboost_python3
SOFLAGS = -shared -fPIC

install: ## Python用ライブラリをコンパイル
	$(CC) $(INCLUDE) $(LDFLAGS) $(SOFLAGS) src/python.cpp src/ihmm/*.cpp src/python/*.cpp -o run/ihmm.so -O3

check_includes:	## Python.hの場所を確認
	python3-config --includes

check_ldflags:	## libpython3の場所を確認
	python3-config --ldflags

.PHONY: test
test: ## LLDB用.
	$(CC) test/test.cpp src/ihmm/*.cpp src/python/*.cpp -o test/test $(INCLUDE) $(LDFLAGS) -O3

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help