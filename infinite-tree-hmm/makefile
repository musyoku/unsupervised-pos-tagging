CC = g++
BOOST = /usr/local/Cellar/boost/1.65.0
INCLUDE = `python3-config --includes` -std=c++11 -I$(BOOST)/include
LDFLAGS = `python3-config --ldflags` -lboost_serialization -lboost_python3 -L$(BOOST)/lib
SOFLAGS = -shared -fPIC -march=native
SOURCES = src/ithmm/*.cpp src/python/*.cpp

install: ## Python用ライブラリをコンパイル
	$(CC) $(INCLUDE) $(SOFLAGS) -o run/ithmm.so src/python.cpp $(SOURCES) $(LDFLAGS) -O3

install_ubuntu: ## Python用ライブラリをコンパイル
	$(CC) -Wl,--no-as-needed -Wno-deprecated $(INCLUDE) $(SOFLAGS) -o run/ithmm.so src/python.cpp $(SOURCES) $(LDFLAGS) -O3

check_includes:	## Python.hの場所を確認
	python3-config --includes

check_ldflags:	## libpython3の場所を確認
	python3-config --ldflags

module_tests: ## 各モジュールのテスト.
	$(CC) test/module_tests/ithmm.cpp src/ithmm/*.cpp -o test/module_tests/ithmm $(INCLUDE) $(LDFLAGS) -O0 -g
	$(CC) test/module_tests/node.cpp src/ithmm/*.cpp -o test/module_tests/node $(INCLUDE) $(LDFLAGS) -O0 -g
	$(CC) test/module_tests/table.cpp src/ithmm/*.cpp -o test/module_tests/table $(INCLUDE) $(LDFLAGS) -O0 -g
	$(CC) test/module_tests/hpylm.cpp src/ithmm/*.cpp -o test/module_tests/hpylm $(INCLUDE) $(LDFLAGS) -O0 -g
	./test/module_tests/hpylm
	./test/module_tests/table
	./test/module_tests/node
	./test/module_tests/ithmm

running_tests:	## テスト
	$(CC) test/running_tests/forward_backward.cpp $(SOURCES) -o test/running_tests/forward_backward $(INCLUDE) $(LDFLAGS) -O0 -g -Wall
	$(CC) test/running_tests/train.cpp $(SOURCES) -o test/running_tests/train $(INCLUDE) $(LDFLAGS) -O3 -Wall
	$(CC) test/running_tests/save.cpp $(SOURCES) -o test/running_tests/save $(INCLUDE) $(LDFLAGS) -O0 -g -Wall

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help