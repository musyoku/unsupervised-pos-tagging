CC = g++
INCLUDE = -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -I /usr/local/include
CFLAGS = -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O3
SOFLAGS = -shared -fPIC -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O3

install: ## Python用ライブラリをビルドします.
	$(CC) src/python.cpp src/bhmm/*.cpp src/python/*.cpp -o run/bhmm.so $(INCLUDE) $(SOFLAGS)

cpp: ## LLDB用.
	$(CC) test/train.cpp src/bhmm/*.cpp src/python/*.cpp -o test/train $(INCLUDE) $(CFLAGS)

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help