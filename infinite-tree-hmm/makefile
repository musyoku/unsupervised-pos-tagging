CC = g++
INCLUDE = -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -I /usr/local/include
CFLAGS = -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O0 -g
CFLAGS_SO = -shared -fPIC -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O2

install: ## Python用ライブラリをビルドします.
	$(CC) ithmm.cpp src/tssb.cpp src/node.cpp src/table.cpp src/hpylm.cpp -o ithmm.so $(INCLUDE) $(CFLAGS_SO)

test: ## LLDB用.
	$(CC) test.cpp src/tssb.cpp src/node.cpp src/table.cpp src/hpylm.cpp $(INCLUDE) $(CFLAGS)

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help