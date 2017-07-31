CC = g++
INCLUDE = -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -I /usr/local/include
CFLAGS = -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O0 -g
CFLAGS_SO = -shared -fPIC -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O2

install: ## Python用ライブラリをビルドします.
	$(CC) python.cpp ithmm/tssb.cpp ithmm/node.cpp ithmm/table.cpp ithmm/hpylm.cpp -o ithmm.so $(INCLUDE) $(CFLAGS_SO)

test: ## LLDB用.
	$(CC) test.cpp ithmm/tssb.cpp ithmm/node.cpp ithmm/table.cpp ithmm/hpylm.cpp $(INCLUDE) $(CFLAGS)

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help