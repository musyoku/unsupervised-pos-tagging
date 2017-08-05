CC = g++
INCLUDE = -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -I /usr/local/include
CFLAGS = -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O0 -g
SOFLAGS = -shared -fPIC -std=c++11 -L/usr/local/lib -lboost_serialization -lboost_python -lpython2.7 -O3

install: ## Python用ライブラリをビルドします.
	$(CC) src/model.cpp src/bhmm/hmm.cpp src/bhmm/sampler.cpp src/bhmm/utils.cpp -o model.so $(INCLUDE) $(SOFLAGS)

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help