#!/bin/sh

set -e
set -x

g++-6 ./test/test.cpp -lstdc++fs -o test.o -D$CPPV
