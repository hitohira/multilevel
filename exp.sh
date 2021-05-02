#!/bin/sh

make
./powerlaw_exp 2> /dev/null > $1
