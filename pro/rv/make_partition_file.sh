#!/bin/bash
#input: startk endk basedir

for ((kk=$1; kk <= $2; kk++))
  do
  for ww in 1.0 4.0 9.0
    do
    grep Partition ${3}fitv${kk}_bestV_V${ww}_sample6_c.log | sed -n 's/[^0-9]*\([0-9]*\)\.\([0-9]*\)/\1.\2/p' >> partition.dat
  done
done
