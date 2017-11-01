#!/bin/bash

for i in {1..5};
do
./main ${20} ${20} ${-0.01 * i}
mv density${-0.01*i}.dat ./isodipole 
done