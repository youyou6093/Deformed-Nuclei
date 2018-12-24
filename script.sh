#!/bin/bash

for i in {0..199}
do
	./main "$i" "20" "28" "Ca48" "5"
done


for i in {0..199}
do
	./main "$i" "28" "24" "Ni52" "5"
done

for i in {0..199}
do
	./main "$i" "24" "28" "Cr52" "5"
done


for i in {0..199}
do
	./main "$i" "28" "26" "Ni54" "5"
done

for i in {0..199}
do
	./main "$i" "26" "28" "Fe54" "5"
done


