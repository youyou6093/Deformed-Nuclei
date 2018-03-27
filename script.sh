for i in {2..5};
do
a=0.02
echo "$a*$i"|bc
c="$a*$i"|bc
echo $[]c
./main 20 20
#mv density${-0.01*i}.dat ./isodipole
done
