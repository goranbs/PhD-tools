#!/bin/bash
rm \#*
editconf -f $1 -o conf1.gro
editconf -f $2 -o conf2.gro

head -n 1 conf1.gro  > head1.tmp
expr $(cat conf1.gro | wc -l ) + $(cat conf2.gro | wc -l ) - 6 > head2.tmp
tail -n 1 conf2.gro > tail.tmp 

sed -e '$d' conf1.gro > conf11.tmp
sed -e '1d' conf11.tmp > conf12.tmp
sed -e '1d' conf12.tmp > conf1.tmp

sed -e '$d' conf2.gro > conf22.tmp
sed -e '1d' conf22.tmp > conf23.tmp
sed -e '1d' conf23.tmp > conf2.tmp

cat head1.tmp head2.tmp conf1.tmp conf2.tmp tail.tmp > temp.gro
genconf -renumber -f temp.gro -o $3
rm *.tmp \#* conf1.gro conf2.gro temp.gro

