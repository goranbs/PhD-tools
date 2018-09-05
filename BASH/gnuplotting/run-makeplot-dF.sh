#!/bin/bash

### -- Run makeplot-dF.sh


for i in `seq 13 25` ; do

	sed -i -e '4c\''sep='${i}'' makeplot-dF.sh
	./makeplot-dF.sh

done
