#!/bin/bash
###########################################
#
# Set the surface charge of slabs according
# to resname ID's in the data files
###########################################

FILES=./gold-q22y/gold-q22-vac*.data

q=0.0625

for file in $FILES ; do
	awk -v CONVFMT=%.16g '{if ($10=="SRF") gsub($4,0.0625); print}' $file > ${file}.out
done


