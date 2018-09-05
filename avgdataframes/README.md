# avgdataframes.py

## python script for averaging over several data frames using pandas

```
python avgdataframes.py -h
```

NB! The script removes any columns with header "std", assuming data in these columns are standard deviations.

NB! If several columns have same header name, this is also problematic. Which is the reason for removing (pop) all the columns with name "std" by hard code.
