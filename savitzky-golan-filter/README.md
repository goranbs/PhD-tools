# smoothing.py
Apply Savitzky-Golan filter to noisy data

I have used this technique to filter output from umbrella samplings.

##  Useage:

`smoothing.py <inputfile.dat> colx coly [wsize] [polyn]`

* The `inputfile` is assumed to contain columnar data.
* Comment lines in the input file is assumed to start with '#'
* `colx` is the column number containing x data
* `coly` is the column number containing f(x) data
* Column counting starts at 1 for the first column
* `wsize` is the window size of the SG filter (optional)
* `polyn` is the polynomial degree in the SG filter (optional)

