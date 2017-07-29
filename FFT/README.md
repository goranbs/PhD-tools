=== FFT.py


Python script for performing Fast Fourier Transform of input signal.

Useage:
'''
>> FFT.py -i inputfile.dat -c [1,2] -p
'''

Provide input file through -i flag.
Provide columns to use through -c flag.
Columns start at 1 being the first column in the data file.
Assumes that the second column is the amplitude f(t) at time t, given in the first column: [t,f(t)].
Show plots if the flag: -p or --plot is provided.
