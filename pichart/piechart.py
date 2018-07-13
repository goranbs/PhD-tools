#!/bin/python
"""
===============
Basic pie chart
===============

Demo of a basic pie chart plus a few additional features.

In addition to the basic pie chart, this demo shows a few optional features:

    * slice labels
    * auto-labeling the percentage
    * offsetting a slice with "explode"
    * drop-shadow
    * custom start angle

Note about the custom start angle:

The default ``startangle`` is 0, which would start the "Frogs" slice on the
positive x-axis. This example sets ``startangle = 90`` such that everything is
rotated counter-clockwise by 90 degrees, and the frog slice starts on the
positive y-axis.
"""

import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'size'   : 14}
        
wedgeprops = {'linewidth' : 0.0
                            }
                            
plt.rc('font', **font)

# Pie chart, where the slices will be ordered and plotted counter-clockwise:
labels = 'Frogs', 'Hogs', 'Dogs', 'Logs', 'Cats', 'Socks'
colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'paleturquoise', 'wheat']
sizes = [5, 12, 18, 15, 20, 30]
explode = (-0.03, -0.02, -0.01, 0, 0.1, 0.2)  # only "explode" the 2nd slice (i.e. 'Hogs')

fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, labeldistance=1.1, autopct='%1.0f%%', pctdistance=0.8,
        shadow=True, startangle=90, colors=colors, wedgeprops=wedgeprops)
ax1.set_alpha(0.7)        
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.savefig("Piechart.eps",format="eps")
plt.show()