#!/bin/bash

# Take in list of bed files
# Slop 3' end of each entry by 50bp
# Intersect all entries -> If peak appears in multiple entires, call it SOM, RPM > 1 is clonal
# Should be a bedtools option to calculate how many file overlaps there are for a given region