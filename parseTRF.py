#!/usr/local/bin/python
# parses Tandem Repeat Finder output and applies several constraints (e.g. removes redundancy in the output)
import sys
import test
import os
import re
import subprocess
from Bio import SeqIO
from collections import deque

repeats = deque()
filtered_repeats = deque()

class Repeat:
    """Stores one line of Tandem Repeat Finder output"""
    def __init__(self, start, end, unit, repeat, cycles):
        self.start = start
        self.end = end
        self.unit = unit
        self.repeat = repeat
        self.cycles=cycles
        self.length = len(repeat)

def parseTRF(file):
    for line in open(file):
        li=line.strip()
        fields=li.split()
        
        if (len(fields)==15):
            print fields
            repeats.append(Repeat(start=fields[0], end=fields[1], unit=fields[13], repeat=fields[14], cycles=fields[3]))
        else:
            print ("Skipping header.")

trf_file = sys.argv[1]
parseTRF(trf_file)

while (len(repeats)>1):
    for j in repeats:
        sys.stdout.write(j.unit + " ")
    print ("")

    r=repeats.popleft()
    rnext=repeats.popleft()
    print ("---")
    print (r.start + " " + r.end + " " + r.unit)
    print (rnext.start + " " + rnext.end + " " + rnext.unit)

    removal=False
    keepFirst=False #we want to remove first element by default; unless second element has bigger unit and same coordinates
    keepSecond=True

    #Rule1 if start and end are the same, keep only repeat with smaller unit size
    if ((r.start==rnext.start) and (r.end==rnext.end)):
        print ("Same.")
        if (len(r.unit)<len(rnext.unit)):
            #second one has the same unit (multiplied by some coefficient and hence should be removed
            keepSecond=False
            keepFirst=True
            removal=True
    else:
        #Rule2 if start and end are close to each ther (one unit apart), keep only repeat with smaller unit size
        distance=abs(int(rnext.start)-int(r.start))
        minimal_unit=min(len(r.unit),len(rnext.unit))
        #print ("distance: " + str(distance) + " minimal_unit: " + str(minimal_unit))

        if (distance<=minimal_unit):
            print ("Repeats too close to each other.")
            #repeats are too close to each other, the one with larger unit should be removed
            if(len(r.unit)<len(rnext.unit)):
                #remove second
                print ("remove second")
                keepSecond=False
                keepFirst=True
                removal=True
            else:
                print ("remove first")
                keepSecond=True
                keepFirst=False
                removal=True

    if (keepSecond==True):
        #append to the left so that it's at the beginning for the next iteration
        repeats.appendleft(rnext)
    if (keepFirst==True):
        #needs to stay in the loop
        repeats.appendleft(r)

    if (removal==False):
        #if nothing was removed in this cycle of the loop, that means we can safely add element on the very left to the filtered set
        filtered_repeats.append(r)

#add last element
filtered_repeats.append(repeats[0])

#output filtered repeats
print ("**********")
for r in filtered_repeats:
    print (r.start + " " + r.end + " " + r.unit + " " + str(r.length))


