#!/usr/local/bin/python
# parses Tandem Repeat Finder output and applies several constraints (e.g. removes redundancy in the output)
import sys
import test
import os
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import deque

repeats = deque()
filtered_repeats = deque()

def RotateMe(text,mode=0,steps=1):
    # function from http://www.how2code.co.uk/2014/05/how-to-rotate-the-characters-in-a-text-string/
    # Takes a text string and rotates
    # the characters by the number of steps.
    # mode=0 rotate right
    # mode=1 rotate left
    length=len(text)
     
    for step in range(steps):
    # repeat for required steps
        if mode==0:
            # rotate right
            text=text[length-1] + text[0:length-1]
        else:
            # rotate left
            text=text[1:length] + text[0]
    return text

def SmallestRotation(seq):
    smallest=seq
    for i in range(0,len(seq)):
        actual=RotateMe(seq,0,i)
        #print ("*" + actual)
        if (actual<smallest):
            #found new minimum
            smallest=actual
    return smallest

def lexicographicallySmallestRotation(seq):
    my_seq=Seq(seq)
    reverse_complement=my_seq.reverse_complement()
    reverse_complement=str(reverse_complement)

    smrt_seq=SmallestRotation(seq)
    smrt_rev_compl_seq=SmallestRotation(reverse_complement)

    #lexicographically smallest rotation is either one of the rotations of the sequence or its reverse complement
    if (smrt_seq < smrt_rev_compl_seq):
        return smrt_seq
    else:
        return smrt_rev_compl_seq

class Repeat:
    """Stores one line of Tandem Repeat Finder output"""
    def __init__(self, start, end, unit, repeat, cycles):
        self.start = start
        self.end = end
        self.unit = lexicographicallySmallestRotation(unit) #we want to store only lexicographically smallest rotation
        self.original = unit
        self.repeat = repeat
        self.cycles=cycles
        self.length = len(repeat)

def analyzeRepeatsFromOneSequence():
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
            if (len(r.unit)<=len(rnext.unit)):
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
            print "The repeat passed the filters."

    #add last element
    filtered_repeats.append(repeats[0])

    #output filtered repeats
    print ("**********")
    for r in filtered_repeats:
        print (r.start + " " + r.end + " " + r.unit + " " + r.original + " " + str(len(r.unit)) + " " + str(r.cycles) + " " + str(r.length))
        f.write(r.start + " " + r.end + " " + r.unit + " " + r.original + " " + str(len(r.unit)) + " " + str(r.cycles) + " " + str(r.length) + "\n")


def parseTRF(file):
    number_of_sequences=0
    for line in open(file):
        li=line.strip()
        fields=li.split()
        
        if ("Sequence:" in line):
            number_of_sequences=number_of_sequences+1
            if (number_of_sequences>1):

                if (len(repeats)>0): #tandem repeat finder must have found something
                    analyzeRepeatsFromOneSequence() #analyze repeats from the previous sequence
                    repeats.clear()
                    filtered_repeats.clear()

        #new sequence starts
        if (len(fields)==15):
            print fields
            repeats.append(Repeat(start=fields[0], end=fields[1], unit=fields[13], repeat=fields[14], cycles=fields[3]))
        else:
            print line
            print ("Skipping header.")
    #analyze last sequence
    if (len(repeats)>0):
        analyzeRepeatsFromOneSequence()
        repeats.clear()
        filtered_repeats.clear()

trf_file = sys.argv[1]
f = open(trf_file + "_output.txt","w")
f.write("start end unit original unit_length cycles length\n")
parseTRF(trf_file)
f.close()



