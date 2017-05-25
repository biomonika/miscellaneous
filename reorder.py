import sys
import os
import os.path
from os.path import basename
import itertools
from itertools import izip
import re
import collections
from collections import defaultdict


######################################
#python reorder.py motifFileWithIPds motifFileWithErrors
motifFileWithIPds=str(sys.argv[1])
motifFileWithErrors=str(sys.argv[2])
######################################

print("motifFileWithIPds: " + motifFileWithIPds + "; " + "motifFileWithErrors: " + motifFileWithErrors)

#print motifFile
errorDict=defaultdict()

outputEmpty=motifFileWithErrors.replace("merged_","ordered_")
outputEmpty=basename(outputEmpty)
outputEmpty=("test" + outputEmpty)

if (os.path.exists(outputEmpty)):
	print("File " + outputEmpty + " already exists. Quit.")
	sys.exit()

fout = open(outputEmpty, 'w') #write results here

def formatWhitespaces(text):
	text=text.rstrip() #remove newline
	text=re.sub(' +',' ',text) #remove multiple whitespaces
	text=re.sub('^ ','',text) #remove leading whitespace
	return text

with open(motifFileWithErrors) as f:
	for line in f:
		array=line.rstrip().replace("\t", " ").split(" ") #replace tabs by spaces and split
		#print array
		key=str(array[0:3])
		values=array[3:]
		errorDict[key]=values
f.close()

#print("printing dictionary")
#for k in errorDict.keys():
#	print k, errorDict[k]

with open(motifFileWithIPds) as f:
	i=0
	for line in f:
		i=i+1 #line number
		if (i>30000):
			fout.close()
			quit() #file too large, let's quit
		array=line.rstrip().replace("\t", " ").split(" ") #replace tabs by spaces and split
		key=array[0:3]
		if str(key) in errorDict.keys():
			values=errorDict[str(key)]
		else:
			values=[' ',' ',' ',' ']
		#print values
		merged=key+values
		merged=reduce(lambda key, values: key+" "+values, merged)
		#print(str(merged))
		#print(merged)
		fout.write(merged+"\n")
fout.close()

