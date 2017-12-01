'''
author Bhanu Prakash (cs15btech11037@iith.ac.in)
'''

import os, sys

import dill
import random
import pandas as pd
import csv
import numpy as np
from sklearn.neighbors import KNeighborsClassifier	#for library implementation

k = 1
train_folder = 'train_data/'
opt_flags = ["-loop-vectorize","-loop-unswitch","-loop-unroll","-loop-sink","-inline","-loop-reroll","-loop-reduce","-loop-load-elim","-licm","-loop-instsimplify","-loop-distribute"]

if len(sys.argv)==1:
	print "assuming k=1 in kNN classification"
else:
	if len(sys.argv)==2:
		k = int(sys.argv[1])
	else:
		print "too many parameters ..exiting"
		sys.exit()

for i in range(len(opt_flags)):
	with open(train_folder+'train_'+str(i+1)+'.csv', 'rb') as f:	#read train dataset
		reader = csv.reader(f)
		trainData = list(reader)

	trainData = [[int(x) for x in v] for v in trainData]

	X = [x[:-1] for x in trainData]
	y = [x[-1] for x in trainData]

	kNN = KNeighborsClassifier(n_neighbors=k)	#generating the model
	kNN.fit(X, y)		#train or fit the train points


	#serialize model
	dill.dump(kNN, open(train_folder+'model_'+str(i+1)+'.model', "w"))



