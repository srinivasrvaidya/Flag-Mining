'''
author Bhanu Prakash (cs15btech11037@iith.ac.in)
'''

import os, platform, subprocess, re, sys

import dill
import random
import pandas as pd
import csv
import numpy as np
from sklearn.neighbors import KNeighborsClassifier	#for library implementation

source_dir = 'test/'
benchmark = ''
#cmdopt_benchmark = '~/CPUSPEC2017/cpu2017-1_0_1/benchspec/CPU/500.perlbench_r/run/run_base_refrate_mytest-m64.0000/checkspam.pl 2500 5 25 11 150 1 1 1 1'
cmdopt_benchmark = ''
build_dir = '/home/bhanu/llvm-build/build/'
lib_timer = 'lib/LLVMTimer.so'


train_folder = 'train_data/'
opt_flags = ["-loop-vectorize","-loop-unswitch","-loop-unroll","-loop-sink","-inline","-loop-reroll","-loop-reduce","-loop-load-elim","-licm","-loop-instsimplify","-loop-distribute"]

clang = 'clang-5.0'
llvmdis = 'llvm-dis'
opt = 'opt'

opt_basic_flags = '-mem2reg -loop-simplify'
opt_timer_flags = '-count-loops -insertstr -inserttimer -print-func-name'

lib_fe = 'lib/LoopFeaturesExtraction.so'
fe_pass = '-feature-extraction'
useful = []
weights = []
flag_weights = []

def find_useful(file):
	global useful
	total = 0.0
	with open(file,"r") as f:
		frontier_flag = False # Sahil: This is the flag which tells whether to ignore or to analyze the line
		for line in f:
			'''
			Sahil: I have modified the LLVMTimer  such that it insert the line "Profiler output"
			before the timer analysis, we check whether we get this line, if yes, flip
			the boolean flag. Otherwise, until we reach there, we ignore the line.
			'''
			if "Profiler Output" in line:
				frontier_flag = True
				continue
			if not frontier_flag or line == "\n":
				continue
			else:
				x = [float(word.strip()) for word in line.split(',')]
				total+=x[2]

	if total == 0:
		return total

	with open(file,"r") as f:
		frontier_flag = False # Sahil: This is the flag which tells whether to ignore or to analyze the line
		i=0
		for line in f:
			if "Profiler Output" in line:
				frontier_flag = True
				continue
			if not frontier_flag or line == "\n":
				continue
			else:
				x = [float(word.strip()) for word in line.split(',')]
				if x[2]>0.05*total:
					useful.append(i)
					weights.append(x[2]/total)
				i+=1
	return total


def predict_flags(file):
	global useful
	global train_folder
	global weights
	global flag_weights
	
	testData = []
	with open(file,"r") as f:
		##assuming no header
		i=0
		j=0
		for line in f:
			if j<len(useful) and i==useful[j]:
				#add loop sample to testData
				testData.append([float(k) for k in line[:-1].split(',')])
				j+=1
			i+=1

	preds = []
	models = []
	for i in range(len(opt_flags)):
		with open(train_folder+'model_'+str(i+1)+'.model') as f:
			models.append(dill.load(f))

	for j in range(len(opt_flags)):
		preds.append(models[j].predict(testData))

	for j in range(len(opt_flags)):
		flag_weights.append(0.0)
	for i in range(len(preds)):
		for j in range(len(preds[0])):
			if preds[i][j] == 1:
				flag_weights[i] += weights[j]

	w = [(flag_weights[i],i) for i in range(len(flag_weights))]
	w.sort(key=lambda x: x[0],reverse=True)

	x = []
	for i in range(len(w)):
		if w[i][0] > 0.0:
			x.append(w[i][1])
		else:
			break
	return x


if len(sys.argv)==1:
	print "need a program (.ll file) to predict flags"
else:
	if len(sys.argv)==2:
		benchmark = sys.argv[1]
	else:
		print "too many parameters ..exiting"
		sys.exit()

new_dir = source_dir+"temp/"
os.system("mkdir "+new_dir)
os.system(opt+" "+opt_basic_flags+" "+source_dir+benchmark+".ll -o "+new_dir+benchmark+"_base.bc")
os.system(llvmdis+" "+new_dir+benchmark+"_base.bc -o "+new_dir+benchmark+"_base.ll")

#find features
os.system(opt+' -load '+build_dir+lib_fe+' '+fe_pass+' '+new_dir+benchmark+"_base.ll")
##this features are saved as Output.txt in current directory
##later this features are used to predict the flags

os.system(opt+" -load "+build_dir+lib_timer+" "+opt_timer_flags+" "+new_dir+benchmark+"_base.ll"+" -o "+new_dir+benchmark+"_orig.bc")
os.system(clang+" "+new_dir+benchmark+"_orig.bc -o "+new_dir+benchmark+"_orig.out")

#remove unnecessary files
os.system('rm '+new_dir+'*.bc')

#execute default
os.system("./"+new_dir+benchmark+"_orig.out "+cmdopt_benchmark+" > "+new_dir+benchmark+"_orig.txt")

total = find_useful(new_dir+benchmark+"_orig.txt")

if total !=0:
	preds = predict_flags('Output.txt')

	print flag_weights

	print "Predicted Flags (in the order of importance):-"
	for i in range(len(preds)):
		print opt_flags[preds[i]]
else:
	print "Error! default run-time is zero"

print ''