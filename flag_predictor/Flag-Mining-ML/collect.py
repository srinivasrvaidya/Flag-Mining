'''
author Bhanu Prakash (cs15btech11037@iith.ac.in)
'''

import os, platform, subprocess, re, sys

#directory that contains the benchmark
source_dir = 'test/'
benchmark = ''
#cmdopt_benchmark = '~/CPUSPEC2017/cpu2017-1_0_1/benchspec/CPU/500.perlbench_r/run/run_base_refrate_mytest-m64.0000/checkspam.pl 2500 5 25 11 150 1 1 1 1'
cmdopt_benchmark = ''
build_dir = '/home/bhanu/llvm-build/build/'
lib_timer = 'lib/LLVMTimer.so'
opt_basic_flags = '-mem2reg -loop-simplify' 
opt_timer_flags = '-count-loops -insertstr -inserttimer -print-func-name'
clang = 'clang-5.0'
opt = 'opt'
llvmdis = 'llvm-dis'
#you may want to add some more flags for timer
opt_flags = ["-loop-vectorize","-loop-unswitch","-loop-unroll","-loop-sink","-inline","-loop-reroll","-loop-reduce","-loop-load-elim","-licm","-loop-instsimplify","-loop-distribute"]
remove_llbcfiles = False

train = True
train_folder = 'train_data/'
lib_fe = 'lib/LoopFeaturesExtraction.so'
fe_pass = '-feature-extraction'
useful = []
#train files are written as flag_0.csv, flag_1.csv so on

#use this or manually write the clock speed
def get_processor_name():
	if platform.system() == "Windows":
		return platform.processor()
	elif platform.system() == "Darwin":
		os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
		command ="sysctl -n machdep.cpu.brand_string"
		return subprocess.check_output(command).strip()
	elif platform.system() == "Linux":
		command = "cat /proc/cpuinfo"
		all_info = subprocess.check_output(command, shell=True).strip()
		for line in all_info.split("\n"):
			if "model name" in line:
				return re.sub( ".*model name.*:", "", line,1)
	return ""


def gettimes(file):
	global useful
	clock = None
	total = 0.0
	with open(file,"r") as f:
		times = []
		#processor clock speed (works for linux)
		processor = get_processor_name()
		clock = float([word.strip() for word in processor.split(' ')][-1][:-3])*1e9	#assuming it is GHz
		
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
				print(x)
				times.append(x)
				total+=times[-1][2]

	if file[-9:]=='_orig.txt':
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
					i+=1

	total /= clock
	return total


def add_samples(train_id, f_feat, label):
	global useful
	global train_folder
	f_w = open(train_folder+'train_'+str(train_id)+'.csv','a')
	with open(f_feat,"r") as f:
		##assuming no header
		i=0
		j=0
		for line in f:
			if j<len(useful) and i==useful[j]:
				#add training sample
				f_w.write(line[:-1]+', '+str(label)+'\n')
				j+=1
			i+=1


def createExec():
	global source_dir
	global benchmark
	global opt_flags
	global opt_basic_flags
	global opt_timer_flags
	global opt
	global clang
	global llvmdis
	global remove_llbcfiles
	global build_dir
	global lib_fe
	global fe_pass
	global lib_timer
	new_dir = source_dir+"temp/"
	os.system("mkdir "+new_dir)
	os.system(opt+" "+opt_basic_flags+" "+source_dir+benchmark+".ll -o "+new_dir+benchmark+"_base.bc")
	#below call maynot be required, if so use _base.bc for further steps
	os.system(llvmdis+" "+new_dir+benchmark+"_base.bc -o "+new_dir+benchmark+"_base.ll")

	#find features
	os.system(opt+' -load '+build_dir+lib_fe+' '+fe_pass+' '+new_dir+benchmark+"_base.ll")
	##this features are saved as Output.txt in current directory
	##later this features are associated with flags, to get the training samples


	#below is for creating copies of new ll files with opt flags
	os.system(opt+" -load "+build_dir+lib_timer+" "+opt_timer_flags+" "+new_dir+benchmark+"_base.ll"+" -o "+new_dir+benchmark+"_orig.bc")
	os.system(clang+" "+new_dir+benchmark+"_orig.bc -o "+new_dir+benchmark+"_orig.out")
	for i in range(len(opt_flags)):
		os.system(opt+" -load "+build_dir+lib_timer+" "+opt_flags[i]+" "+opt_timer_flags+" "+new_dir+benchmark+"_base.ll"+" -o "+new_dir+benchmark+"_"+str(i+1)+".bc")
		os.system(clang+" "+new_dir+benchmark+"_"+str(i+1)+".bc -o "+new_dir+benchmark+"_"+str(i+1)+".out")

	if remove_llbcfiles:
		os.system("rm "+new_dir+"*.bc")
		os.system("rm "+new_dir+benchmark+"_base.ll")


def execAll():
	global source_dir
	global benchmark
	global opt_flags
	global cmdopt_benchmark
	global useful
	new_dir = source_dir+"temp/"
	#below line needs some modification (executing the binary)
	os.system("./"+new_dir+benchmark+"_orig.out "+cmdopt_benchmark+" > "+new_dir+benchmark+"_orig.txt")
	#assuming the llvm timer output is in new_dir/benchmark_orig.txt
	def_time = gettimes(new_dir+benchmark+"_orig.txt")


	print useful
	zzz = raw_input("continue...")

	times = []
	for i in range(len(opt_flags)):
		#below line needs some modification (executing the binary)
		os.system("./"+new_dir+benchmark+"_"+str(i+1)+".out "+cmdopt_benchmark+" > "+new_dir+benchmark+"_"+str(i+1)+".txt")
		#assuming the llvm timer output is in new_dir/benchmark_<optnum>.txt
		times.append((gettimes(new_dir+benchmark+"_"+str(i+1)+".txt"),i))

	# print(times)
	times.sort(key=lambda x: x[0])
	print "\n\n****Flag Mining - Result****\n"
	print "default time = "+str(def_time*1000)+" ms"
	print "\nGood Flags:- \n"

	if def_time != 0:
		for i in range(len(times)):
			if times[i][0] <= def_time*0.9999:
				print opt_flags[times[i][1]]+" - "+str(times[i][0]*1000)+" ms"
				####also add new training_samples
				'''here it goes'''
				add_samples(i+1,'Output.txt',1)
			else:
				add_samples(i+1,'Output.txt',0)


		#printing the top 5 flags
		print "\n"
		print "top 5 flags:- \n"
		for i in range(5):
			sys.stdout.write(opt_flags[times[i][1]]+" ")
		print ''
	else:
		print "Default time is already 0!"


if len(sys.argv)==1:
	print "one more parameter benchmark needed"
else:
	if len(sys.argv)==2:
		benchmark = sys.argv[1]
	else:
		print "too many parameters ..exiting"
		sys.exit()

createExec()
execAll()
#remove the temp directory
if remove_llbcfiles:
	os.system('rm '+source_dir+'temp/*')
	os.system('rmdir '+source_dir+'temp')
