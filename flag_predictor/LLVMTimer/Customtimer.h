#ifndef LLVM_INSERT_TIMER_H
#define LLVM_INSERT_TIMER_H

#include "llvm/ADT/Statistic.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Transforms/Utils/ModuleUtils.h"
#include "llvm/Pass.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Debug.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/GlobalVariable.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/IR/Instructions.h"
#include "llvm/Analysis/ScalarEvolution.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Analysis/LoopIterator.h"
#include <llvm/IR/Instructions.h>
#include <llvm/Analysis/LoopInfo.h>
#include <llvm/Analysis/LoopAccessAnalysis.h>
#include <llvm/Analysis/ScalarEvolution.h>
#include "llvm/Analysis/DependenceAnalysis.h"
#include "llvm/IR/InstrTypes.h"
#include "llvm/IR/BasicBlock.h"
#include <llvm/MC/MCInstrDesc.h>
#include <llvm/MC/MCInst.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <string>
using  namespace llvm;
namespace {
    static cl::opt<std::string>LoopType("profile-loop",cl::desc("Specify Which Loops You Want to Instrument"));
    uint64_t totalLoops=0;
    std::unordered_map<uint64_t,std::string>symbolTable;
    raw_fd_ostream* out;
    std::map<int, long long int> TripMultiple;
    std::map<int, long long int> parallelP;
    std::map<int, long long int> maxIns;
    std::map<int, long long int> LoopDepth;
    std::map<int, long long int> LoopNestLevel;
    std::map<int, long long int> numDependence;
    std::map<int, long long int> inDependence;
    std::map<int, long long int> outDependence;
    std::map<int, long long int> flowDependence;
    std::map<int, long long int> antiDependence;
    std::map<int, long long int> numInst;
    std::map<int, long long int> numBr;
    std::map<int, long long int> numDef;
    std::map<int, long long int> numUse;
    std::map<int, long long int> numMemOp;
    std::map<int, long long int> numCallSites;
    std::map<int, long long int> noSubLoops;
    std::map<int, long long int> uniqueInst;
    std::map<BasicBlock*,int> headerToLoopID;
    std::map<std::string,int> twineToLoopID;
    namespace{
	static cl::opt<std::string> outputFile ("output-file", cl::desc("Set the filename for the output of all loops features"));
  	class LoopFeaturesExtraction : public FunctionPass{
	public:
	    static char ID;
	    int LoopID;
	    LoopFeaturesExtraction() : FunctionPass(ID){
			std::error_code EC;
      			out = new raw_fd_ostream(outputFile, EC, sys::fs::F_None);
      			LoopID=0;
		}
            void printAllPath(BasicBlock*s,BasicBlock*d,Loop*L);
            void printAllPathUtil(BasicBlock*s,BasicBlock*d, std::unordered_map<BasicBlock*,bool > &visited, std::vector<BasicBlock*>&allbb,std::map<BasicBlock*, int> &maxInst,std::map<BasicBlock*, int> &parallelPath, std::map<BasicBlock*, bool> &done);
	    void calc_inst(Loop *L, LoopInfo *LI);
	    void getAnalysisUsage(AnalysisUsage &AU) const override;
	    int depth(Loop *L);
	    int runOnLoop(Loop *L,ScalarEvolution *SE, LoopInfo *LI, LoopAccessLegacyAnalysis *LAI);
	    void memoryDependencies(Loop *L,LoopAccessLegacyAnalysis *LAI);
	    bool runOnFunction(Function &F) override;
	    bool doFinalization(Module &M) override;
	};
    }

   namespace {
       class CountLoops:public LoopPass{

       public:
           static char ID;
           CountLoops():LoopPass(ID){}
           void getAnalysisUsage(AnalysisUsage&AU)const ;
           bool runOnLoop(Loop *L,LPPassManager&) override;

       };
    }

    namespace {

        struct InsertTimer : public LoopPass {

            static char ID;
            uint64_t iter = 0;

            InsertTimer() : LoopPass(ID) {}

            void getAnalysisUsage(AnalysisUsage &AU) const override;

            bool runOnLoop(Loop *L, LPPassManager &) override ;

        };

    }


    namespace {
        struct InsertStr : public ModulePass {
            static char ID;

            InsertStr() : ModulePass(ID) {}

            bool runOnModule(Module &M) override ;
        };
    }


    namespace {
        struct PrintFuncName : FunctionPass {
            static char ID;
            uint64_t g = 0;

            PrintFuncName() : FunctionPass(ID) {}


            bool runOnFunction(Function &F) override ;
        };
    }

}

#endif
