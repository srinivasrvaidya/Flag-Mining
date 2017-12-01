// This is an analysis pass, in LLVM v4.0 for feature extraction
// Authors Sahil Yerawar, Bhanu Prakash

#include "llvm/ADT/Statistic.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Transforms/Utils/ModuleUtils.h"
#include "llvm/Pass.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/Analysis/LoopIterator.h"
#include <llvm/IR/Instructions.h>
#include <llvm/Analysis/LoopInfo.h>
#include <llvm/Analysis/LoopAccessAnalysis.h>
#include <llvm/Analysis/ScalarEvolution.h>
#include "llvm/Analysis/DependenceAnalysis.h"
#include "llvm/IR/InstrTypes.h"
#include <llvm/MC/MCInstrDesc.h>
#include <llvm/MC/MCInst.h>
#include <math.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
// #include <pair>
using namespace llvm;
static cl::opt<std::string> outputFile ("output-file", cl::desc("Set the filename for the output of all loops features"));
namespace{
  struct LoopFeatureExtraction : public FunctionPass{
    static char ID;
    int LoopID;
    raw_fd_ostream* out;
    std::map<int, int> TripMultiple;
    std::map<int, int> parallelP;
    std::map<int, int> maxIns;
    std::map<int, int> LoopDepth;
    std::map<int, int> LoopNestLevel;
    std::map<int, int> numDependence;
    std::map<int, int> inDependence;
    std::map<int, int> outDependence;
    std::map<int, int> flowDependence;
    std::map<int, int> antiDependence;
    std::map<int, int> numInst;
    std::map<int, int> numBr;
    std::map<int, int> numDef;
    std::map<int, int> numUse;
    std::map<int, int> numMemOp;
    std::map<int, int> numCallSites;
    std::map<int, int> noSubLoops;
    std::map<int, int> uniqueInst;
    //int counter1 = 0;
    //int counter2 = 0;
    // std::map<int, int> uniquePredicate;
    // std::map<int, int> minDependence;
    // std::map<int, int> maxDependence;
    // std::map<int, double> avgDependence;
    LoopFeatureExtraction() : FunctionPass(ID) {
      std::error_code EC;
      out = new raw_fd_ostream(outputFile, EC, sys::fs::F_None);
      // LoopID=0;
    }


    void printAllPath(BasicBlock*s,BasicBlock*d,Loop*L) {
        std::vector<BasicBlock*> allbb;

        std::unordered_map<BasicBlock*,bool> visited;
        std::map<BasicBlock*, bool> done;
        std::map<BasicBlock *, int> maxInst; 
        std::map<BasicBlock *, int> parallelPath;
        for (auto&i : L->getBlocks()) {
            allbb.push_back(i);
            maxInst[i] = 0;
            parallelPath[i]=1;
            visited[i] = false;
            done[i] = false;
        }
        printAllPathUtil(s, d, visited, allbb, maxInst, parallelPath, done);

        // errs() << "Number of paths : " << parallelPath << "\n";
        // *out << parallelPath << ",";
        parallelP.insert(std::pair<int, int>(LoopID,parallelPath[s]));
        // errs() << "Max Instructions in total path : " << maxInst << "\n";
        // *out << maxInst << ",";
        maxIns.insert(std::pair<int, int>(LoopID,maxInst[s]));
    }
    
    void printAllPathUtil(BasicBlock*s,BasicBlock*d, std::unordered_map<BasicBlock*,bool > &visited, std::vector<BasicBlock*>&allbb,std::map<BasicBlock*, int> &maxInst,std::map<BasicBlock*, int> &parallelPath, std::map<BasicBlock*, bool> &done) {
        visited[s]= true;
        int parallel = 0;
        int maxy = 0;
        if(s == d) {
        	parallelPath[d] = 1;
            maxInst[d] = d->getInstList().size();
        }
        else {
            TerminatorInst *adj=s->getTerminator();
            for (unsigned int i = 0; i < adj->getNumSuccessors() ; ++i) {
                if(std::find(allbb.begin(),allbb.end(),adj->getSuccessor(i))==allbb.end()) {
                	parallel++;
                }
                else if (!visited[adj->getSuccessor(i)]) {
                    printAllPathUtil(adj->getSuccessor(i),d,visited,allbb,maxInst,parallelPath, done);
                    parallel += parallelPath[adj->getSuccessor(i)];
                    if (maxy < maxInst[adj->getSuccessor(i)]) {
                    	maxy = maxInst[adj->getSuccessor(i)];
                    }
                }
                else if (!done[adj->getSuccessor(i)]) {
        			//do nothing
        			//parallel++;
                }
                else {
                	parallel += parallelPath[adj->getSuccessor(i)];
                	if (maxy < maxInst[adj->getSuccessor(i)]) {
                    	maxy = maxInst[adj->getSuccessor(i)];
                    }
                }
            }
            parallelPath[s] = parallel;
            maxInst[s] = maxy+s->getInstList().size();
        }
        done[s]=true;
    }

    void getAnalysisUsage(AnalysisUsage &AU) const override {
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<ScalarEvolutionWrapperPass>();
        AU.addRequired<DependenceAnalysisWrapperPass>();
        AU.setPreservesAll();
    }
    int depth(Loop *L){
      return L->getLoopDepth();
    }

    int runOnLoop(Loop *L,DependenceInfo *DI, ScalarEvolution *SE, LoopInfo *LI){
      int local_LoopID = LoopID;
      //counter1 = 0, counter2=0;
      numCallSites.insert(std::make_pair(local_LoopID, 0));
      numInst.insert(std::make_pair(local_LoopID, 0));
      numBr.insert(std::make_pair(local_LoopID, 0));
      numDef.insert(std::make_pair(local_LoopID, 0));
      numUse.insert(std::make_pair(local_LoopID, 0));
      numMemOp.insert(std::make_pair(local_LoopID, 0));
      uniqueInst.insert(std::make_pair(local_LoopID, 0));
      if(L->getSubLoops().size() == 0){
        // errs()<<depth(L)<<"\n";
        // errs()<<"Loop Nest Level= "<<1<<"\n";
        memoryDependencies(L,DI);
        BasicBlock*start = L->getHeader();
        BasicBlock*dest = L->getLoopLatch();
        printAllPath(start, dest, L);
        // uniquePred(L);
        // errs()<<depth(L)<<"\n";
        noSubLoops.insert(std::pair<int, int>(LoopID, 0));
        TripMultiple.insert(std::pair<int, int>(LoopID,SE->getSmallConstantTripMultiple(L)));
        LoopDepth.insert(std::pair<int,int>(LoopID,depth(L)));
        LoopNestLevel.insert(std::pair<int,int>(LoopID,1));

        for (auto *BB : L->getBlocks()) {
          for (auto &I : *BB) {
            uniqueInst[local_LoopID]++;
            numInst[local_LoopID]++;
            if (I.mayReadOrWriteMemory()) numMemOp[local_LoopID]++;
            if (! cast<Value>(I).getType()->isVoidTy()) {
              numDef[local_LoopID]++;
              numUse[local_LoopID] += cast<Value>(I).getNumUses();
            }
            if (isa<CallInst>(I)) {
              numCallSites[local_LoopID]++;
            }
          }
          TerminatorInst *term=BB->getTerminator();
          if (term->getNumSuccessors() == 2) {
            numBr[local_LoopID]++;
          }
        }
        return 1;
      }
      else {
        // errs()<<depth(L)<<"\n";
        std::vector<int> nest_levels;
        memoryDependencies(L,DI);
        BasicBlock*start = L->getHeader();
        BasicBlock*dest = L->getLoopLatch();
        printAllPath(start, dest, L);
        for(Loop *subLoop : L->getSubLoops()){
          LoopID++;
          int subloopid = LoopID;
          int x = runOnLoop(subLoop,DI,SE,LI);
          nest_levels.push_back(x);

          numBr[local_LoopID] += TripMultiple[subloopid]*numBr[subloopid];
          numCallSites[local_LoopID] += TripMultiple[subloopid]*numCallSites[subloopid];
          numInst[local_LoopID] += TripMultiple[subloopid]*numInst[subloopid];
          numDef[local_LoopID] += TripMultiple[subloopid]*numDef[subloopid];
          numUse[local_LoopID] += TripMultiple[subloopid]*numUse[subloopid];
          numMemOp[local_LoopID] += TripMultiple[subloopid]*numMemOp[subloopid];

          for (auto *BB : subLoop->getBlocks()) {
            for (auto &I : *BB) {
              uniqueInst[local_LoopID]--;
              numInst[local_LoopID]--;
              if (I.mayReadOrWriteMemory()) numMemOp[local_LoopID]--;
              if (! cast<Value>(I).getType()->isVoidTy()) {
                numDef[local_LoopID]--;
                numUse[local_LoopID] -= cast<Value>(I).getNumUses();
              }
              if (isa<CallInst>(I)) {
                numCallSites[local_LoopID]--;
              }
            }
            TerminatorInst *term=BB->getTerminator();
            if (term->getNumSuccessors() == 2) {
              numBr[local_LoopID]--;
            }
          }
        }
        for (auto *BB : L->getBlocks()) {
          for (auto &I : *BB) {
            uniqueInst[local_LoopID]++;
            numInst[local_LoopID]++;
            if (I.mayReadOrWriteMemory()) numMemOp[local_LoopID]++;
            if (! cast<Value>(I).getType()->isVoidTy()) {
              numDef[local_LoopID]++;
              numUse[local_LoopID] += cast<Value>(I).getNumUses();
            }
            if (isa<CallInst>(I)) {
              numCallSites[local_LoopID]++;
            }
          }
          TerminatorInst *term=BB->getTerminator();
          if (term->getNumSuccessors() == 2) {
            numBr[local_LoopID]++;
          }
        }
        // uniquePred(L);
        int u = *std::max_element(nest_levels.begin(),nest_levels.end());
        // errs()<<"Loop Nest Level= "<<u+1<<"\n";
        noSubLoops.insert(std::pair<int, int>(LoopID, L->getSubLoops().size()));
        TripMultiple.insert(std::pair<int, int>(LoopID,SE->getSmallConstantTripMultiple(L)));
        LoopDepth.insert(std::pair<int,int>(local_LoopID,depth(L)));
        LoopNestLevel.insert(std::pair<int,int>(local_LoopID,u+1));
        return u+1;
      }
    }

    // Not Handling unique predicates for now
    // void uniquePred(Loop *L) {
    //   int uniPred = 0;
    //   std::vector<CmpInst::Predicate> preds;
    //   for(BasicBlock *BB: L -> getBlocks()) {
    //     for(auto &I:*BB){
    //        if (CmpInst *cmpInst = dyn_cast<CmpInst>(&I)) {
    //       if(std::find(preds.begin(),preds.end(),cmpInst->getPredicate()) == preds.end()) {
    //           uniPred++;
    //           preds.push_back(cmpInst->getPredicate());
    //           }
    //        }
    //      }
    //    }
    //    uniquePredicate.insert(std::pair<int,int>(LoopID,uniPred));
    //   //  return uniPred;
    // }

    void memoryDependencies(Loop *L, DependenceInfo *DI){
      SmallVector<Instruction *, 8> memory_inst;
      for(auto *BB : L->getBlocks()){
        for(Instruction &I: *BB){
          LoadInst *load = dyn_cast<LoadInst>(&I);
          StoreInst *store = dyn_cast<StoreInst>(&I);
          if(!store && !load)continue;
          memory_inst.push_back(&I);
        }
      }
      SmallVector<Instruction *, 8>::iterator I,IEnd,J,JEnd;
      double averageMemDependence = 0.0;
      int maxDependenceHeightMem  = 0, minDependenceHeightMem = 11111111;
      int totalMemDependences = 0, flowDependences = 0, antiDependences = 0, outDependences = 0, inDependences = 0;
      for (I = memory_inst.begin(), IEnd = memory_inst.end(); I!=IEnd; ++I){
          for (J = I, JEnd = memory_inst.end(); J!=JEnd; ++J){
            Instruction *Source = dyn_cast<Instruction>(*I);
            Instruction *Destination = dyn_cast<Instruction>(*J);
            if(Source == Destination)continue;
            if(isa<LoadInst>(Source) && isa<LoadInst>(Destination))continue;
            //Source->dump();
            //Destination->dump();
            auto D = DI->depends(Source, Destination, true);
            if(D==NULL)continue;
            totalMemDependences++;
            if(D->isInput())inDependences++;
            if(D->isOutput())outDependences++;
            if(!D->isFlow() && !D->isAnti())continue;
            if(D->isFlow())flowDependences++;
            if(D->isAnti())antiDependences++;
            // Not handling dependence distance for now
            // int distance = J-I;
            // if(maxDependenceHeightMem < distance){
            //     maxDependenceHeightMem = distance;
            // }
            // averageMemDependence += distance;
            // if(minDependenceHeightMem > distance) {
            //     minDependenceHeightMem = distance;
            // }

          }
      }

      if(totalMemDependences>0) averageMemDependence /= (1.0 * totalMemDependences);
      else averageMemDependence = 1.0;
      // minDependence.insert(std::pair<int,int>(LoopID,minDependenceHeightMem));
      // maxDependence.insert(std::pair<int,int>(LoopID,maxDependenceHeightMem));
      // avgDependence.insert(std::pair<int,double>(LoopID,averageMemDependence));
      numDependence.insert(std::pair<int,int>(LoopID,totalMemDependences));
      inDependence.insert(std::pair<int,int>(LoopID,inDependences));
      outDependence.insert(std::pair<int,int>(LoopID,outDependences));
      flowDependence.insert(std::pair<int,int>(LoopID,flowDependences));
      antiDependence.insert(std::pair<int,int>(LoopID,antiDependences));
    }
    void clearVectors(){
      // minDependence.clear();
      // maxDependence.clear();
      // avgDependence.clear();
      // uniquePredicate.clear();
      numInst.clear();
      numMemOp.clear();
      numBr.clear();
      numCallSites.clear();
      parallelP.clear();
      maxIns.clear();
      numDependence.clear();
      inDependence.clear();
      outDependence.clear();
      flowDependence.clear();
      antiDependence.clear();
      LoopDepth.clear();
      LoopNestLevel.clear();
      TripMultiple.clear();
      noSubLoops.clear();
    }
    bool runOnFunction(Function &F) override {
      #define DEBUG_TYPE "foo"
      DEBUG(errs()<<"First time debug mode!\n");
      #undef DEBUG_TYPE
     //LoopID = 0;
    // *out<<F.getName()<<"\n";
     DependenceInfo *DI = &getAnalysis<DependenceAnalysisWrapperPass>().getDI();
     LoopInfo *LI = &getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
    //  LoopInfo *LI1;
    //  *LI1= LI;
     ScalarEvolution *SE = &getAnalysis<ScalarEvolutionWrapperPass>().getSE();
     for (Loop *L : *LI) {
       //  errs()<<L->getLoopDepth()<<"\n";
       LoopID++;
       runOnLoop(L,DI,SE, LI);
     }
    //  *out<<LoopID<<"\n";
  //    for(int i = 1;i<=LoopID;i++){
	// // if(LoopDepth[i] == 0 && LoopNestLevel[i] == 0 && numDependence[i] == 0 && inDependence[i] == 0 && outDependence[i] == 0 && flowDependence[i] == 0 && antiDependence[i] == 0)continue;
  //     // else{
  //      *out<<i<<", "<<TripMultiple[i]<<", "<<LoopDepth[i]<<", "<<LoopNestLevel[i]<<", "<<TripMultiple[i]<<", "<<numInst[i]<<", "<<numMemOp[i]<<", "<<numBr[i]<<", "<<numCallSites[i]<<", "<<parallelP[i]<<", "<<maxIns[i]<<", "<<noSubLoops[i]<<", "<<inDependence[i]<<", "<<outDependence[i]<<", "<<flowDependence[i]<<", "<<antiDependence[i]<<"\n";
  //      out->flush();
  //   //  }
  //    }
    //  clearVectors();
    }
    bool doFinalization(Module &M) override{
      // errs()<<"\0";
      // errs().flush();
      for(int i = 1;i<=LoopID;i++){
	// if(LoopDepth[i] == 0 && LoopNestLevel[i] == 0 && numDependence[i] == 0 && inDependence[i] == 0 && outDependence[i] == 0 && flowDependence[i] == 0 && antiDependence[i] == 0)continue;
      // else{
       *out<<i<<", "<<TripMultiple[i]<<", "<<LoopDepth[i]<<", "<<LoopNestLevel[i]<<", "<<TripMultiple[i]<<", "<<numInst[i]<<", "<<numMemOp[i]<<", "<<numBr[i]<<", "<<numCallSites[i]<<", "<<parallelP[i]<<", "<<maxIns[i]<<", "<<noSubLoops[i]<<", "<<inDependence[i]<<", "<<outDependence[i]<<", "<<flowDependence[i]<<", "<<antiDependence[i]<<"\n";
       out->flush();
    //  }
     }
      out->close();
    }
  };

}

char LoopFeatureExtraction::ID = 0;
static RegisterPass<LoopFeatureExtraction> X("feature-extraction", "Analysis pass to extract loop features from the module");
