/*
  LLVMTimer(modified by Sahil Yerawar)
  It has combined my loop feature extractor and Shalini's LLVMTimer
  Prerequisite passes (actual sequence): -mem2reg -loop-simplify
  Passes Order: -custom-feature-extraction -output-file=<Filename> -count-loops -insertstr -inserttimer -print-func-name
*/

#include "llvm/Transforms/Customtimer.h"

using namespace llvm;
namespace {
	namespace {
	void LoopFeaturesExtraction::printAllPath(BasicBlock*s,BasicBlock*d,Loop*L) {
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

    void LoopFeaturesExtraction::printAllPathUtil(BasicBlock*s,BasicBlock*d, std::unordered_map<BasicBlock*,bool > &visited, std::vector<BasicBlock*>&allbb,std::map<BasicBlock*, int> &maxInst,std::map<BasicBlock*, int> &parallelPath, std::map<BasicBlock*, bool> &done) {
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
void LoopFeaturesExtraction::calc_inst(Loop *L, LoopInfo *LI){
      std::vector<BasicBlock*>path1;
      LoopBlocksDFS loopBlocksDFS(L);
      LoopBlocksTraversal loopBlocksTraversal(loopBlocksDFS,LI);
      loopBlocksDFS.perform(LI);
      int total=0;
      int numBranch=0;
      int memoryInst=0;
      int callSite=0;
      for(LoopBlocksDFS::POIterator it=loopBlocksDFS.beginPostorder();it!=loopBlocksDFS.endPostorder();it++) {
          const TerminatorInst *TInst=it.operator*()->getTerminator();
          path1.push_back(*it);
          if(TInst->getNumSuccessors() == 2)
              numBranch++;
          total+=it.operator*()->getInstList().size();
          for(BasicBlock::iterator i = it.operator*()->begin() , i2 = it.operator*()->end() ; i != i2 ; i++) {
              std::string temp_str = i->getOpcodeName();
              // errs() << "operand : " << temp_str << "\n";
              if (temp_str == "load" || temp_str == "store" || temp_str == "alloca") {
                  memoryInst++;
              }
              if (isa<CallInst>(i)) {
                callSite++;
              }
          }
      }

      // errs() << "Total Instructions : " << total << "\n";
      // errs() << "Total number of branches : " << numBranch << "\n";
      // errs() << "Total Memory Instructions : " << memoryInst << "\n";
      // *out << total << ",";
      numInst.insert(std::pair<int, int>(LoopID,total));
      // *out << numBranch << ",";
      numBr.insert(std::pair<int, int>(LoopID,numBranch));
      // *out << memoryInst << ",";
      numMemOp.insert(std::pair<int, int>(LoopID,memoryInst));
      // *out << callSite <<",";
      numCallSites.insert(std::pair<int, int>(LoopID,callSite));
    }

    void LoopFeaturesExtraction::getAnalysisUsage(AnalysisUsage &AU) const{
        AU.addRequired<LoopInfoWrapperPass>();
				AU.addRequired<LoopAccessLegacyAnalysis>();
        AU.addRequired<ScalarEvolutionWrapperPass>();
        // AU.addRequired<DependenceAnalysisWrapperPass>();
        AU.setPreservesAll();
    }

    int LoopFeaturesExtraction::depth(Loop *L){
      return L->getLoopDepth();
    }

    int LoopFeaturesExtraction::runOnLoop(Loop *L,ScalarEvolution *SE, LoopInfo *LI,LoopAccessLegacyAnalysis *LAI){
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
        memoryDependencies(L,LAI);
        BasicBlock*start = L->getHeader();
        BasicBlock*dest = L->getLoopLatch();
        printAllPath(start, dest, L);
        // uniquePred(L);
        // errs()<<depth(L)<<"\n";
        noSubLoops.insert(std::pair<int, int>(LoopID, 0));
        TripMultiple.insert(std::pair<int, int>(LoopID,SE->getSmallConstantTripMultiple(L)));
        LoopDepth.insert(std::pair<int,int>(LoopID,depth(L)));
        LoopNestLevel.insert(std::pair<int,int>(LoopID,1));
				twineToLoopID.insert(std::pair<std::string,int>(Twine(L->getHeader()->getParent()->getName()).concat("|").concat(
								L->getName()).str(),LoopID));

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
        memoryDependencies(L,LAI);
        BasicBlock*start = L->getHeader();
        BasicBlock*dest = L->getLoopLatch();
        printAllPath(start, dest, L);
				twineToLoopID.insert(std::pair<std::string,int>(Twine(L->getHeader()->getParent()->getName()).concat("|").concat(
								L->getName()).str(),LoopID));
        for(Loop *subLoop : L->getSubLoops()){
          LoopID++;
          int subloopid = LoopID;
          int x = runOnLoop(subLoop,SE,LI,LAI);
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
    void LoopFeaturesExtraction::memoryDependencies(Loop *L, /*DependenceInfo *DI,*/LoopAccessLegacyAnalysis *LAI){
      // SmallVector<Instruction *, 8> memory_inst;
      // for(auto *BB : L->getBlocks()){
      //   for(Instruction &I: *BB){
      //     LoadInst *load = dyn_cast<LoadInst>(&I);
      //     StoreInst *store = dyn_cast<StoreInst>(&I);
      //     if(!store && !load)continue;
      //     memory_inst.push_back(&I);
      //   }
      // }
      // SmallVector<Instruction *, 8>::iterator I,IEnd,J,JEnd;
      // double averageMemDependence = 0.0;
      // int maxDependenceHeightMem  = 0, minDependenceHeightMem = 11111111;
      // int totalMemDependences = 0, flowDependences = 0, antiDependences = 0, outDependences = 0, inDependences = 0;
      // for (I = memory_inst.begin(), IEnd = memory_inst.end(); I!=IEnd; ++I){
      //     for (J = I, JEnd = memory_inst.end(); J!=JEnd; ++J){
      //       Instruction *Source = dyn_cast<Instruction>(*I);
      //       Instruction *Destination = dyn_cast<Instruction>(*J);
      //       if(Source == Destination)continue;
      //       if(isa<LoadInst>(Source) && isa<LoadInst>(Destination))continue;
      //       //Source->dump();
      //       //Destination->dump();
      //       auto D = DI->depends(Source, Destination, true);
      //       if(D==NULL)continue;
      //       totalMemDependences++;
      //       if(D->isInput())inDependences++;
      //       if(D->isOutput())outDependences++;
      //       if(!D->isFlow() && !D->isAnti())continue;
      //       if(D->isFlow())flowDependences++;
      //       if(D->isAnti())antiDependences++;
      //       // Not handling dependence distance for now
      //       // int distance = J-I;
      //       // if(maxDependenceHeightMem < distance){
      //       //     maxDependenceHeightMem = distance;
      //       // }
      //       // averageMemDependence += distance;
      //       // if(minDependenceHeightMem > distance) {
      //       //     minDependenceHeightMem = distance;
      //       // }
			//
      //     }
      // }
			//
      // if(totalMemDependences>0) averageMemDependence /= (1.0 * totalMemDependences);
      // else averageMemDependence = 1.0;
      // minDependence.insert(std::pair<int,int>(LoopID,minDependenceHeightMem));
      // maxDependence.insert(std::pair<int,int>(LoopID,maxDependenceHeightMem));
      // avgDependence.insert(std::pair<int,double>(LoopID,averageMemDependence));
			if ((((LAI->getInfo(L)).getDepChecker()).getDependences()) != nullptr){

				numDependence.insert(std::pair<int,int>(LoopID,(((LAI->getInfo(L)).getDepChecker()).getDependences())->size()));
			}
			else{
				numDependence.insert(std::pair<int,int>(LoopID,-1));
			}
			// inDependence.insert(std::pair<int,int>(LoopID,inDependences));
      // outDependence.insert(std::pair<int,int>(LoopID,outDependences));
      // flowDependence.insert(std::pair<int,int>(LoopID,flowDependences));
      // antiDependence.insert(std::pair<int,int>(LoopID,antiDependences));
    }

bool LoopFeaturesExtraction::runOnFunction(Function &F){
	#define DEBUG_TYPE "foo"
	DEBUG(errs()<<"First time debug mode!\n");
	#undef DEBUG_TYPE
 //LoopID = 0;
// *out<<F.getName()<<"\n";
 // DependenceInfo *DI = &getAnalysis<DependenceAnalysisWrapperPass>().getDI();
 LoopInfo *LI = &getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
 auto *LAI = &getAnalysis<LoopAccessLegacyAnalysis>();

//  LoopInfo *LI1;
//  *LI1= LI;
 ScalarEvolution *SE = &getAnalysis<ScalarEvolutionWrapperPass>().getSE();
 for (Loop *L : *LI) {
	 //  errs()<<L->getLoopDepth()<<"\n";
	 LoopID++;
	 runOnLoop(L,SE, LI,LAI);
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
bool LoopFeaturesExtraction::doFinalization(Module &M){
	// errs()<<"\0";
	// errs().flush();
	for(int i = 1;i<=LoopID;i++){
// if(LoopDepth[i] == 0 && LoopNestLevel[i] == 0 && numDependence[i] == 0 && inDependence[i] == 0 && outDependence[i] == 0 && flowDependence[i] == 0 && antiDependence[i] == 0)continue;
	// else{
	 *out<<i<<", "<<TripMultiple[i]<<", "<<LoopDepth[i]<<", "<<LoopNestLevel[i]<<", "<<TripMultiple[i]<<", "<<numInst[i]<<", "<<numMemOp[i]<<", "<<numBr[i]<<", "<<numCallSites[i]<<", "<<parallelP[i]<<", "<<maxIns[i]<<", "<<noSubLoops[i]<<", "<<numDependence[i]<<"\n";
	 out->flush();
//  }
 }
	out->close();
  errs()<<"Done feature extraction"<<"\n";
}

    char LoopFeaturesExtraction::ID = 0;
    static RegisterPass<LoopFeaturesExtraction> X("custom-feature-extraction", "Analysis pass to extract loop features from the module");

    }



    namespace {

        void CountLoops::getAnalysisUsage(AnalysisUsage &AU) const {
            AU.addRequired<LoopInfoWrapperPass>();
            AU.addRequired<ScalarEvolutionWrapperPass>();
            AU.setPreservesAll();
        }

        bool CountLoops::runOnLoop(Loop *L, LPPassManager &) {
//              if (L->getSubLoops().size() > 0) {
//                  return false;
//              }

            if(LoopType.compare("inner")==0) {
                if (L->getSubLoops().size() > 0) {
                    return false;
                }
            }
						// errs()<<"Done "<<totalLoops<<"\n";
            else if(LoopType.compare("outer")==0)
            {
                if(L->getLoopDepth()!=1)
                {
                    return false;
                }
            }
                symbolTable[totalLoops] = Twine(L->getHeader()->getParent()->getName()).concat("|").concat(
                        L->getName()).str();

                totalLoops++;
                errs()<<totalLoops<<"\n";
            return false;
        }

    }
    char CountLoops::ID = 0;
    static RegisterPass<CountLoops> XYZ("count-loops", "For Counting Loops");

    namespace {


        void InsertTimer::getAnalysisUsage(AnalysisUsage &AU) const {
            AU.addRequired<LoopInfoWrapperPass>();
            AU.addRequired<ScalarEvolutionWrapperPass>();
            AU.setPreservesAll();
        }

        bool InsertTimer::runOnLoop(Loop *L, LPPassManager &) {

//                if (L->getSubLoops().size() > 0) {
//                    return false;
//                }


            if(LoopType.compare("inner")==0) {
                if(L->getSubLoops().size()>0)
                    return false;


            }
            else if(LoopType.compare("outer")==0)
            {
                if(L->getLoopDepth()!=1)
                {
                    return false;
                }
            }


            std::string test_string = Twine(L->getHeader()->getParent()->getName()).concat("|").concat(
                    L->getName()).str();
            int ref_loopID = twineToLoopID[test_string]-1;
            // errs() << iter << ":" << symbolTable[iter] << "\n";
            errs() << iter <<":"<<symbolTable[iter]<<"\n";
            Function &F = *L->getHeader()->getParent();
            Module *mod = F.getParent();


            //Getting All Exit Blocks of Loop in exits
            SmallVector<BasicBlock *, 8> exits;
            std::vector<Instruction *> exitsFirstInst;
            L->getExitBlocks(exits);

            //iterate over each exit block and insert them in a vector

            for (auto e:exits) {
                if (e->isLandingPad())
                    continue;
                exitsFirstInst.push_back(e->getFirstNonPHI());
            }

            //adding Loop Counter
            GlobalVariable *c = mod->getGlobalVariable("cd");
            std::vector<Constant *> constIndic;
            ConstantInt *constantInt0 = ConstantInt::get(mod->getContext(), APInt(32, StringRef("0"), 10));
            ConstantInt *constantInti = ConstantInt::get(mod->getContext(), APInt(32, iter));
            constIndic.push_back(constantInt0);
            constIndic.push_back(constantInti);

            //define the global variables
            GlobalVariable *str1 = mod->getGlobalVariable(".str1");
            GlobalVariable *t1 = mod->getGlobalVariable(".t1");
            GlobalVariable *t2 = mod->getGlobalVariable(".t2");

//            float clck=CLOCKS_PER_SEC;
//            ConstantFP *constantFP = ConstantFP::get(mod->getContext(),APFloat(clck));

            //for t1:: To be inserted in the 'entry' block
            // %0 = call i64 @clock()
            // %1 = trunc i64 %0 to i32
            // store i32 %1, i32* @t1
            auto p = L->getLoopPreheader();

            CallInst *call_clock_t1 = CallInst::Create(mod->getFunction("clock"), "", p->getTerminator());
            call_clock_t1->setCallingConv(CallingConv::C);
            call_clock_t1->setTailCall(false);
            // CastInst *castInst_clock_t1 = new TruncInst(call_clock_t1,IntegerType::get(mod->getContext(),64),"",p->getTerminator());
          //  StoreInst *storeInstClockBefore = new StoreInst(call_clock_t1, t1, p->getTerminator());

            //for t2:: to be inserted in the exit blocks
            //ie. "for.end"
            for (auto &e : exitsFirstInst) {

                //adding counter
                // for.end:
                // %9 = load i64, i64* getelementptr inbounds ([300 x i64], [300 x i64]* @cd, i32 0, i32 0)
                // %10 = add i64 %9, 1
                // store i64 %10, i64* getelementptr inbounds ([300 x i64], [300 x i64]* @cd, i32 0, i32 0)
                Constant *counter = ConstantExpr::getGetElementPtr(
                        ArrayType::get(IntegerType::get(mod->getContext(), 64), totalLoops), c, constIndic);
                LoadInst *loadCounter = new LoadInst(counter, "", e);
                BinaryOperator *increment = BinaryOperator::CreateAdd(loadCounter,
                                                                      ConstantInt::get(mod->getContext(),
                                                                                       APInt(64, StringRef("1"),
                                                                                             10)), "", e);
                StoreInst *storeCounter = new StoreInst(increment, counter, e);


                //adding call to clock
                // %11 = call i64 @clock()
                // %12 = trunc i64 %11 to i32
                // store i32 %12, i32* @t2
                CallInst *call_clock_t2 = CallInst::Create(mod->getFunction("clock"), "", e);
                call_clock_t2->setCallingConv(CallingConv::C);
                call_clock_t2->setTailCall(false);
                //CastInst *castInst_clock_t2 = new TruncInst(call_clock_t2, IntegerType::get(mod->getContext(), 64), "", e);
               // StoreInst *storeInstClockAfter = new StoreInst(call_clock_t2, t2, false, e);

                // %13 = load i32, i32* @t2
                // %14 = sitofp i32 %13 to float
              //  LoadInst *loadInst_t2 = new LoadInst(t2, "", false, e);
                // CastInst *float_t2 = new SIToFPInst(loadInst_t2, Type::getFloatTy(mod->getContext()), "", e);

                // %15 = load i32, i32* @t1
                // %16 = sitofp i32 %15 to float
              //  LoadInst *loadInst_t1 = new LoadInst(t1, "", false, e);
                // CastInst *float_t1 = new SIToFPInst(loadInst_t1, Type::getFloatTy(mod->getContext()), "", e);

                // %17 = fsub float %14, %16
                // %18 = fdiv float %17, 1.000000e+06
                BinaryOperator *sub = BinaryOperator::Create(Instruction::Sub, call_clock_t2, call_clock_t1, "", e);

                GlobalVariable *t = mod->getGlobalVariable("t");
                Constant *timer_click = ConstantExpr::getGetElementPtr(
                        ArrayType::get(IntegerType::get(mod->getContext(), 64), totalLoops), t, constIndic);
                LoadInst *timer_load = new LoadInst(timer_click, "", e);
                BinaryOperator *add_ans = BinaryOperator::CreateAdd(timer_load, sub, "", e);
                StoreInst *storeTimer = new StoreInst(add_ans, timer_click, e);

            }

            iter++;
            return true;
        }




}
    char InsertTimer::ID = 0;
    static RegisterPass<InsertTimer> X("inserttimer", "Calculating time to executing Internal Loop");


    namespace {


            bool InsertStr::runOnModule(Module &M)  {


                //Declare Global Variable to Store Counter
                ArrayType *cc = ArrayType::get(IntegerType::get(M.getContext(), 64), totalLoops);
                GlobalVariable *c = new GlobalVariable(M, cc, false, GlobalValue::ExternalLinkage, 0, "cd");
                ConstantAggregateZero *const_array_2 = ConstantAggregateZero::get(cc);
                c->setInitializer(const_array_2);

                //Declare Global Variable to Store time
                ArrayType *tc = ArrayType::get(Type::getFloatTy(M.getContext()), totalLoops);
                GlobalVariable *t = new GlobalVariable(M, cc, false, GlobalValue::ExternalLinkage, 0, "t");
                ConstantAggregateZero *consgt = ConstantAggregateZero::get(cc);
                t->setInitializer(consgt);


                //Declaring Printf Format for Counter
                ArrayType *ArrayTy_0 = ArrayType::get(IntegerType::get(M.getContext(), 8), 4);
                GlobalVariable *time = new GlobalVariable(M, ArrayTy_0, true, GlobalValue::ExternalLinkage, 0,
                                                          ".counter");
                Constant *format = ConstantDataArray::getString(M.getContext(), "%d,", true);
                time->setInitializer(format);


                //Declare Global Variable for printf Format for printing Time Taken
                ArrayType *arrayType = ArrayType::get(IntegerType::get(M.getContext(), 8), 4);
                GlobalVariable *str = new GlobalVariable(M, arrayType, true, GlobalValue::ExternalLinkage, 0, ".str1");
                Constant *format2 = ConstantDataArray::getString(M.getContext(), "%d\x0A", true);
                str->setInitializer(format2);


								ArrayType *arrayType3 = ArrayType::get(IntegerType::get(M.getContext(), 8), 18);
                GlobalVariable *str3 = new GlobalVariable(M, arrayType3, true, GlobalValue::ExternalLinkage, 0, ".str3");
                Constant *format3 = ConstantDataArray::getString(M.getContext(), "\nProfiler Output\n", true);
                str3->setInitializer(format3);


                //Declare Global Variable for putting t1 before the loop and t2 after the loop
                ConstantInt *constantInt1 = ConstantInt::get(M.getContext(), APInt(64, StringRef("0"), 10));
                GlobalVariable *t1 = new GlobalVariable(M, IntegerType::get(M.getContext(), 64), false,
                                                        GlobalValue::ExternalLinkage, 0, ".t1");
                t1->setInitializer(constantInt1);
                GlobalVariable *t2 = new GlobalVariable(M, IntegerType::get(M.getContext(), 64), false,
                                                        GlobalValue::ExternalLinkage, 0, ".t2");
                t2->setInitializer(constantInt1);

                //declare clock Function
                std::vector<Type *> clock_func_args;
                FunctionType *clock_func_type = FunctionType::get(IntegerType::get(M.getContext(), 64), clock_func_args,
                                                                  false);
                Function *clock_func = Function::Create(clock_func_type, GlobalValue::ExternalLinkage, "clock", &M);
                clock_func->setCallingConv(CallingConv::C);


                //declare printf Function
                PointerType *PointerTy_6 = PointerType::get(IntegerType::get(M.getContext(), 8), 0);
                std::vector<Type *> printf_func_args;
                printf_func_args.push_back(PointerTy_6);
                FunctionType *printf_func_type = FunctionType::get(IntegerType::get(M.getContext(), 32),
                                                                   printf_func_args, true);
                Function *printf_func = Function::Create(printf_func_type, GlobalValue::ExternalLinkage, "printf", &M);
                printf_func->setCallingConv(CallingConv::C);
								errs()<<"Finished printing module"<<"\n";
            }
    }
    char InsertStr::ID = 0;
    static RegisterPass<InsertStr> Y("insertstr", "For Adding Some Global variable");


    namespace {


            bool PrintFuncName::runOnFunction(Function &F)  {
                Module *mod = F.getParent();
                Function *FF = mod->getFunction("exit");
                //Vector of Two Constant Zeroes used to access first element with
                std::vector<Constant *> zeroVector;
                //StringRef gives the string equivalent of the integer. here that is 0.
                ConstantInt *constantInt0 = ConstantInt::get(mod->getContext(),APInt(32, StringRef("0"), 10));
                zeroVector.push_back(constantInt0);               // zeroVector contains two ConstantInts, both having with value 0
                zeroVector.push_back(constantInt0);

                //Pointer to Global Counter Array name 'cd'
                GlobalVariable *c = mod->getGlobalVariable("cd");
                //Pointer to timer Array
                GlobalVariable *t = mod->getGlobalVariable("t");
                //Accessing printf Format for counter
                GlobalVariable *printfFormatCounter = mod->getGlobalVariable(".counter");
                Constant *counterFormatGEP = ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 8), 4), printfFormatCounter,zeroVector);




                //Accesing print Format for Time
                GlobalVariable *pFormat = mod->getGlobalVariable(".str1");
                Constant *timerFormatGEP = ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 8), 4), pFormat, zeroVector);




                uint64_t i = g;
                if(FF)
                {
                    for (BasicBlock &BB : F)
                    {
                        for (Instruction &e : BB)
                        {
                            if(CallInst::classof(&e))
                            {
                                CallInst *x  = dyn_cast<CallInst>(&e);
                                if(x->getCalledFunction() == FF)
                                {
                                    // errs() << "Exit instruction\n";
                                    //For all loops
																		// GlobalVariable *pFormat_init = mod->getGlobalVariable(".str3");
																		// Constant *initi_load = ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 18), 4), pFormat_init, zeroVector);
																		//


																		std::vector<Constant *> constIndicg0_init;
																		ConstantInt *constantInt0_00 = ConstantInt::get(mod->getContext(), APInt(32, 0));
																		constIndicg0_init.push_back(constantInt0);
																		constIndicg0_init.push_back(constantInt0_00);
																		std::vector<Value *> print_iniit_arg;
																		print_iniit_arg.push_back(ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 8), 18), mod->getGlobalVariable(".str3"), constIndicg0_init));
																		CallInst *call_print_init = CallInst::Create(mod->getFunction("printf"), print_iniit_arg, "", &e);
																		call_print_init->setCallingConv(CallingConv::C);
																		call_print_init->setTailCall(false);



                                    while (i < totalLoops) {




																				//for printing loop id
                                        std::vector<Value *> loop_id_argss;
                                        loop_id_argss.push_back(counterFormatGEP);
                                        loop_id_argss.push_back(ConstantInt::get(mod->getContext(), APInt(32, i)));
                                        CallInst *call_printf01 = CallInst::Create(mod->getFunction("printf"), loop_id_argss,"", &e);
                                        call_printf01->setCallingConv(CallingConv::C);
                                        call_printf01->setTailCall(false);


                                        //For Counter Purposes
                                        // %32 = load i64, i64* getelementptr inbounds ([300 x i64], [300 x i64]* @cd, i32 0, i32 0)
                                        std::vector<Constant *> constIndicg0;
                                        ConstantInt *constantInt0i = ConstantInt::get(mod->getContext(), APInt(32, i));
                                        constIndicg0.push_back(constantInt0);
                                        constIndicg0.push_back(constantInt0i);
                                        Constant *counterg0 = ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 64), totalLoops), c, constIndicg0);
                                        LoadInst *valC0 = new LoadInst(counterg0, "", &e);

                                        // %33 = call i32 (i8*, ...) @printf(i8* getelementptr inbounds ([4 x i8], [4 x i8]* @.counter, i32 0, i32 0), i64 %32)
                                        std::vector<Value *> printf_args0;
                                        printf_args0.push_back(counterFormatGEP);
                                        printf_args0.push_back(valC0);
                                        CallInst *call_printf02 = CallInst::Create(mod->getFunction("printf"), printf_args0, "", &e);
                                        call_printf02->setCallingConv(CallingConv::C);
                                        call_printf02->setTailCall(false);

                                        //For Time Differences Purposes
                                        Constant *timer = ConstantExpr::getGetElementPtr(ArrayType::get(Type::getInt64Ty(mod->getContext()), totalLoops), t, constIndicg0);
                                        LoadInst *timerLoad = new LoadInst(timer, "", &e);
                                        //CastInst *ans = new FPExtInst(timerLoad, Type::getDoubleTy(mod->getContext()), "", e);

                                        std::vector<Value *> p_args0;
                                        p_args0.push_back(timerFormatGEP);
                                        p_args0.push_back(timerLoad);
                                        CallInst *call_print_func0 = CallInst::Create(mod->getFunction("printf"), p_args0, "", &e);
                                        call_print_func0->setCallingConv(CallingConv::C);
                                        call_print_func0->setTailCall(false);
                                        i++;
                                    }

                                }
                            }
                        }
                    }
                }


                StringRef m = "main";
                if (F.getName().compare(m) == 0) {

                    for (auto &BB:F) {
                        for (auto &e:BB) {
                            if (ReturnInst::classof(&e)) {
															std::vector<Constant *> constIndicg0_init;
															ConstantInt *constantInt0_00 = ConstantInt::get(mod->getContext(), APInt(32, 0));
															constIndicg0_init.push_back(constantInt0);
															constIndicg0_init.push_back(constantInt0_00);
															std::vector<Value *> print_iniit_arg;
															print_iniit_arg.push_back(ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 8), 18), mod->getGlobalVariable(".str3"), constIndicg0_init));
															CallInst *call_print_init = CallInst::Create(mod->getFunction("printf"), print_iniit_arg, "", &e);
															call_print_init->setCallingConv(CallingConv::C);
															call_print_init->setTailCall(false);

                                while (g < totalLoops) {

                                    //for printing loop id
                                    std::vector<Value *> loop_id_args;
                                    loop_id_args.push_back(counterFormatGEP);
                                    loop_id_args.push_back(ConstantInt::get(mod->getContext(), APInt(32, g)));
                                    CallInst *call_printf = CallInst::Create(mod->getFunction("printf"), loop_id_args,
                                                                             "", &e);
                                    call_printf->setCallingConv(CallingConv::C);
                                    call_printf->setTailCall(false);

                                    //For Counter Purposes
                                    // %32 = load i64, i64* getelementptr inbounds ([300 x i64], [300 x i64]* @cd, i32 0, i32 0)
                                    std::vector<Constant *> constIndicg;
                                    ConstantInt *constantInti = ConstantInt::get(mod->getContext(), APInt(32, g));
                                    constIndicg.push_back(constantInt0);
                                    constIndicg.push_back(constantInti);
                                    Constant *counterg = ConstantExpr::getGetElementPtr(ArrayType::get(IntegerType::get(mod->getContext(), 64), totalLoops), c, constIndicg);
                                    LoadInst *valC = new LoadInst(counterg, "", &e);

                                    // %33 = call i32 (i8*, ...) @printf(i8* getelementptr inbounds ([4 x i8], [4 x i8]* @.counter, i32 0, i32 0), i64 %32)
                                    std::vector<Value *> printf_args;
                                    printf_args.push_back(counterFormatGEP);
                                    printf_args.push_back(valC);
                                    CallInst *call_printf1 = CallInst::Create(mod->getFunction("printf"), printf_args,
                                                                              "", &e);
                                    call_printf1->setCallingConv(CallingConv::C);
                                    call_printf1->setTailCall(false);

                                    //For Time Differences Purposes
                                    Constant *timer = ConstantExpr::getGetElementPtr(
                                            ArrayType::get(Type::getInt64Ty(mod->getContext()), totalLoops), t, constIndicg);
                                    LoadInst *timerLoad = new LoadInst(timer, "", &e);
                                    //CastInst *ans = new FPExtInst(timerLoad, Type::getDoubleTy(mod->getContext()), "", e);

                                    std::vector<Value *> p_args;
                                    p_args.push_back(timerFormatGEP);
                                    p_args.push_back(timerLoad);
                                    CallInst *call_print_func = CallInst::Create(mod->getFunction("printf"), p_args, "",
                                                                                 &e);
                                    call_print_func->setCallingConv(CallingConv::C);
                                    call_print_func->setTailCall(false);


                                    g++;
                                }

                                return true;
                            }
                        }
                    }
                }
                return false;
            }
    }
    char PrintFuncName::ID = 0;
    static RegisterPass<PrintFuncName> Z("print-func-name", "It prints all Function name persent into module");
}
