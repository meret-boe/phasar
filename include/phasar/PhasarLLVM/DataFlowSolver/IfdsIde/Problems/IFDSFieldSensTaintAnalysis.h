/**
 * @author Sebastian Roland <seroland86@gmail.com>
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IFDSFIELDSENSTAINTANALYSIS_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IFDSFIELDSENSTAINTANALYSIS_H

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Value.h"

#include "phasar/PhasarLLVM/ControlFlow/LLVMBasedICFG.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSFieldSensTaintAnalysis/Stats/TraceStats.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/LLVMZeroValue.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"
#include "phasar/PhasarLLVM/Domain/ExtendedValue.h"
#include "phasar/PhasarLLVM/Utils/TaintConfiguration.h"
#include "phasar/Utils/LLVMShorthands.h"

namespace llvm {
class Value;
class Function;
class StructType;
} // namespace llvm

namespace psr {

class LLVMBasedICFG;
class LLVMTypeHierarchy;
class LLVMPointsToInfo;

struct IFDSFieldSensTaintAnalysisDomain : public LLVMAnalysisDomainDefault {
  using d_t = ExtendedValue;
};

class IFDSFieldSensTaintAnalysis
    : public IFDSTabulationProblem<IFDSFieldSensTaintAnalysisDomain> {
public:
  using ConfigurationTy = TaintConfiguration<ExtendedValue>;

  IFDSFieldSensTaintAnalysis(
      const ProjectIRDB *IRDB, const LLVMTypeHierarchy *TH,
      const LLVMBasedICFG *ICF, LLVMPointsToInfo *PT,
      const TaintConfiguration<ExtendedValue> &TaintConfig,
      std::set<std::string> EntryPoints = {"main"});
  ~IFDSFieldSensTaintAnalysis() override = default;

  FlowFunctionPtrType
  getNormalFlowFunction(const llvm::Instruction *Curr,
                        const llvm::Instruction *Succ) override;

  FlowFunctionPtrType
  getCallFlowFunction(const llvm::Instruction *CallSite,
                      const llvm::Function *DestFun) override;

  FlowFunctionPtrType
  getRetFlowFunction(const llvm::Instruction *CallSite,
                     const llvm::Function *CalleeFun,
                     const llvm::Instruction *ExitInst,
                     const llvm::Instruction *RetSite) override;

  FlowFunctionPtrType
  getCallToRetFlowFunction(const llvm::Instruction *CallSite,
                           const llvm::Instruction *RetSite,
                           std::set<const llvm::Function *> Callees) override;

  FlowFunctionPtrType
  getSummaryFlowFunction(const llvm::Instruction *CallSite,
                         const llvm::Function *DestFun) override;

  std::map<const llvm::Instruction *, std::set<ExtendedValue>>
  initialSeeds() override;

  void
  emitTextReport(const SolverResults<const llvm::Instruction *, ExtendedValue,
                                     BinaryDomain> &SolverResults,
                 std::ostream &OS = std::cout) override;

  [[nodiscard]] ExtendedValue createZeroValue() const override {
    // create a special value to represent the zero value!
    return ExtendedValue(LLVMZeroValue::getInstance());
  }

  [[nodiscard]] bool isZeroValue(ExtendedValue Ev) const override {
    return LLVMZeroValue::getInstance()->isLLVMZeroValue(Ev.getValue()); //NOLINT
  }

  void printNode(std::ostream &Os, const llvm::Instruction *N) const override {
    Os << llvmIRToString(N);
  }

  void printDataFlowFact(std::ostream &Os, ExtendedValue Ev) const override {
    Os << llvmIRToString(Ev.getValue()) << "\n";
    for (const auto *const MemLocationPart : Ev.getMemLocationSeq()) {
      Os << "A:\t" << llvmIRToString(MemLocationPart) << "\n";
    }
    if (!Ev.getEndOfTaintedBlockLabel().empty()) {
      Os << "L:\t" << Ev.getEndOfTaintedBlockLabel() << "\n";
    }
    if (Ev.isVarArg()) {
      Os << "VT:\t" << Ev.isVarArgTemplate() << "\n";
      for (const auto *const VaListMemLocationPart : Ev.getVaListMemLocationSeq()) {
        Os << "VLA:\t" << llvmIRToString(VaListMemLocationPart) << "\n";
      }
      Os << "VI:\t" << Ev.getVarArgIndex() << "\n";
      Os << "CI:\t" << Ev.getCurrentVarArgIndex() << "\n";
    }
  }

  void printFunction(std::ostream &Os, const llvm::Function *M) const override {
    Os << M->getName().str();
  }

private:
  TaintConfiguration<ExtendedValue> TaintConfig;

  TraceStats TraceStats;
};

} // namespace psr

#endif // PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IFDSFIELDSENSTAINTANALYSIS_H
