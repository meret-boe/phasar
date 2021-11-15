/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSLINEARCONSTANTANALYSIS_H_
#define PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSLINEARCONSTANTANALYSIS_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"

// Forward declaration of types for which we only use its pointer or ref type
namespace llvm {
class Instruction;
class Function;
class StructType;
class Value;
} // namespace llvm

namespace psr {

class LLVMBasedICFG;
class LLVMTypeHierarchy;
class LLVMPointsToInfo;

// A small pair data type to encode data flow facts for this LCA
struct LCAPair {
  const llvm::Value *First;
  int Second;
  LCAPair();
  LCAPair(const llvm::Value *V, int I);
  friend bool operator==(const LCAPair &Lhs, const LCAPair &Rhs);
  friend bool operator!=(const LCAPair &Lhs, const LCAPair &Rhs);
  friend bool operator<(const LCAPair &Lhs, const LCAPair &Rhs);
};

} // namespace psr

// Specialize hash to be used in containers like std::unordered_map
namespace std {
template <> struct hash<psr::LCAPair> {
  std::size_t operator()(const psr::LCAPair &K) const;
};
} // namespace std

namespace psr {

struct IFDSLinearConstantAnalysisDomain : public LLVMAnalysisDomainDefault {
  using d_t = LCAPair;
};

class IFDSLinearConstantAnalysis
    : public IFDSTabulationProblem<IFDSLinearConstantAnalysisDomain> {
public:
  IFDSLinearConstantAnalysis(const ProjectIRDB *IRDB,
                             const LLVMTypeHierarchy *TH,
                             const LLVMBasedICFG *ICF, LLVMPointsToInfo *PT,
                             std::set<std::string> EntryPoints = {"main"});

  ~IFDSLinearConstantAnalysis() override = default;

  FlowFunctionPtrType getNormalFlowFunction(n_t Curr, n_t Succ) override;

  FlowFunctionPtrType getCallFlowFunction(n_t CallSite, f_t DestFun) override;

  FlowFunctionPtrType getRetFlowFunction(n_t CallSite, f_t CalleeFun,
                                         n_t ExitInst, n_t RetSite) override;

  FlowFunctionPtrType getCallToRetFlowFunction(
      n_t CallSite, n_t RetSite,
      std::set<IFDSLinearConstantAnalysis::f_t> Callees) override;

  FlowFunctionPtrType getSummaryFlowFunction(n_t CallSite,
                                             f_t DestFun) override;

  std::map<n_t, std::set<d_t>> initialSeeds() override;

  [[nodiscard]] d_t createZeroValue() const override;

  [[nodiscard]] bool isZeroValue(d_t D) const override;

  void printNode(std::ostream &Os, n_t N) const override;

  void printDataFlowFact(std::ostream &Os, d_t D) const override;

  void printFunction(std::ostream &Os, f_t M) const override;
};

} // namespace psr

#endif
