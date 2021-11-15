/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSTAINTANALYSIS_H_
#define PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSTAINTANALYSIS_H_

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"
#include "phasar/PhasarLLVM/Utils/TaintConfiguration.h"

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

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
struct HasNoConfigurationType;

/**
 * This analysis tracks data-flows through a program. Data flows from
 * dedicated source functions, which generate tainted values, into
 * dedicated sink functions. A leak is reported once a tainted value
 * reached a sink function.
 *
 * @see TaintConfiguration on how to specify your own
 * taint-sensitive source and sink functions.
 */
class IFDSTaintAnalysis
    : public IFDSTabulationProblem<LLVMAnalysisDomainDefault> {
private:
  const TaintConfiguration<const llvm::Value *> &SourceSinkFunctions;

public:
  // Setup the configuration type
  using ConfigurationTy = TaintConfiguration<const llvm::Value *>;

  /// Holds all leaks found during the analysis
  std::map<n_t, std::set<d_t>> Leaks;

  /**
   *
   * @param icfg
   * @param TSF
   * @param EntryPoints
   */
  IFDSTaintAnalysis(const ProjectIRDB *IRDB, const LLVMTypeHierarchy *TH,
                    const LLVMBasedICFG *ICF, LLVMPointsToInfo *PT,
                    const TaintConfiguration<const llvm::Value *> &TSF,
                    std::set<std::string> EntryPoints = {"main"});

  ~IFDSTaintAnalysis() override = default;

  FlowFunctionPtrType getNormalFlowFunction(n_t Curr, n_t Succ) override;

  FlowFunctionPtrType getCallFlowFunction(n_t CallSite, f_t DestFun) override;

  FlowFunctionPtrType getRetFlowFunction(n_t CallSite, f_t CalleeFun,
                                         n_t ExitInst, n_t RetSite) override;

  FlowFunctionPtrType getCallToRetFlowFunction(n_t CallSite, n_t RetSite,
                                               std::set<f_t> Callees) override;

  FlowFunctionPtrType getSummaryFlowFunction(n_t CallSite,
                                             f_t DestFun) override;

  std::map<n_t, std::set<d_t>> initialSeeds() override;

  [[nodiscard]] d_t createZeroValue() const override;

  bool isZeroValue(d_t D) const override;

  void printNode(std::ostream &Os, n_t N) const override;

  void printDataFlowFact(std::ostream &Os, d_t D) const override;

  void printFunction(std::ostream &Os, f_t M) const override;

  void emitTextReport(const SolverResults<n_t, d_t, BinaryDomain> &SR,
                      std::ostream &OS = std::cout) override;
};
} // namespace psr

#endif
