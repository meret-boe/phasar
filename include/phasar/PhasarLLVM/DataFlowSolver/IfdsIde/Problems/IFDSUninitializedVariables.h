/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSUNINITIALIZEDVARIABLES_H_
#define PHASAR_PHASARLLVM_IFDSIDE_PROBLEMS_IFDSUNINITIALIZEDVARIABLES_H_

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"

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

class IFDSUninitializedVariables
    : public IFDSTabulationProblem<LLVMIFDSAnalysisDomainDefault> {
private:
  struct UninitResult {
    UninitResult() = default;
    unsigned line_nr = 0;

    std::string func_name;
    std::string file_path;

    std::string src_code;

    std::map<std::string, l_t> variableToValue;
    std::vector<std::string> var_names;
    std::vector<IFDSUninitializedVariables::n_t> ir_trace;
    bool empty() const;
    void print(std::ostream &os);
  };
      
  std::map<n_t, std::set<d_t>> UndefValueUses;

public:
  using uninit_results_t = std::map<std::string, std::map<unsigned, UninitResult>>;
  
  static void stripBottomResults(std::unordered_map<d_t, l_t> &Res);

  static const l_t TOP;
  static const l_t BOTTOM;

  IFDSUninitializedVariables(const ProjectIRDB *IRDB,
                             const LLVMTypeHierarchy *TH,
                             const LLVMBasedICFG *ICF, LLVMPointsToInfo *PT,
                             std::set<std::string> EntryPoints = {"main"});

  ~IFDSUninitializedVariables() override = default;

  FlowFunctionPtrType getNormalFlowFunction(n_t curr, n_t succ) override;

  FlowFunctionPtrType getCallFlowFunction(n_t callSite, f_t destFun) override;

  FlowFunctionPtrType getRetFlowFunction(n_t callSite, f_t calleeFun,
                                         n_t exitInst, n_t retSite) override;

  FlowFunctionPtrType getCallToRetFlowFunction(n_t callSite, n_t retSite,
                                               std::set<f_t> callees) override;

  FlowFunctionPtrType getSummaryFlowFunction(n_t callSite,
                                             f_t destFun) override;

  InitialSeeds<n_t, d_t, l_t> initialSeeds() override;

  d_t createZeroValue() const override;

  bool isZeroValue(d_t d) const override;

  void printNode(std::ostream &os, n_t n) const override;

  void printDataFlowFact(std::ostream &os, d_t d) const override;

  void printFunction(std::ostream &os, f_t m) const override;

  void emitTextReport(const SolverResults<n_t, d_t, l_t> &Results,
                      std::ostream &OS = std::cout) override;

  const std::map<n_t, std::set<d_t>> &getAllUndefUses() const;

  std::vector<UninitResult> aggregateResults();

  uninit_results_t getUnitializedResults(SolverResults<n_t, d_t, l_t> SR);
};

} // namespace psr

#endif
