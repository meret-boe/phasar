/******************************************************************************
 * Copyright (c) 2018 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDETYPESTATEANALYSIS_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDETYPESTATEANALYSIS_H

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

#include "llvm/IR/InstrTypes.h"

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctionComposer.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IDETabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/TypeStateDescriptions/TypeStateDescription.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"

namespace llvm {
class CallBase;
class Instruction;
class Function;
class Value;
} // namespace llvm

namespace psr {

class LLVMBasedICFG;
class LLVMTypeHierarchy;
class LLVMPointsToInfo;

struct IDETypeStateAnalysisDomain : public LLVMAnalysisDomainDefault {
  using l_t = int;
};

class IDETypeStateAnalysis
    : public IDETabulationProblem<IDETypeStateAnalysisDomain> {
public:
  using IDETabProblemType = IDETabulationProblem<IDETypeStateAnalysisDomain>;
  using typename IDETabProblemType::d_t;
  using typename IDETabProblemType::f_t;
  using typename IDETabProblemType::i_t;
  using typename IDETabProblemType::l_t;
  using typename IDETabProblemType::n_t;
  using typename IDETabProblemType::t_t;
  using typename IDETabProblemType::v_t;

  using ConfigurationTy = TypeStateDescription;

private:
  const TypeStateDescription &TSD;
  std::map<const llvm::Value *, std::unordered_set<const llvm::Value *>>
      PointsToCache;
  std::map<const llvm::Value *, std::set<const llvm::Value *>>
      RelevantAllocaCache;

  /**
   * @brief Returns all alloca's that are (indirect) aliases of V.
   *
   * Currently PhASAR's points-to information does not include alloca
   * instructions, since alloca instructions, i.e. memory locations, are of
   * type T* for a target type T. Thus they do not alias directly. Therefore,
   * for each alias of V we collect related alloca instructions by checking
   * load and store instructions for used alloca's.
   */
  std::set<d_t> getRelevantAllocas(d_t V);

  /**
   * @brief Returns whole-module aliases of V.
   *
   * This function retrieves whole-module points-to information. We store
   * already computed points-to information in a cache to prevent expensive
   * recomputation since the whole module points-to graph can be huge. This
   * might become unnecessary once PhASAR's PointsToGraph starts using a cache
   * itself.
   */
  std::set<d_t> getWMPointsToSet(d_t V);

  /**
   * @brief Provides whole module aliases and relevant alloca's of V.
   */
  std::set<d_t> getWMAliasesAndAllocas(d_t V);

  /**
   * @brief Provides local aliases and relevant alloca's of V.
   */
  std::set<d_t> getLocalAliasesAndAllocas(d_t V, const std::string &Fname);

  /**
   * @brief Checks if the type machtes the type of interest.
   */
  bool hasMatchingType(d_t V);

public:
  const l_t TOP;
  const l_t BOTTOM;

  IDETypeStateAnalysis(const ProjectIRDB *IRDB, const LLVMTypeHierarchy *TH,
                       const LLVMBasedICFG *ICF, LLVMPointsToInfo *PT,
                       const TypeStateDescription &TSD,
                       std::set<std::string> EntryPoints = {"main"});

  ~IDETypeStateAnalysis() override = default;

  // start formulating our analysis by specifying the parts required for IFDS

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

  // in addition provide specifications for the IDE parts

  std::shared_ptr<EdgeFunction<l_t>>
  getNormalEdgeFunction(n_t Curr, d_t CurrNode, n_t Succ,
                        d_t SuccNode) override;

  std::shared_ptr<EdgeFunction<l_t>>
  getCallEdgeFunction(n_t CallSite, d_t SrcNode, f_t DestinationFunction,
                      d_t DestNode) override;

  std::shared_ptr<EdgeFunction<l_t>>
  getReturnEdgeFunction(n_t CallSite, f_t CalleeFunction, n_t ExitInst,
                        d_t ExitNode, n_t ReSite, d_t RetNode) override;

  std::shared_ptr<EdgeFunction<l_t>>
  getCallToRetEdgeFunction(n_t CallSite, d_t CallNode, n_t RetSite,
                           d_t RetSiteNode, std::set<f_t> Callees) override;

  std::shared_ptr<EdgeFunction<l_t>>
  getSummaryEdgeFunction(n_t CallSite, d_t CallNode, n_t RetSite,
                         d_t RetSiteNode) override;

  l_t topElement() override;

  l_t bottomElement() override;

  /**
   * We have a lattice with BOTTOM representing all information
   * and TOP representing no information. The other lattice elements
   * are defined by the type state description, i.e. represented by the
   * states of the finite state machine.
   *
   * @note Only one-level lattice's are handled currently
   */
  l_t join(l_t Lhs, l_t Rhs) override;

  std::shared_ptr<EdgeFunction<l_t>> allTopFunction() override;

  void printNode(std::ostream &Os, n_t N) const override;

  void printDataFlowFact(std::ostream &Os, d_t D) const override;

  void printFunction(std::ostream &Os, f_t M) const override;

  void printEdgeFact(std::ostream &Os, l_t L) const override;

  void emitTextReport(const SolverResults<n_t, d_t, l_t> &SR,
                      std::ostream &OS = std::cout) override;

  // customize the edge function composer
  class TSEdgeFunctionComposer : public EdgeFunctionComposer<l_t> {
  private:
    l_t BotElement;

  public:
    TSEdgeFunctionComposer(std::shared_ptr<EdgeFunction<l_t>> F,
                           std::shared_ptr<EdgeFunction<l_t>> G, l_t Bot)
        : EdgeFunctionComposer<l_t>(std::move(F), std::move(G)), BotElement(Bot){};
    std::shared_ptr<EdgeFunction<l_t>>
    joinWith(std::shared_ptr<EdgeFunction<l_t>> OtherFunction) override;
  };

  class TSEdgeFunction : public EdgeFunction<l_t>,
                         public std::enable_shared_from_this<TSEdgeFunction> {
  protected:
    const TypeStateDescription &TSD;
    // Do not use a reference here, since LLVM's StringRef's (obtained by str())
    // might turn to nullptr for whatever reason...
    const std::string Token;
    const llvm::CallBase *CallSite;

  public:
    TSEdgeFunction(const TypeStateDescription &Tsd, std::string  Tok,
                   const llvm::CallBase *Cb)
        : TSD(Tsd), Token(std::move(Tok)), CallSite(Cb){};

    l_t computeTarget(l_t Source) override;

    std::shared_ptr<EdgeFunction<l_t>>
    composeWith(std::shared_ptr<EdgeFunction<l_t>> SecondFunction) override;

    std::shared_ptr<EdgeFunction<l_t>>
    joinWith(std::shared_ptr<EdgeFunction<l_t>> OtherFunction) override;

    bool equalTo(std::shared_ptr<EdgeFunction<l_t>> Other) const override;

    void print(std::ostream &OS, bool IsForDebug = false) const override;
  };
  class TSConstant : public EdgeFunction<l_t>,
                     public std::enable_shared_from_this<TSConstant> {
    const TypeStateDescription &TSD;
    l_t State;

  public:
    TSConstant(const TypeStateDescription &TSD, l_t State);

    l_t computeTarget(l_t Source) override;

    std::shared_ptr<EdgeFunction<l_t>>
    composeWith(std::shared_ptr<EdgeFunction<l_t>> SecondFunction) override;

    std::shared_ptr<EdgeFunction<l_t>>
    joinWith(std::shared_ptr<EdgeFunction<l_t>> OtherFunction) override;

    bool equalTo(std::shared_ptr<EdgeFunction<l_t>> Other) const override;

    void print(std::ostream &OS, bool IsForDebug = false) const override;
  };
};

} // namespace psr

#endif
