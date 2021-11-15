/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * IFDSSummaryGenerator.h
 *
 *  Created on: 17.05.2017
 *      Author: philipp
 */

#ifndef PHASAR_PHASARLLVM_IFDSIDE_SOLVER_IFDSSUMMARYGENERATOR_H_
#define PHASAR_PHASARLLVM_IFDSIDE_SOLVER_IFDSSUMMARYGENERATOR_H_

#include <iostream> // std::cout
#include <set>
#include <vector>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/FlowFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/LLVMZeroValue.h"
#include "phasar/PhasarLLVM/Utils/SummaryStrategy.h"

namespace psr {

template <typename N, typename D, typename F, typename I,
          typename ConcreteTabulationProblem, typename ConcreteSolver>
class IFDSSummaryGenerator {
protected:
  const F ToSummarize;
  const I Icfg;
  const SummaryGenerationStrategy CTXStrategy;

  virtual std::vector<D> getInputs() = 0;
  virtual std::vector<bool> generateBitPattern(const std::vector<D> &Inputs,
                                               const std::set<D> &Subset) = 0;

  class CTXFunctionProblem : public ConcreteTabulationProblem {
  public:
    const N Start;
    std::set<D> Facts;

    CTXFunctionProblem(N Start, std::set<D> Facts, I Icfg)
        : ConcreteTabulationProblem(Icfg), Start(Start), Facts(Facts) {
      this->solver_config.followReturnsPastSeeds = false;
      this->solver_config.autoAddZero = true;
      this->solver_config.computeValues = true;
      this->solver_config.recordEdges = false;
      this->solver_config.computePersistedSummaries = false;
    }

     std::map<N, std::set<D>> initialSeeds() override {
      std::map<N, std::set<D>> Seeds;
      Seeds.insert(make_pair(Start, Facts));
      return Seeds;
    }
  };

public:
  IFDSSummaryGenerator(F Function, I Icfg, SummaryGenerationStrategy Strategy)
      : ToSummarize(Function), Icfg(Icfg), CTXStrategy(Strategy) {}
  virtual ~IFDSSummaryGenerator() = default;
  virtual std::set<
      std::pair<std::vector<bool>, std::shared_ptr<FlowFunction<D>>>>
  generateSummaryFlowFunction() {
    std::set<std::pair<std::vector<bool>, std::shared_ptr<FlowFunction<D>>>>
        Summary;
    std::vector<D> Inputs = getInputs();
    std::set<D> Inputset;
    Inputset.insert(Inputs.begin(), Inputs.end());
    std::set<std::set<D>> InputCombinations;
    // initialize the input combinations that should be considered
    switch (CTXStrategy) {
    case SummaryGenerationStrategy::always_all:
      InputCombinations.insert(Inputset);
      break;
    case SummaryGenerationStrategy::always_none:
      InputCombinations.insert(std::set<D>());
      break;
    case SummaryGenerationStrategy::all_and_none:
      InputCombinations.insert(Inputset);
      InputCombinations.insert(std::set<D>());
      break;
    case SummaryGenerationStrategy::powerset:
      InputCombinations = computePowerSet(Inputset);
      break;
    case SummaryGenerationStrategy::all_observed:
      // TODO here we have to track what we have already observed first!
      break;
    }
    for (auto Subset : InputCombinations) {
      std::cout << "Generate summary for specific context: "
                << generateBitPattern(Inputs, Subset) << "\n";
      CTXFunctionProblem FunctionProblem(
          *Icfg.getStartPointsOf(ToSummarize).begin(), Subset, Icfg);
      ConcreteSolver Solver(FunctionProblem, true);
      Solver.solve();
      // get the result at the end of this function and
      // create a flow function from this set using the GenAll class
      std::set<N> LastInsts = Icfg.getExitPointsOf(ToSummarize);
      std::set<D> Results;
      for (auto Fact : Solver.resultsAt(*LastInsts.begin())) {
        Results.insert(Fact.first);
      }
      Summary.insert(make_pair(
          generateBitPattern(Inputs, Subset),
          std::make_shared<GenAll<D>>(Results, LLVMZeroValue::getInstance())));
    }
    return Summary;
  }
};

} // namespace psr

#endif
