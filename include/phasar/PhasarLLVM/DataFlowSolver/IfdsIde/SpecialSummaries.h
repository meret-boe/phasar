/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * SpecialSummaries.h
 *
 *  Created on: 05.05.2017
 *      Author: philipp
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SPECIALSUMMARIES_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SPECIALSUMMARIES_H

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "phasar/Config/Configuration.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/FlowFunctions.h"
#include "phasar/PhasarLLVM/Utils/BinaryDomain.h"
#include "phasar/Utils/IO.h" // readFile

#include "llvm/IR/Function.h"

namespace psr {

template <typename D, typename V = BinaryDomain> class SpecialSummaries {
  using FlowFunctionType = FlowFunction<D>;
  using FlowFunctionPtrType = std::shared_ptr<FlowFunction<D>>;

private:
  std::map<std::string, FlowFunctionPtrType> SpecialFlowFunctions;
  std::map<std::string, std::shared_ptr<EdgeFunction<V>>> SpecialEdgeFunctions;
  std::vector<std::string> SpecialFunctionNames;

  // Constructs the SpecialSummaryMap such that it contains all glibc,
  // llvm.intrinsics and C++'s new, new[], delete, delete[] with identity
  // flow functions.
  SpecialSummaries() {
    // insert default flow and edge functions
    for (auto FunctionName : //NOLINT
         PhasarConfig::getPhasarConfig().specialFunctionNames()) {
      SpecialFlowFunctions.insert(
          std::make_pair(FunctionName, Identity<D>::getInstance()));
      SpecialEdgeFunctions.insert(
          std::make_pair(FunctionName, EdgeIdentity<V>::getInstance()));
    }
  }

public:
  SpecialSummaries(const SpecialSummaries &) = delete;
  SpecialSummaries &operator=(const SpecialSummaries &) = delete;
  SpecialSummaries(SpecialSummaries &&) = delete;
  SpecialSummaries &operator=(SpecialSummaries &&) = delete;
  ~SpecialSummaries() = default;

  static SpecialSummaries<D, V> &getInstance() {
    static SpecialSummaries<D, V> Instance;
    return Instance;
  }

  // Returns true, when an existing function is overwritten, false otherwise.
  bool provideSpecialSummary(const std::string &Name,
                             FlowFunctionPtrType Flowfunction) {
    bool Override = containsSpecialSummary(Name);
    SpecialFlowFunctions[Name] = Flowfunction;
    return Override;
  }

  // Returns true, when an existing function is overwritten, false otherwise.
  bool provideSpecialSummary(const std::string &Name,
                             FlowFunctionPtrType Flowfunction,
                             std::shared_ptr<EdgeFunction<V>> Edgefunction) {
    bool Override = containsSpecialSummary(Name);
    SpecialFlowFunctions[Name] = Flowfunction;
    SpecialEdgeFunctions[Name] = Edgefunction;
    return Override;
  }

  bool containsSpecialSummary(const llvm::Function *Function) {
    return containsSpecialSummary(Function->getName().str());
  }

  bool containsSpecialSummary(const std::string &Name) {
    return SpecialFlowFunctions.count(Name);
  }

  FlowFunctionPtrType
  getSpecialFlowFunctionSummary(const llvm::Function *Function) {
    return getSpecialFlowFunctionSummary(Function->getName().str());
  }

  FlowFunctionPtrType getSpecialFlowFunctionSummary(const std::string &Name) {
    return SpecialFlowFunctions[Name];
  }

  std::shared_ptr<EdgeFunction<V>>
  getSpecialEdgeFunctionSummary(const llvm::Function *Function) {
    return getSpecialEdgeFunctionSummary(Function->getName().str());
  }

  std::shared_ptr<EdgeFunction<V>>
  getSpecialEdgeFunctionSummary(const std::string &Name) {
    return SpecialEdgeFunctions[Name];
  }

  friend std::ostream &operator<<(std::ostream &Os,
                                  const SpecialSummaries<D> &Ss) {
    Os << "SpecialSummaries:\n";
    for (auto &Entry : Ss.SpecialFunctionNames) {
      Os << Entry.first << " ";
    }
    return Os;
  }
};
} // namespace psr

#endif
