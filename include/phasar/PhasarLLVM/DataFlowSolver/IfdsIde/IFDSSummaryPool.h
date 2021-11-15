/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * IFDSSummaryPool.h
 *
 *  Created on: 08.05.2017
 *      Author: philipp
 */

#ifndef PHASAR_PHASARLLVM_IFDSIDE_IFDSSUMMARYPOOL_H_
#define PHASAR_PHASARLLVM_IFDSIDE_IFDSSUMMARYPOOL_H_

#include <algorithm>
#include <iostream> // Suppress the cout as soon as to possible and get rid of this header
#include <map>
#include <set>
#include <vector>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSSummary.h"

namespace psr {

template <typename D, typename N> class IFDSSummaryPool {
private:
  /// Stores the summary that starts at a given node.
  std::map<N, std::map<std::vector<bool>, IFDSSummary<D, N>>> SummaryMap;

public:
  IFDSSummaryPool() = default;
  ~IFDSSummaryPool() = default;

  void insertSummary(N StartNode, std::vector<bool> Context,
                     IFDSSummary<D, N> Summary) {
    SummaryMap[StartNode][Context] = Summary;
  }

  bool containsSummary(N StartNode) { return SummaryMap.count(StartNode); }

  std::set<D> getSummary(N StartNode, std::vector<bool> Context) {
    return SummaryMap[StartNode][Context];
  }

  void print() {
    std::cout << "DynamicSummaries:\n";
    for (auto &Entry : SummaryMap) {
      std::cout << "Function: " << Entry.first << "\n";
      for (auto &ContextSummaries : Entry.second) {
        std::cout << "Context: ";
        for_each(ContextSummaries.first.begin(), ContextSummaries.first.end(),
                 [](bool B) { std::cout << B; });
        std::cout << "\n";
        std::cout << "Beg results:\n";
        for (auto &Result : ContextSummaries.second) {
          // result->dump();
          std::cout << "fixme\n";
        }
        std::cout << "End results!\n";
      }
    }
  }
};

} // namespace psr

#endif
