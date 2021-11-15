/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_IDESUMMARIES_H_
#define PHASAR_PHASARLLVM_IFDSIDE_IDESUMMARIES_H_

#include <memory>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/Utils/Table.h"

namespace psr {

template <typename N, typename D, typename F, typename V> class IDESummaries {
private:
  Table<N, D, Table<N, D, std::shared_ptr<EdgeFunction<V>>>> Summaries;

public:
  void
  addSummaries(Table<N, D, Table<N, D, std::shared_ptr<EdgeFunction<V>>>> Sum) {
    Summaries.insert(Sum);
  }
  Table<N, D, Table<N, D, std::shared_ptr<EdgeFunction<V>>>> getSummaries() {
    return Summaries;
  }
};

} // namespace psr

#endif
