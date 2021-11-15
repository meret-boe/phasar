/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_EDGEVALUESET_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_EDGEVALUESET_H

#include <initializer_list>
#include <unordered_set>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/EdgeValue.h"

namespace psr {

class EdgeValueSet {
  std::unordered_set<EdgeValue> Underlying;

public:
  EdgeValueSet();
  template <typename Iter>
  EdgeValueSet(Iter Beg, Iter Ed) : Underlying(Beg, Ed) {}
  EdgeValueSet(std::initializer_list<EdgeValue> Ilist);
  auto begin() -> decltype(Underlying.begin());
  auto end() -> decltype(Underlying.end());
  auto begin() const -> decltype(Underlying.begin());
  auto end() const -> decltype(Underlying.end());
  int count(const EdgeValue &Ev) const;
  auto find(const EdgeValue &Ev) -> decltype(Underlying.find(Ev));
  auto find(const EdgeValue &Ev) const -> decltype(Underlying.find(Ev));

  size_t size() const;
  auto insert(const EdgeValue &Ev) -> decltype(Underlying.insert(Ev));
  auto insert(EdgeValue &&Ev) -> decltype(Underlying.insert(Ev));
  bool empty() const;
  bool operator==(const EdgeValueSet &Other) const;
  bool operator!=(const EdgeValueSet &Other) const;
};

} // namespace psr

#endif
