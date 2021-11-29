/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * PathEdge.h
 *
 *  Created on: 16.08.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SOLVER_PATHEDGE_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SOLVER_PATHEDGE_H

#include <ostream>

namespace psr {

template <typename N, typename D> class PathEdge {
private:
  const N Target;
  const D DSource;
  const D DTarget;

public:
  PathEdge(D DSource, N Target, D DTarget)
      : Target(Target), DSource(DSource), DTarget(DTarget) {}

  ~PathEdge() = default;

  PathEdge(const PathEdge &) = default;

  PathEdge &operator=(const PathEdge &) = default;

  PathEdge(PathEdge &&)  noexcept = default;

  PathEdge &operator=(PathEdge &&)  noexcept = default;

  [[nodiscard]] N getTarget() const { return Target; }

  [[nodiscard]] D factAtSource() const { return DSource; }

  [[nodiscard]] D factAtTarget() const { return DTarget; }

  friend std::ostream &operator<<(std::ostream &Os, const PathEdge &PathEdge) {
    return Os << "<" << PathEdge.DSource << "> -> <" << PathEdge.Target << ","
              << PathEdge.DTarget << ">";
  }
};

} // namespace psr

#endif
