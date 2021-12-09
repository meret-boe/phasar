/******************************************************************************
 * Copyright (c) 2019 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_WPDS_SOLVERCONFIGURATION_H_
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_WPDS_SOLVERCONFIGURATION_H_

#include <iosfwd>

#include "phasar/PhasarLLVM/DataFlowSolver/WPDS/WPDSOptions.h"

namespace psr {

struct WPDSSolverConfig {
  WPDSSolverConfig() = default;
  WPDSSolverConfig(bool RecordWitnesses, WPDSSearchDirection SearchDirection,
                   WPDSType Wpdsty);
  ~WPDSSolverConfig() = default;
  WPDSSolverConfig(const WPDSSolverConfig &) = default;
  WPDSSolverConfig &operator=(const WPDSSolverConfig &) = default;
  WPDSSolverConfig(WPDSSolverConfig &&) = default;
  WPDSSolverConfig &operator=(WPDSSolverConfig &&) = default;
  bool RecordWitnesses = false;
  WPDSSearchDirection SearchDirection = WPDSSearchDirection::FORWARD;
  WPDSType Wpdsty = WPDSType::FWPDS;
  friend std::ostream &operator<<(std::ostream &Os, const WPDSSolverConfig &Sc);
};

} // namespace psr

#endif
