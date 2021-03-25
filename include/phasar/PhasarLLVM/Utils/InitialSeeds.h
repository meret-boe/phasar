/******************************************************************************
 * Copyright (c) 2021 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * InitialSeeds.h
 *
 *  Created on: 25.03.2021
 *      Author: Alexander Meinhold
 */

#ifndef PHASAR_PHASARLLVM_UTILS_InitialSeeds_H_
#define PHASAR_PHASARLLVM_UTILS_InitialSeeds_H_

#include <map>
#include <set>
#include <utility>
#include <memory>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/Utils/BinaryDomain.h"

namespace psr {

template <typename N, typename D, typename L = BinaryDomain> struct InitialSeeds {
public:
  InitialSeeds(std::map<N, std::set<std::pair<D, std::shared_ptr<EdgeFunction<L>>>>> Seeds) : Seeds(Seeds) {};

  InitialSeeds(std::map<N, std::set<D>> Seeds) {
    // TODO: implement
  };
private:
  std::map<N, std::set<std::pair<D, std::shared_ptr<EdgeFunction<L>>>>>
      Seeds;
};

} // namespace psr

#endif
