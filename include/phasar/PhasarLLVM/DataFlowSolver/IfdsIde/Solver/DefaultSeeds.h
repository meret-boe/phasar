/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * DefaultSeeds.h
 *
 *  Created on: 03.11.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_PHASARLLVM_IFDSIDE_SOLVER_DEFAULTSEEDS_H_
#define PHASAR_PHASARLLVM_IFDSIDE_SOLVER_DEFAULTSEEDS_H_

#include <map>
#include <set>
#include <vector>

namespace psr {

class DefaultSeeds {
public:
  template <typename N, typename D>
  static std::map<N, std::set<D>> make(std::vector<N> Instructions,
                                       D ZeroNode) {
    std::map<N, std::set<D>> Res;
    for (const N &Instruction : Instructions) {
      Res.insert({Instruction, std::set<D>{ZeroNode}});
}
    return Res;
  }
};
} // namespace psr

#endif
