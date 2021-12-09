/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_WPDS_JOINLATTICETOSEMIRINGELEM_H_
#define PHASAR_PHASARLLVM_WPDS_JOINLATTICETOSEMIRINGELEM_H_

#include <iosfwd>
#include <memory>

#include "wali/SemElem.hpp"

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/JoinLattice.h"

namespace psr {
/**
 * An idempotent semi-ring is one whose addition is idempotent
 *
 *    a + a = a
 *
 * that is (R, +, 0) - is a join-semilattice with zero.
 */
template <typename V> class JoinLatticeToSemiRingElem : public wali::SemElem {
public:
  std::shared_ptr<EdgeFunction<V>> F;
  JoinLattice<V> &L;
  V ValV;

  JoinLatticeToSemiRingElem(std::shared_ptr<EdgeFunction<V>> F,
                            JoinLattice<V> &L)
      : wali::SemElem(), F(F), L(L) {}
   ~JoinLatticeToSemiRingElem() override = default;

  std::ostream &print(std::ostream &Os) const override { return Os << *F; }

  [[nodiscard]] wali::sem_elem_t one() const override {
    // std::cout << "JoinLatticeToSemiRingElem::one()" << std::endl;
    return wali::ref_ptr<JoinLatticeToSemiRingElem<V>>(
        new JoinLatticeToSemiRingElem(
            std::make_shared<AllBottom<V>>(L.bottomElement()), L));
  }

  [[nodiscard]] wali::sem_elem_t zero() const override {
    // std::cout << "JoinLatticeToSemiRingElem::zero()" << std::endl;
    return wali::ref_ptr<JoinLatticeToSemiRingElem<V>>(
        new JoinLatticeToSemiRingElem(
            std::make_shared<AllTop<V>>(L.topElement()), L));
  }

  wali::sem_elem_t extend(SemElem *Se) override {
    // std::cout << "JoinLatticeToSemiRingElem::extend()" << std::endl;
    auto ThisF = static_cast<JoinLatticeToSemiRingElem *>(this);
    auto ThatF = static_cast<JoinLatticeToSemiRingElem *>(Se);
    return wali::ref_ptr<JoinLatticeToSemiRingElem<V>>(
        new JoinLatticeToSemiRingElem(ThisF->F->composeWith(ThatF->F), L));
  }

  wali::sem_elem_t combine(SemElem *Se) override {
    // std::cout << "JoinLatticeToSemiRingElem::combine()" << std::endl;
    auto ThisF = static_cast<JoinLatticeToSemiRingElem *>(this);
    auto ThatF = static_cast<JoinLatticeToSemiRingElem *>(Se);
    return wali::ref_ptr<JoinLatticeToSemiRingElem<V>>(
        new JoinLatticeToSemiRingElem(ThisF->F->joinWith(ThatF->F), L));
  }

  bool equal(SemElem *Se) const override {
    // std::cout << "JoinLatticeToSemiRingElem::equal()" << std::endl;
    auto ThisF = static_cast<const JoinLatticeToSemiRingElem *>(this);
    auto ThatF = static_cast<const JoinLatticeToSemiRingElem *>(Se);
    return ThisF->F->equal_to(ThatF->F);
  }
};

template <typename V>
std::ostream &operator<<(std::ostream &Os,
                         const JoinLatticeToSemiRingElem<V> &ETS) {
  ETS.print(Os);
  return Os;
}

} // namespace psr

#endif
