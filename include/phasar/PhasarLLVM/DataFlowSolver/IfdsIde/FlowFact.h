/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_FLOWFACT_H_
#define PHASAR_PHASARLLVM_IFDSIDE_FLOWFACT_H_

#include <iosfwd>
#include <type_traits>

namespace psr {

/// A common superclass of dataflow-facts used by non-template tabulation
/// problems (for example in plugins)
class FlowFact {
public:
  virtual ~FlowFact() = default;
  virtual void print(std::ostream &OS) const = 0;

  /// An abbreviation of an unsafe cast to T. Please use this only, if you know
  /// by 100% that this FlowFact is of type T
  template <typename T> T *as() {
    static_assert(std::is_base_of_v<FlowFact, T>);
    return reinterpret_cast<T *>(this);
  }
  /// An abbreviation of an unsafe cast to T. Please use this only, if you know
  /// by 100% that this FlowFact is of type T
  template <typename T> const T *as() const {
    static_assert(std::is_base_of_v<FlowFact, T>);
    return reinterpret_cast<const T *>(this);
  }
};

/*
static inline std::ostream &operator<<(std::ostream &OS, const FlowFact &F) {
  F.print(OS);
  return OS;
}
*/

} // namespace psr

#endif
