/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * CallString.h
 *
 *  Created on: 06.06.2017
 *      Author: philipp
 */

#ifndef PHASAR_PHASARLLVM_MONO_CALLSTRING_H_
#define PHASAR_PHASARLLVM_MONO_CALLSTRING_H_

#include <algorithm>
#include <deque>
#include <initializer_list>
#include <iosfwd>
#include <iterator>
#include <stdexcept>

namespace psr {

template <typename T, unsigned K> class CallString {
private:
  std::deque<T> Cs;
  static const unsigned ValK = K;

public:
  CallString() = default;
  CallString(std::initializer_list<T> Ilist) : Cs(Ilist) {
    if (Ilist.size() > ValK) {
      throw std::runtime_error(
          "initial call string length exceeds maximal length K");
    }
  }
  void push(T S) {
    if (Cs.size() > ValK - 1) {
      Cs.pop_front();
    }
    Cs.push_back(S);
  }
  T returnSite() {
    if (!Cs.empty()) {
      return Cs.back();
}
    return nullptr;
  }
  void pop() {
    if (!Cs.empty()) {
      Cs.pop_back();
    }
  }
  size_t size() { return Cs.size(); }
  std::deque<T> getInternalCS() const { return Cs; }
  friend bool operator==(const CallString &Lhs, const CallString &Rhs) {
    return Lhs.Cs == Rhs.Cs;
  }
  friend bool operator!=(const CallString &Lhs, const CallString &Rhs) {
    return !(Lhs == Rhs);
  }
  friend bool operator<(const CallString &Lhs, const CallString &Rhs) {
    return Lhs.Cs < Rhs.Cs;
  }
  friend std::ostream &operator<<(std::ostream &Os, const CallString &C) {
    std::copy(C.Cs.begin(), --C.Cs.end(), std::ostream_iterator<T>(Os, " * "));
    Os << C.Cs.back();
    return Os;
  }
};

} // namespace psr

#endif
