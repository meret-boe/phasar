#ifndef _PHASAR_PHASARLLVM_MONO_CALLSTRINGCTX_H_
#define _PHASAR_PHASARLLVM_MONO_CALLSTRINGCTX_H_

#include <deque>
#include <functional>
#include <initializer_list>

#include "boost/functional/hash.hpp"
#include "phasar/Utils/LLVMShorthands.h"

namespace psr {

template <typename N, unsigned K> class CallStringCTX {
protected:
  std::deque<N> Cs;
  static const unsigned ValK = K;
  friend struct std::hash<psr::CallStringCTX<N, K>>;

public:
  CallStringCTX() = default;

  CallStringCTX(std::initializer_list<N> Ilist) : Cs(Ilist) {
    if (Ilist.size() > ValK) {
      throw std::runtime_error(
          "initial call std::string length exceeds maximal length K");
    }
  }

  void pushBack(N ValN) {
    if (Cs.size() > ValN - 1) {
      Cs.pop_front();
    }
    Cs.push_back(ValN);
  }

  N popBack() {
    if (!Cs.empty()) {
      N ValN = Cs.back();
      Cs.pop_back();
      return ValN;
    }
    return N{};
  }

  bool isEqual(const CallStringCTX &Rhs) const { return Cs == Rhs.Cs; }

  bool isDifferent(const CallStringCTX &Rhs) const { return !isEqual(Rhs); }

  friend bool operator==(const CallStringCTX<N, K> &Lhs,
                         const CallStringCTX<N, K> &Rhs) {
    return Lhs.isEqual(Rhs);
  }

  friend bool operator!=(const CallStringCTX<N, K> &Lhs,
                         const CallStringCTX<N, K> &Rhs) {
    return !Lhs.isEqual(Rhs);
  }

  friend bool operator<(const CallStringCTX<N, K> &Lhs,
                        const CallStringCTX<N, K> &Rhs) {
    return Lhs.Cs < Rhs.Cs;
  }

  void print(std::ostream &Os) const {
    Os << "Call string: [ ";
    for (auto C : Cs) {
      Os << llvmIRToString(C);
      if (C != Cs.back()) {
        Os << " * ";
      }
    }
    Os << " ]";
  }

  friend std::ostream &operator<<(std::ostream &Os,
                                  const CallStringCTX<N, K> &C) {
    C.print(Os);
    return Os;
  }

  [[nodiscard]] bool empty() const { return Cs.empty(); }

  [[nodiscard]] std::size_t size() const { return Cs.size(); }
};

} // namespace psr

namespace std {

template <typename N, unsigned K> struct hash<psr::CallStringCTX<N, K>> {
  size_t operator()(const psr::CallStringCTX<N, K> &CS) const noexcept {
    boost::hash<std::deque<N>> HashDeque;
    std::hash<unsigned> HashUnsigned;
    size_t U = HashUnsigned(K);
    size_t H = HashDeque(CS.cs);
    return U ^ (H << 1);
  }
};

} // namespace std

#endif