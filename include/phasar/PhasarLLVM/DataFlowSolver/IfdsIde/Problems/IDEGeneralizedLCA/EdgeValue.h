/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel, Alexander Meinhold and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_EDGEVALUE_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_EDGEVALUE_H

#include <iostream>
#include <memory>
#include <unordered_set>
#include <variant>

#include "llvm/ADT/APSInt.h"
#include "llvm/ADT/Twine.h"
#include "llvm/IR/Constant.h"
#include "llvm/IR/Instructions.h"

namespace psr {

enum class Ordering { Less, Greater, Equal, Incomparable };

class EdgeValue {
public:
  enum Type { Top, Integer, String, FloatingPoint };

private:
  std::variant<llvm::APInt, llvm::APFloat, std::string, std::nullptr_t> Value =
      nullptr;
  Type Typ;

public:
  EdgeValue(const llvm::Value *Val);
  EdgeValue(const EdgeValue &Ev);
  EdgeValue(llvm::APInt &&Vi);
  EdgeValue(const llvm::APInt &Vi);
  EdgeValue(llvm::APFloat &&Vf);
  EdgeValue(long long Vi);
  EdgeValue(int Vi);
  EdgeValue(double D);
  EdgeValue(float D);
  EdgeValue(std::string &&Vs);
  EdgeValue(std::nullptr_t);
  ~EdgeValue();
  const static EdgeValue To;
  bool tryGetInt(uint64_t &Res) const;
  bool tryGetFP(double &Res) const;
  bool tryGetString(std::string &Res) const;
  [[nodiscard]] bool isTop() const;
  [[nodiscard]] bool isNumeric() const;
  [[nodiscard]] bool isString() const;
  [[nodiscard]] enum Type getKind() const;
  // std::unique_ptr<ObjectLLVM> asObjLLVM(llvm::LLVMContext &ctx) const;
  [[nodiscard]] bool sqSubsetEq(const EdgeValue &Other) const;
  [[nodiscard]] EdgeValue performBinOp(llvm::BinaryOperator::BinaryOps Op,
                         const EdgeValue &Other) const;
  [[nodiscard]] EdgeValue typecast(enum Type Dest, unsigned Bits) const;
  EdgeValue &operator=(const EdgeValue &Ev);

  operator bool();
  friend bool operator==(const EdgeValue &V1, const EdgeValue &V2);

  // binary operators
  friend EdgeValue operator+(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator-(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator*(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator/(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator%(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator&(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator|(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator^(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator<<(const EdgeValue &V1, const EdgeValue &V2);
  friend EdgeValue operator>>(const EdgeValue &V1, const EdgeValue &V2);
  static int compare(const EdgeValue &V1, const EdgeValue &V2);

  // unary operators
  EdgeValue operator-() const;
  EdgeValue operator~() const;
  friend std::ostream &operator<<(std::ostream &Os, const EdgeValue &Ev);
  static std::string typeToString(enum Type Ty);
};
class EdgeValueSet;
using ev_t = EdgeValueSet;

ev_t performBinOp(llvm::BinaryOperator::BinaryOps Op, const ev_t &V1,
                  const ev_t &V2, size_t MaxSize);
ev_t performTypecast(const ev_t &Ev, EdgeValue::Type Dest, unsigned Bits);
Ordering compare(const ev_t &V1, const ev_t &V2);
ev_t join(const ev_t &V1, const ev_t &V2, size_t MaxSize);
/// \brief implements square subset equal
bool operator<(const ev_t &V1, const ev_t &V2);
bool isTopValue(const ev_t &V);
std::ostream &operator<<(std::ostream &Os, const ev_t &V);

} // namespace psr

namespace std {

template <> struct hash<psr::EdgeValue> {
  hash() = default;
  size_t operator()(const psr::EdgeValue &Val) const {
    auto Hc = hash<int>()(Val.getKind());
    uint64_t AsInt;
    double AsFloat;
    string AsString;
    if (Val.tryGetInt(AsInt)) {
      return hash<uint64_t>()(AsInt) * 31 + Hc;
    } 
    if (Val.tryGetFP(AsFloat)) {
      return hash<double>()(round(AsFloat)) * 31 + Hc;
    } 
    if (Val.tryGetString(AsString)) {
      return hash<string>()(AsString) * 31 + Hc;
    }
    return Hc;
  }
};

} // namespace std

#endif
