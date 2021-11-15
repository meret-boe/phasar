/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel, Alexander Meinhold and others
 *****************************************************************************/

#include <cassert>

#include "llvm/ADT/APFloat.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/GlobalVariable.h"
#include "llvm/Support/raw_ostream.h"

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/EdgeValue.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/EdgeValueSet.h"

namespace psr {
using std::get;

std::ostream &printSemantics(const llvm::APFloat &Fl) {
  if (&Fl.getSemantics() == &llvm::APFloat::IEEEdouble()) {
    return std::cout << "IEEEdouble";
  } if (&Fl.getSemantics() == &llvm::APFloat::IEEEhalf()) {
    return std::cout << "IEEEhalf";
  } if (&Fl.getSemantics() == &llvm::APFloat::IEEEquad()) {
    return std::cout << "IEEEquad";
  } if (&Fl.getSemantics() == &llvm::APFloat::IEEEsingle()) {
    return std::cout << "IEEEsingle";
  } if (&Fl.getSemantics() == &llvm::APFloat::PPCDoubleDouble()) {
    return std::cout << "PPCDoubleDouble";
  } if (&Fl.getSemantics() == &llvm::APFloat::x87DoubleExtended()) {
    return std::cout << "x87DoubleExtended";
  } if (&Fl.getSemantics() == &llvm::APFloat::Bogus()) {
    return std::cout << "Bogus";
  }
  return std::cout << "Sth else";
 
}

const EdgeValue EdgeValue::To = EdgeValue(nullptr);

EdgeValue::EdgeValue(const llvm::Value *Val) : Typ(Top) {
  if (const auto *Cnst = llvm::dyn_cast<llvm::Constant>(Val)) {
    if (Cnst->getType()->isIntegerTy()) {
      Typ = Integer;
      Value = llvm::APInt(llvm::cast<llvm::ConstantInt>(Cnst)->getValue());
    } else if (Cnst->getType()->isFloatingPointTy()) {
      Typ = FloatingPoint;
      const auto &CnstFP = llvm::cast<llvm::ConstantFP>(Cnst)->getValueAPF();

      llvm::APFloat Apf(CnstFP);
      bool Unused;
      Apf.convert(llvm::APFloat::IEEEdouble(),
                  llvm::APFloat::roundingMode::NearestTiesToEven, &Unused);
      Value = llvm::APFloat(Apf);
    } else if (llvm::isa<llvm::ConstantPointerNull>(Cnst)) {
      Typ = String;
      Value = std::string();
    } else if (Cnst->getType()->isPointerTy() &&
               Cnst->getType()->getPointerElementType()->isIntegerTy()) {
      Typ = String;
      const auto *Gep = llvm::cast<llvm::ConstantExpr>(
          Cnst); // already checked, hence cast instead of dyn_cast
      if (auto *Glob =
              llvm::dyn_cast<llvm::GlobalVariable>(Gep->getOperand(0))) {
        Value = std::string(
            llvm::cast<llvm::ConstantDataArray>(Glob->getInitializer())
                ->getAsCString()
                .str());
      } else {
        // inttoptr
        Value = nullptr;
        Typ = Top;
      }
    } else {
      Value = nullptr;
      Typ = Top;
    }
  } else {
    Value = nullptr;
    Typ = Top;
  }
}

EdgeValue::EdgeValue(const EdgeValue &Ev) : Typ(Ev.Typ) {
  switch (Typ) {
  case Top:
    Value = nullptr;
    break;
  case Integer:
    Value = std::get<llvm::APInt>(Ev.Value);
    break;
  case FloatingPoint:
    Value = std::get<llvm::APFloat>(Ev.Value);
    break;
  case String:
    Value = std::get<std::string>(Ev.Value);
    break;
  }
}

EdgeValue &EdgeValue::operator=(const EdgeValue &Ev) {
  this->~EdgeValue();
  new (this) EdgeValue(Ev);
  return *this;
}

EdgeValue::~EdgeValue() { Value.~variant(); }

EdgeValue::EdgeValue(llvm::APInt &&Vi) : Typ(EdgeValue::Integer) {
  Value = llvm::APInt(std::move(Vi));
}

EdgeValue::EdgeValue(const llvm::APInt &Vi) : Typ(EdgeValue::Integer) {
  Value = llvm::APInt(Vi);
}

EdgeValue::EdgeValue(llvm::APFloat &&Vf) : Typ(EdgeValue::FloatingPoint) {
  llvm::APFloat Fp = llvm::APFloat(std::move(Vf));
  bool Unused;
  Fp.convert(llvm::APFloat::IEEEdouble(),
             llvm::APFloat::roundingMode::NearestTiesToEven, &Unused);
  Value = Fp;
}

EdgeValue::EdgeValue(long long Vi) : Typ(EdgeValue::Integer) {
  Value = llvm::APInt(llvm::APInt(sizeof(long long) << 3, Vi));
}

EdgeValue::EdgeValue(int Vi) : Typ(EdgeValue::Integer) {
  Value = llvm::APInt(llvm::APInt(sizeof(int) << 3, Vi));
}

EdgeValue::EdgeValue(double D) : Typ(EdgeValue::FloatingPoint) {
  Value = llvm::APFloat(D);
}

EdgeValue::EdgeValue(float D) : Typ(EdgeValue::FloatingPoint) {
  Value = llvm::APFloat(D);
}

EdgeValue::EdgeValue(std::string &&Vs) : Typ(EdgeValue::String) {
  Value = std::string(Vs);
}

EdgeValue::EdgeValue(std::nullptr_t) : Typ(EdgeValue::Top) {}
bool EdgeValue::tryGetInt(uint64_t &Res) const {
  if (Typ != Integer) {
    return false;
}
  Res = std::get<llvm::APInt>(Value).getLimitedValue();
  return true;
}

bool EdgeValue::tryGetFP(double &Res) const {
  if (Typ != FloatingPoint) {
    return false;
}
  Res = std::get<llvm::APFloat>(Value).convertToDouble();
  return true;
}

bool EdgeValue::tryGetString(std::string &Res) const {
  if (Typ != String) {
    return false;
}
  Res = get<std::string>(Value);
  return true;
}

bool EdgeValue::isTop() const { return Typ == Top; }

bool EdgeValue::isNumeric() const {
  return Typ == Integer || Typ == FloatingPoint;
}

bool EdgeValue::isString() const { return Typ == String; }

EdgeValue::Type EdgeValue::getKind() const { return Typ; }

EdgeValue::operator bool() {
  switch (Typ) {
  case Integer:
    return !std::get<llvm::APInt>(Value).isNullValue();
  case FloatingPoint:
    return std::get<llvm::APFloat>(Value).isNonZero();
  case String:
    return !std::get<std::string>(Value).empty();
  default:
    break;
  }
  return false;
}

bool operator==(const EdgeValue &V1, const EdgeValue &V2) {
  // std::cout << "Compare edge values" << std::endl;
  if (V1.Typ != V2.Typ) {
    // std::cout << "Comparing incompatible types" << std::endl;
    return false;
  }
  switch (V1.Typ) {
  case EdgeValue::Top:
    return true;
  case EdgeValue::Integer:
    // if (v1.Value.asInt != v2.Value.asInt)
    //  std::cout << "integer unequal" << std::endl;
    return std::get<llvm::APInt>(V1.Value) == std::get<llvm::APInt>(V2.Value);
  case EdgeValue::FloatingPoint: {
    // std::cout << "compare floating points" << std::endl;
    auto Cp = std::get<llvm::APFloat>(V1.Value).compare(
        std::get<llvm::APFloat>(V2.Value));
    if (Cp == llvm::APFloat::cmpResult::cmpEqual) {
      // std::cout << "FP equal" << std::endl;
      return true;
    }
    auto D1 = std::get<llvm::APFloat>(V1.Value).convertToDouble();
    auto D2 = std::get<llvm::APFloat>(V2.Value).convertToDouble();

    const double Epsilon = 0.000001;
    // std::cout << "Compare " << d1 << " against " << d2 << std::endl;
    return D1 == D2 || D1 - D2 < Epsilon || D2 - D1 < Epsilon;
  }
  case EdgeValue::String:
    // if (v1.Value.asString != v2.Value.asString)
    //  std::cout << "String unequal" << std::endl;
    return std::get<std::string>(V1.Value) == std::get<std::string>(V2.Value);
  default: // will not happen
    std::cerr << "FATAL ERROR" << std::endl;
    return false;
  }
}

bool EdgeValue::sqSubsetEq(const EdgeValue &Other) const {
  return Other.isTop() || Other.Typ == Typ;
}

// binary operators
EdgeValue operator+(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) +
                     std::get<llvm::APInt>(V2.Value));
  case EdgeValue::FloatingPoint:
    return EdgeValue(std::get<llvm::APFloat>(V1.Value) +
                     std::get<llvm::APFloat>(V2.Value));
  case EdgeValue::String:
    return EdgeValue(std::get<std::string>(V1.Value) +
                     std::get<std::string>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator-(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) -
                     std::get<llvm::APInt>(V1.Value));
  case EdgeValue::FloatingPoint:
    // printSemantics(v1.Value.asFP) << " <=> ";
    // printSemantics(v2.Value.asFP) << std::endl;
    return EdgeValue(std::get<llvm::APFloat>(V1.Value) -
                     std::get<llvm::APFloat>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator*(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) *
                     std::get<llvm::APInt>(V2.Value));
  case EdgeValue::FloatingPoint:
    return EdgeValue(std::get<llvm::APFloat>(V1.Value) *
                     std::get<llvm::APFloat>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator/(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(
        std::get<llvm::APInt>(V1.Value).sdiv(std::get<llvm::APInt>(V2.Value)));
  case EdgeValue::FloatingPoint:
    return EdgeValue(std::get<llvm::APFloat>(V1.Value) /
                     std::get<llvm::APFloat>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator%(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(
        std::get<llvm::APInt>(V1.Value).srem(std::get<llvm::APInt>(V2.Value)));
  case EdgeValue::FloatingPoint: {
    llvm::APFloat Fl = std::get<llvm::APFloat>(V1.Value);
    Fl.remainder(std::get<llvm::APFloat>(V2.Value));
    return EdgeValue(std::move(Fl));
  }
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator&(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) &
                     std::get<llvm::APInt>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator|(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) |
                     std::get<llvm::APInt>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator^(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value) ^
                     std::get<llvm::APInt>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator<<(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(std::get<llvm::APInt>(V1.Value)
                     << std::get<llvm::APInt>(V2.Value));
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue operator>>(const EdgeValue &V1, const EdgeValue &V2) {
  if (V1.Typ != V2.Typ) {
    return EdgeValue(nullptr);
}
  switch (V1.Typ) {
  case EdgeValue::Integer:
    return EdgeValue(
        std::get<llvm::APInt>(V1.Value).ashr(std::get<llvm::APInt>(V2.Value)));
  default:
    return EdgeValue(nullptr);
  }
}

// unary operators
EdgeValue EdgeValue::operator-() const {
  if (Typ == Integer) {
    return EdgeValue(-std::get<llvm::APInt>(Value));
}
  return EdgeValue(nullptr);
}

EdgeValue EdgeValue::operator~() const {
  if (Typ == Integer) {
    return EdgeValue(~std::get<llvm::APInt>(Value));
}
  return EdgeValue(nullptr);
}

int EdgeValue::compare(const EdgeValue &V1, const EdgeValue &V2) { //NOLINT
  switch (V1.Typ) {
  case EdgeValue::Integer: {
    auto V1val = std::get<llvm::APInt>(V1.Value).getLimitedValue();
    uint64_t V2val;
    double V2valFp;
    if (V2.tryGetInt(V2val)) {
      return V1val - V2val; //NOLINT
    } if (V2.tryGetFP(V2valFp)) {
      return V1val < V2valFp ? -1 : (V1val > V2valFp ? 1 : 0); //NOLINT
    }
    break;
  }
  case EdgeValue::FloatingPoint: {
    auto V1val = std::get<llvm::APFloat>(V1.Value).convertToDouble();
    uint64_t V2val;
    double V2valFp;
    bool IsInt = V2.tryGetInt(V2val);
    if (IsInt || V2.tryGetFP(V2valFp)) {
      if (IsInt) {
        V2valFp = V2val; //NOLINT
}

      return V1val < V2valFp ? -1 : (V1val > V2valFp ? 1 : 0);
    }

    break;
  }
  case EdgeValue::String: {
    std::string V2val;
    if (V2.tryGetString(V2val)) {
      return std::get<std::string>(V1.Value).compare(V2val);
    }
    break;
  }
  default:
    break;
  }

  return 0;
}

std::ostream &operator<<(std::ostream &Os, const EdgeValue &Ev) {
  switch (Ev.Typ) {
  case EdgeValue::Integer: {
    std::string S;
    llvm::raw_string_ostream Ros(S);
    Ros << std::get<llvm::APInt>(Ev.Value);
    return Os << Ros.str();
  }
  case EdgeValue::String:
    return Os << "\"" << std::get<std::string>(Ev.Value) << "\"";
  case EdgeValue::FloatingPoint: {
    return Os << std::get<llvm::APFloat>(Ev.Value).convertToDouble();
  }
  default:
    return Os << "<TOP>";
  }
}

EdgeValue EdgeValue::typecast(enum Type Dest, unsigned Bits) const {
  switch (Dest) {

  case Integer:
    switch (Typ) {
    case Integer:
      if (std::get<llvm::APInt>(Value).getBitWidth() <= Bits) {
        return *this;
      } else {
        return EdgeValue(std::get<llvm::APInt>(Value) & ((1 << Bits) - 1));
}
    case FloatingPoint: {
      bool Unused;
      llvm::APSInt Ai;
      std::get<llvm::APFloat>(Value).convertToInteger(
          Ai, llvm::APFloat::roundingMode::NearestTiesToEven, &Unused);
      return EdgeValue(Ai);
    }
    default:
      return EdgeValue(nullptr);
    }
  case FloatingPoint:
    switch (Typ) {
    case Integer:
      if (Bits > 32) {
        return EdgeValue((double)std::get<llvm::APInt>(Value).getSExtValue());
      } else {
        return EdgeValue((float)std::get<llvm::APInt>(Value).getSExtValue());
}
    case FloatingPoint:
      return *this;
    default:
      return EdgeValue(nullptr);
    }
  default:
    return EdgeValue(nullptr);
  }
}

EdgeValue EdgeValue::performBinOp(llvm::BinaryOperator::BinaryOps Op,
                                  const EdgeValue &Other) const {
  switch (Op) {
  case llvm::BinaryOperator::BinaryOps::Add:
  case llvm::BinaryOperator::BinaryOps::FAdd:
    return *this + Other;
  case llvm::BinaryOperator::BinaryOps::And:
    return *this & Other;
  case llvm::BinaryOperator::BinaryOps::AShr:
    return *this >> Other;
  case llvm::BinaryOperator::BinaryOps::FDiv:
  case llvm::BinaryOperator::BinaryOps::SDiv:
    return *this / Other;
  case llvm::BinaryOperator::BinaryOps::LShr: {
    if (Typ != Other.Typ) {
      return EdgeValue(nullptr);
}
    switch (Typ) {
    case EdgeValue::Integer:
      return EdgeValue(std::get<llvm::APInt>(Value).lshr(
          std::get<llvm::APInt>(Other.Value)));
    default:
      return EdgeValue(nullptr);
    }
  }
  case llvm::BinaryOperator::BinaryOps::Mul:
  case llvm::BinaryOperator::BinaryOps::FMul:
    return *this * Other;
  case llvm::BinaryOperator::BinaryOps::Or:
    return *this | Other;
  case llvm::BinaryOperator::BinaryOps::Shl:
    return *this << Other;
  case llvm::BinaryOperator::BinaryOps::SRem:
  case llvm::BinaryOperator::BinaryOps::FRem:
    return *this % Other;
  case llvm::BinaryOperator::BinaryOps::Sub:
  case llvm::BinaryOperator::BinaryOps::FSub:
    return *this - Other;
  case llvm::BinaryOperator::BinaryOps::UDiv: {
    if (Typ != Other.Typ) {
      return EdgeValue(nullptr);
}
    switch (Typ) {
    case EdgeValue::Integer:
      return EdgeValue(std::get<llvm::APInt>(Value).udiv(
          std::get<llvm::APInt>(Other.Value)));
    default:
      return EdgeValue(nullptr);
    }
  }
  case llvm::BinaryOperator::BinaryOps::URem: {
    if (Typ != Other.Typ) {
      return EdgeValue(nullptr);
}
    switch (Typ) {
    case EdgeValue::Integer:
      return EdgeValue(std::get<llvm::APInt>(Value).urem(
          std::get<llvm::APInt>(Other.Value)));
    default:
      return EdgeValue(nullptr);
    }
  }
  case llvm::BinaryOperator::BinaryOps::Xor:
    return *this ^ Other;
  default:
    return EdgeValue(nullptr);
  }
}

ev_t performBinOp(llvm::BinaryOperator::BinaryOps Op, const ev_t &V1,
                  const ev_t &V2, size_t MaxSize) {
  // std::cout << "Perform Binop on " << v1 << " and " << v2 << std::endl;

  if (V1.empty() || isTopValue(V1) || V2.empty() || isTopValue(V2)) {
    // std::cout << "\t=> <TOP>" << std::endl;
    return {EdgeValue(nullptr)};
  }
  ev_t Ret({});
  for (const auto &Ev1 : V1) {
    for (const auto &Ev2 : V2) {

      Ret.insert(Ev1.performBinOp(Op, Ev2));
      if (Ret.size() > MaxSize) {
        // std::cout << "\t=> <TOP>" << std::endl;
        return ev_t({EdgeValue(nullptr)});
      }
    }
  }
  // std::cout << "\t=> " << ret << std::endl;
  return Ret;
}

ev_t performTypecast(const ev_t &Ev, EdgeValue::Type Dest, unsigned Bits) {
  if (Ev.empty() || isTopValue(Ev)) {
    // std::cout << "\t=> <TOP>" << std::endl;
    return {EdgeValue(nullptr)};
  }
  ev_t Ret({});
  for (const auto &V : Ev) {
    auto Tc = V.typecast(Dest, Bits);
    if (Tc.isTop()) {
      return ev_t({EdgeValue(nullptr)});
}
    Ret.insert(Tc);
  }
  return Ret;
}

Ordering compare(const ev_t &V1, const ev_t &V2) {
  const auto &Smaller = V1.size() <= V2.size() ? V1 : V2;
  const auto &Larger = V1.size() > V2.size() ? V1 : V2;

  for (const auto &Elem : Smaller) {
    if (!Larger.count(Elem)) {
      return Ordering::Incomparable;
    }
  }
  return V1.size() == V2.size()
             ? Ordering::Equal
             : (&Smaller == &V1 ? Ordering::Less : Ordering::Greater);
}

ev_t join(const ev_t &V1, const ev_t &V2, size_t MaxSize) {
  // std::cout << "Join " << v1 << " and " << v2 << std::endl;
  if (isTopValue(V1) || isTopValue(V2)) {
    // std::cout << "\t=> <TOP>" << std::endl;
    return {EdgeValue(nullptr)};
  }
  ev_t Ret(V1.begin(), V1.end());

  for (const auto &Elem : V2) {
    Ret.insert(Elem);
    if (Ret.size() > MaxSize) {
      // std::cout << "\t=> <TOP>" << std::endl;
      return {EdgeValue(nullptr)};
    }
  }
  // std::cout << "\t=> " << ret << std::endl;

  return Ret;
}

bool isTopValue(const ev_t &V) { return V.size() == 1 && V.begin()->isTop(); }
std::ostream &operator<<(std::ostream &Os, const ev_t &V) {
  Os << "{";
  bool Frst = true;
  for (const auto &Elem : V) {
    if (Frst) {
      Frst = false;
    } else {
      Os << ", ";
}
    Os << Elem;
  }
  return Os << "}";
}

bool operator<(const ev_t &V1, const ev_t &V2) {
  if (V1.size() >= V2.size()) {
    return V1 != V2 && (V1.empty() || V2 == ev_t({EdgeValue::Top}));
  }  
    for (const auto &Elem : V1) { //NOLINT
      if (!V2.count(Elem)) {
        return false;
}
    }
    return true;
 
}

std::string EdgeValue::typeToString(enum Type Ty) {
  switch (Ty) {
  case Integer:
    return "Integer";
  case FloatingPoint:
    return "FloatingPoint";
  case String:
    return "String";
  default:
    return "Top";
  }
}

} // namespace psr
