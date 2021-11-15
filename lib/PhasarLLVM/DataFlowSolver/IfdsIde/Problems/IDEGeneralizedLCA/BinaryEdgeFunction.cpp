/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel and others
 *****************************************************************************/

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/BinaryEdgeFunction.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/LCAEdgeFunctionComposer.h"

namespace psr {

IDEGeneralizedLCA::l_t
BinaryEdgeFunction::computeTarget(IDEGeneralizedLCA::l_t Source) {
  /*auto ret = leftConst ? performBinOp(op, cnst, source, maxSize)
                       : performBinOp(op, source, cnst, maxSize);
  std::cout << "Binary(" << source << ") = " << ret << std::endl;
  return ret;*/
  if (LeftConst) {
<<<<<<< HEAD
    return performBinOp(Op, Cnst, Source, MaxSize);
  }  
    return performBinOp(Op, Source, Cnst, MaxSize);
 
=======
    return performBinOp(Op, Const, Source, MaxSize);
  }
  return performBinOp(Op, Source, Const, MaxSize);
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
}

std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>>
BinaryEdgeFunction::composeWith(
    std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> SecondFunction) {
  if (auto *EI = dynamic_cast<EdgeIdentity<IDEGeneralizedLCA::l_t> *>(
          SecondFunction.get())) {
    return this->shared_from_this();
  }
  if (dynamic_cast<AllBottom<IDEGeneralizedLCA::l_t> *>(SecondFunction.get())) {
    // print(std::cout << "Compose ");
    // std::cout << " with ALLBOT" << std::endl;
    return shared_from_this();
  }
  return std::make_shared<LCAEdgeFunctionComposer>(this->shared_from_this(),
                                                   SecondFunction, MaxSize);
}

std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>>
BinaryEdgeFunction::joinWith(
    std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> OtherFunction) {
  if (OtherFunction.get() == this ||
      OtherFunction->equalTo(this->shared_from_this())) {
    return this->shared_from_this();
  }
  if (auto *AT =
          dynamic_cast<AllTop<IDEGeneralizedLCA::l_t> *>(OtherFunction.get())) {
    return this->shared_from_this();
  }
  return std::make_shared<AllBottom<IDEGeneralizedLCA::l_t>>(
      IDEGeneralizedLCA::l_t({EdgeValue(nullptr)}));
}

bool BinaryEdgeFunction::equalTo(
    std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> Other) const {
  return this == Other.get();
}

<<<<<<< HEAD
void BinaryEdgeFunction::print(std::ostream &OS, bool  /*IsForDebug*/) const {
=======
void BinaryEdgeFunction::print(std::ostream &OS,
                               [[maybe_unused]] bool IsForDebug) const {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  OS << "Binary_" << Op;
}

} // namespace psr
