/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_BINARYEDGEFUNCTION_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_BINARYEDGEFUNCTION_H

#include <utility>

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/EdgeValueSet.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/IDEGeneralizedLCA.h"

namespace psr {

class BinaryEdgeFunction
    : public EdgeFunction<IDEGeneralizedLCA::l_t>,
      public std::enable_shared_from_this<BinaryEdgeFunction> {
  llvm::BinaryOperator::BinaryOps Op;
<<<<<<< HEAD
  const IDEGeneralizedLCA::l_t Cnst;
=======
  const IDEGeneralizedLCA::l_t Const;
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  bool LeftConst;
  size_t MaxSize;

public:
  BinaryEdgeFunction(llvm::BinaryOperator::BinaryOps Op,
<<<<<<< HEAD
                     IDEGeneralizedLCA::l_t Cnst, bool LeftConst,
                     size_t MaxSize)
      : Op(Op), Cnst(std::move(Cnst)), LeftConst(LeftConst), MaxSize(MaxSize) {}
=======
                     const IDEGeneralizedLCA::l_t &Const, bool LeftConst,
                     size_t MaxSize)
      : Op(Op), Const(Const), LeftConst(LeftConst), MaxSize(MaxSize) {}
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506

  IDEGeneralizedLCA::l_t computeTarget(IDEGeneralizedLCA::l_t Source) override;

  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> composeWith(
      std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> SecondFunction)
      override;

  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>>
  joinWith(std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> OtherFunction)
      override;

  bool equalTo(std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> Other)
      const override;

  void print(std::ostream &OS, bool IsForDebug = false) const override;
};

} // namespace psr

#endif
