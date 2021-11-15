/******************************************************************************
 * Copyright (c) 2020 Fabian Schiebel.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Fabian Schiebel and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_JOINEDGEFUNCTION_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_PROBLEMS_IDEGENERALIZEDLCA_JOINEDGEFUNCTION_H

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IDEGeneralizedLCA/IDEGeneralizedLCA.h"

namespace psr {

class JoinEdgeFunction : public EdgeFunction<IDEGeneralizedLCA::l_t>,
                         public std::enable_shared_from_this<JoinEdgeFunction> {
<<<<<<< HEAD
  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> Frst;
  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> Scnd;
=======
  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> First;
  std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> Second;
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  size_t MaxSize;

public:
  JoinEdgeFunction(
<<<<<<< HEAD
      const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &Frst,
      const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &Scnd,
      size_t MaxSize);
=======
      const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &First,
      const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &Second,
      size_t MaxSize);

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
<<<<<<< HEAD
  void print(std::ostream &OS, bool IsForDebug = false) const override;
=======

  void print(std::ostream &OS, bool IsForDebug = false) const override;

>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &getFirst() const;

  const std::shared_ptr<EdgeFunction<IDEGeneralizedLCA::l_t>> &
  getSecond() const;
};

} // namespace psr

#endif
