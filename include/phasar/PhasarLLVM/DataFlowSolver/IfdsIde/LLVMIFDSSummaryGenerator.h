/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * IFDSSummaryGenerator.h
 *
 *  Created on: 03.05.2017
 *      Author: philipp
 */

#ifndef PHASAR_PHASARLLVM_IFDSIDE_LLVMIFDSSUMMARYGENERATOR_H_
#define PHASAR_PHASARLLVM_IFDSIDE_LLVMIFDSSUMMARYGENERATOR_H_

#include <set>
#include <vector>

#include "llvm/IR/Function.h"
#include "llvm/IR/Instruction.h"
#include "llvm/IR/Value.h"

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/IFDSSummaryGenerator.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/LLVMIFDSSolver.h"
#include "phasar/Utils/LLVMShorthands.h"

namespace psr {

template <typename I, typename ConcreteIFDSTabulationProblem>
class LLVMIFDSSummaryGenerator
    : public IFDSSummaryGenerator<const llvm::Instruction *,
                                  const llvm::Value *, const llvm::Function *,
                                  I, ConcreteIFDSTabulationProblem,
                                  LLVMIFDSSolver<const llvm::Value *, I>> {
private:
  virtual std::vector<const llvm::Value *> getInputs() {
    std::vector<const llvm::Value *> Inputs;
    // collect arguments
    for (auto &Arg : this->toSummarize->args()) {
      Inputs.push_back(&Arg);
    }
    // collect global values
    auto Globals = globalValuesUsedinFunction(this->toSummarize);
    Inputs.insert(Inputs.end(), globals.begin(), globals.end());
    return Inputs;
  }

  virtual std::vector<bool>
  generateBitPattern(const std::vector<const llvm::Value *> &Inputs,
                     const std::set<const llvm::Value *> &Subset) {
    // initialize all bits to zero
    std::vector<bool> Bitpattern(Inputs.size(), false);
    if (Subset.empty()) {
      return Bitpattern;
    }
    for (const auto *Elem : Subset) {
      for (size_t VarI = 0; VarI < Inputs.size(); ++VarI) {
        if (Elem == Inputs[VarI]) {
          Bitpattern[VarI] = true;
        }
      }
    }
    return Bitpattern;
  }

public:
  LLVMIFDSSummaryGenerator(const llvm::Function *F, I Icfg,
                           SummaryGenerationStrategy S)
      : IFDSSummaryGenerator<const llvm::Instruction *, const llvm::Value *,
                             const llvm::Function *, I,
                             ConcreteIFDSTabulationProblem,
                             LLVMIFDSSolver<const llvm::Value *, I>>(F, icfg,
                                                                     S) {}

  virtual ~LLVMIFDSSummaryGenerator() = default;
};
} // namespace psr

#endif
