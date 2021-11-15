/**
 * @author Sebastian Roland <seroland86@gmail.com>
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_FLOWFUNCTIONS_PHINODEFLOWFUNCTION_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_FLOWFUNCTIONS_PHINODEFLOWFUNCTION_H

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSFieldSensTaintAnalysis/FlowFunctions/FlowFunctionBase.h"

namespace psr {

class PHINodeFlowFunction : public FlowFunctionBase {
public:
  PHINodeFlowFunction(const llvm::Instruction *CurrentInst,
                      class TraceStats &TraceStats, const ExtendedValue& ZeroValue)
      : FlowFunctionBase(CurrentInst, TraceStats, ZeroValue) {}
  ~PHINodeFlowFunction() override = default;

  std::set<ExtendedValue> computeTargetsExt(ExtendedValue &Fact) override;
};

} // namespace psr

#endif // PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_FLOWFUNCTIONS_PHINODEFLOWFUNCTION_H
