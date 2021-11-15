/**
 * @author Sebastian Roland <seroland86@gmail.com>
 */

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSFieldSensTaintAnalysis/FlowFunctions/GenerateFlowFunction.h"

namespace psr {

std::set<ExtendedValue>
GenerateFlowFunction::computeTargetsExt(ExtendedValue &Fact) {
  TraceStats.add(CurrentInst);

  if (Fact == ZeroValue) {
    return {ExtendedValue(CurrentInst)};
  }

  return {Fact};
}

} // namespace psr
