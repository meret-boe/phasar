/**
 * @author Sebastian Roland <seroland86@gmail.com>
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_STATS_LINENUMBERENTRY_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_STATS_LINENUMBERENTRY_H

#include <functional>

namespace psr {

class LineNumberEntry {
public:
  LineNumberEntry(unsigned int LineNumber) : LineNumber(LineNumber) {}
  ~LineNumberEntry() = default;

  bool operator<(const LineNumberEntry &Rhs) const {
    return std::less<unsigned int>{}(LineNumber, Rhs.LineNumber);
  }

  [[nodiscard]] unsigned int getLineNumber() const { return LineNumber; }

  [[nodiscard]] bool isReturnValue() const { return ReturnValue; }
  void setReturnValue(bool RetVal) { ReturnValue = RetVal; }

private:
  unsigned int LineNumber;
  bool ReturnValue = false;
};

} // namespace psr

#endif // PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_IFDSFIELDSENSTAINTANALYSIS_STATS_LINENUMBERENTRY_H
