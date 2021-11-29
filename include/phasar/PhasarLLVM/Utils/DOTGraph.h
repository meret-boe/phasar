/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * DOTGraph.h
 *
 *  Created on: 31.08.2019
 *      Author: rleer
 */

#ifndef PHASAR_PHASARLLVM_UTILS_DOTGRAPH_H_
#define PHASAR_PHASARLLVM_UTILS_DOTGRAPH_H_

#include <iosfwd>
#include <map>
#include <set>
#include <string>

#include "phasar/Config/Configuration.h"
#include "phasar/Utils/Utilities.h"

namespace psr {

class DOTConfig {
public:
  static  std::string cfNodeAttr() { return CFNode; }
  static  std::string cfIntraEdgeAttr() { return CFIntraEdge; }
  static  std::string cfInterEdgeAttr() { return CFInterEdge; }

  static  std::string factNodeAttr() { return FactNode; }
  static  std::string factIdEdgeAttr() { return FactIDEdge; }
  static  std::string factCrossEdgeAttr() { return FactCrossEdge; }
  static  std::string factInterEdgeAttr() { return FactInterEdge; }

  static  std::string lambdaNodeAttr() { return LambdaNode; }
  static  std::string lambdaIdEdgeAttr() { return LambdaIDEdge; }
  static  std::string lambdaInterEdgeAttr() { return LambdaInterEdge; }

  static void
  importDOTConfig(std::string ConfigPath = PhasarConfig::phasarDirectory());

  static DOTConfig &getDOTConfig();
  ~DOTConfig() = default;
  DOTConfig(const DOTConfig &) = delete;
  DOTConfig(DOTConfig &&) = delete;
  DOTConfig& operator=(DOTConfig Other);
  DOTConfig& operator=(DOTConfig&& Other) noexcept;

private:
  DOTConfig() = default;
  inline static const std::string FontSize = "fontsize=11";
  inline static const std::string ArrowSize = "arrowsize=0.7";

  inline static std::string CFNode = "node [style=filled, shape=record]"; //NOLINT
  inline static std::string CFIntraEdge = "Edge []";                      //NOLINT
  inline static std::string CFInterEdge = "Edge [weight=0.1]";            //NOLINT
  inline static std::string FactNode = "node [style=rounded]";            //NOLINT
  inline static std::string FactIDEdge =                                  //NOLINT
      "Edge [style=dotted, arrowhead=normal, " + FontSize + ", " + ArrowSize +
      ']';
  inline static std::string FactCrossEdge =                               //NOLINT
      "Edge [style=dotted, arrowhead=normal, " + FontSize + ", " + ArrowSize +
      ']';
  inline static std::string FactInterEdge =                               //NOLINT
      "Edge [weight=0.1, style=dashed, " + FontSize + ", " + ArrowSize + ']';
  inline static std::string LambdaNode = "node [style=rounded]";          //NOLINT
  inline static std::string LambdaIDEdge =                                //NOLINT
      "Edge [style=dotted, arrowhead=normal, " + FontSize + ", " + ArrowSize +
      ']';
  inline static std::string LambdaInterEdge =                             //NOLINT
      "Edge [weight=0.1, style=dashed, " + FontSize + ", " + ArrowSize + ']';
};

struct DOTNode {
  /* stmt Id = <func-name>_<stmt-Id>
   * e.g. foo_42 for the statement w/ Id 42 in function foo
   *
   * fact Id = <func-name>_<fact-Id>_<stmt-Id>
   * e.g. fact i in main valid at statement 3 : main_1_3
   * Note that in this example fact i is the very first valid fact
   * encountered so the fact-Id is 1 (zero value has fact-Id 0)
   */
  std::string Id;
  std::string FuncName;
  std::string Label;
  std::string StmtId;
  unsigned FactId{};
  bool IsVisible = true;

  DOTNode() = default;
  DOTNode(std::string FName, std::string L, std::string SId, unsigned FId = 0,
          bool IsStmt = true, bool Isv = true);
  [[nodiscard]] std::string str(const std::string &Indent = "") const;
};

bool operator<(const DOTNode &Lhs, const DOTNode &Rhs);
bool operator==(const DOTNode &Lhs, const DOTNode &Rhs);
std::ostream &operator<<(std::ostream &Os, const DOTNode &Node);

struct DOTEdge {
  DOTNode Source;
  DOTNode Target;
  bool IsVisible;
  std::string EdgeFnLabel;
  std::string ValueLabel;

  DOTEdge(DOTNode Src, DOTNode Tar, bool Isv = true, std::string Efl = "",
          std::string Vl = "");
  [[nodiscard]] std::string str(const std::string &Indent = "") const;
};

bool operator<(const DOTEdge &Lhs, const DOTEdge &Rhs);
std::ostream &operator<<(std::ostream &Os, const DOTEdge &Edge);

struct DOTFactSubGraph {
  // fact subgraph Id = <func-name>_<fact-Id>
  std::string Id;
  unsigned FactId;
  std::string Label;
  // stmt-Id -> fact-node
  std::map<std::string, DOTNode, StringIdLess> Nodes;
  std::set<DOTEdge> Edges;

  [[nodiscard]] std::string str(const std::string &Indent = "") const;
};

std::ostream &operator<<(std::ostream &Os, const DOTFactSubGraph &FactSg);

struct DOTFunctionSubGraph {
  // function subgraph Id = <func-name>
  std::string Id;
  std::set<DOTNode> Stmts;
  // fact-Id -> fact-subgraph
  std::map<unsigned, DOTFactSubGraph> Facts;
  std::set<DOTEdge> IntraCfEdges;
  /// d1 -> d2 where d1 != d2
  std::set<DOTEdge> CrossFactEdges;

  [[nodiscard]] std::string str(const std::string &Indent = "") const;
  DOTFactSubGraph *getOrCreateFactSG(unsigned FactId, std::string &Label);
  // TODO: pass the actual lambda EF name and value as parameter from DOTGraph
  [[nodiscard]] std::string generateLambdaSG(const std::string &Indent = "") const;
  void createLayoutCFNodes();
  void createLayoutFactNodes();
  void createLayoutFactEdges();
};

std::ostream &operator<<(std::ostream &Os,
                         const DOTFunctionSubGraph &FunctionSg);

template <typename D> struct DOTGraph {
  std::string Label;
  std::map<std::string, DOTFunctionSubGraph> Functions;
  std::set<DOTEdge> InterCfEdges;
  std::set<DOTEdge> InterLambdaEdges;
  std::set<DOTEdge> InterFactEdges;

  DOTGraph() = default;

  unsigned getFactID(D Fact) {
    unsigned Id = 0;
    if (DtoFactId.count(Fact)) {
      Id = DtoFactId[Fact];
    } else {
      Id = FactIdCount;
      DtoFactId[Fact] = FactIdCount++;
    }
    return Id;
  }

  bool containsFactSG(const std::string& FName, unsigned FactId) {
    if (Functions.count(FName)) {
      if (Functions[FName].Facts.count(FactId)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] std::string str() const {
    std::string Indent = "  ";
    std::string Str = "digraph {\n" + Indent + "label=\"" + Label + "\"\n";
    // Print function subgraphs
    Str += '\n' + Indent + "// Function sub graphs\n";
    for (auto Fsg : Functions) {
      Fsg.second.createLayoutCFNodes();
      Fsg.second.createLayoutFactNodes();
      Fsg.second.createLayoutFactEdges();
      Str += Fsg.second.str(Indent) + "\n\n";
    }

    // Print inter control flow Edges
    Str += Indent + "// Inter-procedural control flow Edges\n" + Indent +
           DOTConfig::cfInterEdgeAttr() + '\n';
    for (const DOTEdge& E : InterCfEdges) {
      Str += E.str(Indent) + '\n';
    }

    // Print inter lambda Edges
    Str += '\n' + Indent + "// Inter-procedural lambda Edges\n" + Indent +
           DOTConfig::lambdaInterEdgeAttr() + '\n';
    for (const DOTEdge& E : InterLambdaEdges) {
      Str += E.str(Indent) + '\n';
    }

    // Print inter fact Edges
    Str += '\n' + Indent + "// Inter-procedural fact Edges\n" + Indent +
           DOTConfig::factInterEdgeAttr() + '\n';
    for (const DOTEdge& E : InterFactEdges) {
      Str += E.str(Indent) + '\n';
    }
    return Str + '}';
  }

  friend std::ostream &operator<<(std::ostream &Os, const DOTGraph<D> &Graph) {
    return Os << Graph.str();
  }

private:
  // We introduce a fact-ID for data-flow facts D since only statements N have
  // an ID
  unsigned FactIdCount = 1;
  std::map<D, unsigned> DtoFactId;
};

} // namespace psr

#endif
