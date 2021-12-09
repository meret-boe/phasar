/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * IDESolver.h
 *
 *  Created on: 04.08.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SOLVER_IDESOLVER_H
#define PHASAR_PHASARLLVM_DATAFLOWSOLVER_IFDSIDE_SOLVER_IDESOLVER_H

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>

#include "nlohmann/json.hpp"

#include "boost/algorithm/string/trim.hpp"

#include "llvm/Support/raw_ostream.h"

#include "phasar/Config/Configuration.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/FlowEdgeFunctionCache.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/FlowFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IDETabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IFDSTabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/JoinLattice.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IFDSSolverTest.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/IFDSToIDETabulationProblem.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/JoinHandlingNode.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/JumpFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/LinkedNode.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/PathEdge.h"
#include "phasar/PhasarLLVM/Domain/AnalysisDomain.h"
#include "phasar/PhasarLLVM/Utils/DOTGraph.h"
#include "phasar/Utils/LLVMShorthands.h"
#include "phasar/Utils/Logger.h"
#include "phasar/Utils/PAMMMacros.h"
#include "phasar/Utils/Table.h"

namespace psr {

// Forward declare the Transformation
template <typename AnalysisDomainTy, typename Container>
class IFDSToIDETabulationProblem;

struct NoIFDSExtension {};

template <typename AnalysisDomainTy, typename Container> struct IFDSExtension {
  using BaseAnalysisDomain = typename AnalysisDomainTy::BaseAnalysisDomain;

  IFDSExtension(IFDSTabulationProblem<BaseAnalysisDomain, Container> &Problem)
      : TransformedProblem(
            std::make_unique<
                IFDSToIDETabulationProblem<BaseAnalysisDomain, Container>>(
                Problem)) {}

  std::unique_ptr<IFDSToIDETabulationProblem<BaseAnalysisDomain, Container>>
      TransformedProblem;
};

/// Solves the given IDETabulationProblem as described in the 1996 paper by
/// Sagiv, Horwitz and Reps. To solve the problem, call solve(). Results
/// can then be queried by using resultAt() and resultsAt().
template <typename AnalysisDomainTy,
          typename Container = std::set<typename AnalysisDomainTy::d_t>,
          bool = IsAnalysisDomainExtensions<AnalysisDomainTy>::value>
class IDESolver
    : protected std::conditional_t<
          IsAnalysisDomainExtensions<AnalysisDomainTy>::value,
          IFDSExtension<AnalysisDomainTy, Container>, NoIFDSExtension> {
public:
  using ProblemTy = IDETabulationProblem<AnalysisDomainTy, Container>;
  using container_type = typename ProblemTy::container_type;
  using FlowFunctionPtrType = typename ProblemTy::FlowFunctionPtrType;
  using EdgeFunctionPtrType = typename ProblemTy::EdgeFunctionPtrType;

  using l_t = typename AnalysisDomainTy::l_t;
  using n_t = typename AnalysisDomainTy::n_t;
  using i_t = typename AnalysisDomainTy::i_t;
  using d_t = typename AnalysisDomainTy::d_t;
  using f_t = typename AnalysisDomainTy::f_t;
  using t_t = typename AnalysisDomainTy::t_t;
  using v_t = typename AnalysisDomainTy::v_t;

  IDESolver(IDETabulationProblem<AnalysisDomainTy, Container> &Problem)
      : IDEProblem(Problem), ZeroValue(Problem.getZeroValue()),
        ICF(Problem.getICFG()), SolverConfig(Problem.getIFDSIDESolverConfig()),
        CachedFlowEdgeFunctions(Problem), AllTop(Problem.allTopFunction()),
        JumpFn(std::make_shared<JumpFunctions<AnalysisDomainTy, Container>>(
            AllTop, IDEProblem)),
        InitialSeeds(Problem.initialSeeds()) {}

  IDESolver(const IDESolver &) = delete;
  IDESolver &operator=(const IDESolver &) = delete;
  IDESolver(IDESolver &&) = delete;
  IDESolver &operator=(IDESolver &&) = delete;

  virtual ~IDESolver() = default;

  nlohmann::json getAsJson() {
    using TableCell = typename Table<n_t, d_t, l_t>::Cell;
    const static std::string DataFlowID = "DataFlow";
    nlohmann::json J;
    auto Results = this->Valtab.cellSet();
    if (Results.empty()) {
      J[DataFlowID] = "EMPTY";
    } else {
      std::vector<TableCell> Cells(Results.begin(), Results.end());
      sort(Cells.begin(), Cells.end(), [](TableCell A, TableCell B) {
        return A.getRowKey() < B.getRowKey();
      });
      n_t Curr;
      for (unsigned I = 0; I < Cells.size(); ++I) {
        Curr = Cells[I].getRowKey();
        std::string N = IDEProblem.NtoString(Cells[I].getRowKey());
        boost::algorithm::trim(N);
        std::string Node =
            ICF->getFunctionName(ICF->getFunctionOf(Curr)) + "::" + N;
        J[DataFlowID][Node];
        std::string Fact = IDEProblem.DtoString(Cells[I].getColumnKey());
        boost::algorithm::trim(Fact);
        std::string Value = IDEProblem.LtoString(Cells[I].getValue());
        boost::algorithm::trim(Value);
        J[DataFlowID][Node]["Facts"] += {Fact, Value};
      }
    }
    return J;
  }

  /// \brief Runs the solver on the configured problem. This can take some time.
  virtual void solve() {
    PAMM_GET_INSTANCE;
    REG_COUNTER("Gen facts", 0, PAMM_SEVERITY_LEVEL::Core);
    REG_COUNTER("Kill facts", 0, PAMM_SEVERITY_LEVEL::Core);
    REG_COUNTER("Summary-reuse", 0, PAMM_SEVERITY_LEVEL::Core);
    REG_COUNTER("Intra Path Edges", 0, PAMM_SEVERITY_LEVEL::Core);
    REG_COUNTER("Inter Path Edges", 0, PAMM_SEVERITY_LEVEL::Core);
    REG_COUNTER("FF Queries", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("EF Queries", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Value Propagation", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Value Computation", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("SpecialSummary-FF Application", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("SpecialSummary-EF Queries", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("JumpFn Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Process Call", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Process Normal", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Process Exit", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("[Calls] getPointsToSet", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_HISTOGRAM("Data-flow facts", PAMM_SEVERITY_LEVEL::Full);
    REG_HISTOGRAM("Points-to", PAMM_SEVERITY_LEVEL::Full);

    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                      << "IDE solver is solving the specified problem";
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "Submit initial seeds, construct exploded super graph");
    // computations starting here
    START_TIMER("DFA Phase I", PAMM_SEVERITY_LEVEL::Full);
    // We start our analysis and construct exploded supergraph
    submitInitialSeeds();
    STOP_TIMER("DFA Phase I", PAMM_SEVERITY_LEVEL::Full);
    if (SolverConfig.computeValues()) {
      START_TIMER("DFA Phase II", PAMM_SEVERITY_LEVEL::Full);
      // Computing the final values for the edge functions
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Compute the final values according to the edge functions");
      computeValues();
      STOP_TIMER("DFA Phase II", PAMM_SEVERITY_LEVEL::Full);
    }
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO) << "Problem solved");
    if constexpr (PammCurrSevLevel >= PammSeverityLevel::Core) {
      computeAndPrintStatistics();
    }
    if (SolverConfig.emitESG()) {
      emitESGAsDot();
    }
  }

  /// Returns the L-type result for the given value at the given statement.
  [[nodiscard]] virtual l_t resultAt(n_t Stmt, d_t Value) {
    return Valtab.get(Stmt, Value);
  }

  /// Returns the L-type result at the given statement for the given data-flow
  /// fact while respecting LLVM's SSA semantics.
  ///
  /// An example: when a value is loaded and the location loaded from, here
  /// variable 'i', is a data-flow fact that holds, then the loaded value '%0'
  /// will usually be generated and also holds. However, due to the underlying
  /// theory (and respective implementation) this load instruction causes the
  /// loaded value to be generated and thus, it will be valid only AFTER the
  /// load instruction, i.e., at the successor instruction.
  ///
  ///   %0 = load i32, i32* %i, align 4
  ///
  /// This result accessor function returns the results at the successor
  /// instruction(s) reflecting that the expression on the left-hand side holds
  /// if the expression on the right-hand side holds.
  template <typename NTy = n_t>
  [[nodiscard]] typename std::enable_if_t<
      std::is_same_v<std::remove_reference_t<NTy>, llvm::Instruction *>, l_t>
  resultAtInLLVMSSA(NTy Stmt, d_t Value) {
    if (Stmt->getType()->isVoidTy()) {
      return Valtab.get(Stmt, Value);
    }
    assert(Stmt->getNextNode() && "Expected to find a valid successor node!");
    return Valtab.get(Stmt->getNextNode(), Value);
  }

  /// Returns the resulting environment for the given statement.
  /// The artificial zero value can be automatically stripped.
  /// TOP values are never returned.
  [[nodiscard]] virtual std::unordered_map<d_t, l_t>
  resultsAt(n_t Stmt, bool StripZero = false) /*TODO const*/ {
    std::unordered_map<d_t, l_t> Result = Valtab.row(Stmt);
    if (StripZero) {
      for (auto It = Result.begin(); It != Result.end();) {
        if (IDEProblem.isZeroValue(It->first)) {
          It = Result.erase(It);
        } else {
          ++It;
        }
      }
    }
    return Result;
  }

  /// Returns the data-flow results at the given statement while respecting
  /// LLVM's SSA semantics.
  ///
  /// An example: when a value is loaded and the location loaded from, here
  /// variable 'i', is a data-flow fact that holds, then the loaded value '%0'
  /// will usually be generated and also holds. However, due to the underlying
  /// theory (and respective implementation) this load instruction causes the
  /// loaded value to be generated and thus, it will be valid only AFTER the
  /// load instruction, i.e., at the successor instruction.
  ///
  ///   %0 = load i32, i32* %i, align 4
  ///
  /// This result accessor function returns the results at the successor
  /// instruction(s) reflecting that the expression on the left-hand side holds
  /// if the expression on the right-hand side holds.
  template <typename NTy = n_t>
  [[nodiscard]] typename std::enable_if_t<
      std::is_same_v<std::remove_reference_t<NTy>, llvm::Instruction *>,
      std::unordered_map<d_t, l_t>>
  resultsAtInLLVMSSA(NTy Stmt, bool StripZero = false) {
    std::unordered_map<d_t, l_t> Result = [this, Stmt]() {
      if (Stmt->getType()->isVoidTy()) {
        return Valtab.row(Stmt);
      }
      return Valtab.row(Stmt->getNextNode());
    }();
    if (StripZero) {
      // TODO: replace with std::erase_if (C++20)
      for (auto It = Result.begin(); It != Result.end();) {
        if (IDEProblem.isZeroValue(It->first)) {
          It = Result.erase(It);
        } else {
          ++It;
        }
      }
    }
    return Result;
  }

  virtual void emitTextReport(std::ostream &OS = std::cout) {
    IDEProblem.emitTextReport(getSolverResults(), OS);
  }

  virtual void emitGraphicalReport(std::ostream &OS = std::cout) {
    IDEProblem.emitGraphicalReport(getSolverResults(), OS);
  }

  virtual void dumpResults(std::ostream &OS = std::cout) {
    PAMM_GET_INSTANCE;
    START_TIMER("DFA IDE Result Dumping", PAMM_SEVERITY_LEVEL::Full);
    OS << "\n***************************************************************\n"
       << "*                  Raw IDESolver results                      *\n"
       << "***************************************************************\n";
    auto Cells = this->Valtab.cellVec();
    if (Cells.empty()) {
      OS << "No results computed!" << std::endl;
    } else {
      LlvmValueIdLess LlvmIdLess;
      std::sort(
          Cells.begin(), Cells.end(),
          [&LlvmIdLess](const auto &A, const auto &B) {
            if constexpr (std::is_same_v<n_t, const llvm::Instruction *>) {
              return LlvmIdLess(A.getRowKey(), B.getRowKey());
            } else {
              // If non-LLVM IR is used
              return A.getRowKey() < B.getRowKey();
            }
          });
      n_t Prev = n_t{};
      n_t Curr = n_t{};
      f_t PrevFn = f_t{};
      f_t CurrFn = f_t{};
      for (unsigned I = 0; I < Cells.size(); ++I) {
        Curr = Cells[I].getRowKey();
        CurrFn = ICF->getFunctionOf(Curr);
        if (PrevFn != CurrFn) {
          PrevFn = CurrFn;
          OS << "\n\n============ Results for function '" +
                    ICF->getFunctionName(CurrFn) + "' ============\n";
        }
        if (Prev != Curr) {
          Prev = Curr;
          std::string NString = IDEProblem.NtoString(Curr);
          std::string Line(NString.size(), '-');
          OS << "\n\nN: " << NString << "\n---" << Line << '\n';
        }
        OS << "\tD: " << IDEProblem.DtoString(Cells[I].getColumnKey())
           << " | V: " << IDEProblem.LtoString(Cells[I].getValue()) << '\n';
      }
    }
    OS << '\n';
    STOP_TIMER("DFA IDE Result Dumping", PAMM_SEVERITY_LEVEL::Full);
  }

  void dumpAllInterPathEdges() {
    std::cout << "COMPUTED INTER PATH EDGES" << std::endl;
    auto Interpe = this->ComputedInterPathEdges.cellSet();
    for (const auto &Cell : Interpe) {
      std::cout << "FROM" << std::endl;
      IDEProblem.printNode(std::cout, Cell.getRowKey());
      std::cout << "TO" << std::endl;
      IDEProblem.printNode(std::cout, Cell.getColumnKey());
      std::cout << "FACTS" << std::endl;
      for (const auto &Fact : Cell.getValue()) {
        std::cout << "fact" << std::endl;
        IDEProblem.printFlowFact(std::cout, Fact.first);
        std::cout << "produces" << std::endl;
        for (const auto &Out : Fact.second) {
          IDEProblem.printFlowFact(std::cout, Out);
        }
      }
    }
  }

  void dumpAllIntraPathEdges() {
    std::cout << "COMPUTED INTRA PATH EDGES" << std::endl;
    auto Intrape = this->ComputedIntraPathEdges.cellSet();
    for (auto &Cell : Intrape) {
      std::cout << "FROM" << std::endl;
      IDEProblem.printNode(std::cout, Cell.getRowKey());
      std::cout << "TO" << std::endl;
      IDEProblem.printNode(std::cout, Cell.getColumnKey());
      std::cout << "FACTS" << std::endl;
      for (auto &Fact : Cell.getValue()) {
        std::cout << "fact" << std::endl;
        IDEProblem.printFlowFact(std::cout, Fact.first);
        std::cout << "produces" << std::endl;
        for (auto &Out : Fact.second) {
          IDEProblem.printFlowFact(std::cout, Out);
        }
      }
    }
  }

  SolverResults<n_t, d_t, l_t> getSolverResults() {
    return SolverResults<n_t, d_t, l_t>(this->Valtab,
                                        IDEProblem.getZeroValue());
  }

protected:
  // have a shared point to allow for a copy constructor of IDESolver
  IDETabulationProblem<AnalysisDomainTy, Container> &IDEProblem;
  d_t ZeroValue;
  const i_t *ICF;
  IFDSIDESolverConfig &SolverConfig;
  unsigned PathEdgeCount = 0;

  FlowEdgeFunctionCache<AnalysisDomainTy, Container> CachedFlowEdgeFunctions;

  Table<n_t, n_t, std::map<d_t, Container>> ComputedIntraPathEdges;

  Table<n_t, n_t, std::map<d_t, Container>> ComputedInterPathEdges;

  EdgeFunctionPtrType AllTop;

  std::shared_ptr<JumpFunctions<AnalysisDomainTy, Container>> JumpFn;

  std::map<std::tuple<n_t, d_t, n_t, d_t>, std::vector<EdgeFunctionPtrType>>
      IntermediateEdgeFunctions;

  // stores summaries that were queried before they were computed
  // see CC 2010 paper by Naeem, Lhotak and Rodriguez
  Table<n_t, d_t, Table<n_t, d_t, EdgeFunctionPtrType>> Endsummarytab;

  // edges going along calls
  // see CC 2010 paper by Naeem, Lhotak and Rodriguez
  Table<n_t, d_t, std::map<n_t, Container>> Incomingtab;

  // stores the return sites (inside callers) to which we have unbalanced
  // returns if SolverConfig.followReturnPastSeeds is enabled
  std::set<n_t> UnbalancedRetSites;

  std::map<n_t, std::set<d_t>> InitialSeeds;

  Table<n_t, d_t, l_t> Valtab;

  std::map<std::pair<n_t, d_t>, size_t> FSummaryReuse;

  // When transforming an IFDSTabulationProblem into an IDETabulationProblem,
  // we need to allocate dynamically, otherwise the objects lifetime runs out
  // - as a modifiable r-value reference created here that should be stored in
  // a modifiable l-value reference within the IDESolver implementation leads
  // to (massive) undefined behavior (and nightmares):
  // https://stackoverflow.com/questions/34240794/understanding-the-warning-binding-r-value-to-l-value-reference
  template <typename IFDSAnalysisDomainTy,
            typename = std::enable_if_t<
                IsAnalysisDomainExtensions<AnalysisDomainTy>::value,
                IFDSAnalysisDomainTy>>
  IDESolver(IFDSTabulationProblem<IFDSAnalysisDomainTy, Container> &Problem)
      : IFDSExtension<AnalysisDomainTy, Container>(Problem),
        IDEProblem(*this->TransformedProblem),
        ZeroValue(IDEProblem.getZeroValue()), ICF(IDEProblem.getICFG()),
        SolverConfig(IDEProblem.getIFDSIDESolverConfig()),
        CachedFlowEdgeFunctions(IDEProblem),
        AllTop(IDEProblem.allTopFunction()),
        JumpFn(std::make_shared<JumpFunctions<AnalysisDomainTy, Container>>(
            AllTop, IDEProblem)),
        InitialSeeds(IDEProblem.initialSeeds()) {}

  /// Lines 13-20 of the algorithm; processing a call site in the caller's
  /// context.
  ///
  /// For each possible callee, registers incoming call edges.
  /// Also propagates call-to-return flows and summarized callee flows within
  /// the caller.
  ///
  /// 	The following cases must be considered and handled:
  ///		1. Process as usual and just process the call
  ///		2. Create a new summary for that function (which shall be done
  ///       by the problem)
  ///		3. Just use an existing summary provided by the problem
  ///		4. If a special function is called, use a special summary
  ///       function
  ///
  /// @param edge an edge whose target node resembles a method call
  ///
  virtual void processCall(const PathEdge<n_t, d_t> Edge) { //NOLINT
    PAMM_GET_INSTANCE;
    INC_COUNTER("Process Call", 1, PAMM_SEVERITY_LEVEL::Full);
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Process call at target: "
                  << IDEProblem.NtoString(Edge.getTarget()));
    d_t D1 = Edge.factAtSource();
    n_t N = Edge.getTarget(); // a call node; line 14...
    d_t D2 = Edge.factAtTarget();
    EdgeFunctionPtrType F = jumpFunction(Edge);
    const std::set<n_t> ReturnSiteNs = ICF->getReturnSitesOfCallAt(N);
    const std::set<f_t> Callees = ICF->getCalleesOfCallAt(N);

    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG) << "Possible callees:";
        for (auto callee
             : Callees) {
          BOOST_LOG_SEV(lg::get(), DEBUG) << "  " << callee->getName().str();
        } BOOST_LOG_SEV(lg::get(), DEBUG)
        << "Possible return sites:";
        for (auto ret
             : ReturnSiteNs) {
          BOOST_LOG_SEV(lg::get(), DEBUG) << "  " << IDEProblem.NtoString(ret);
        } BOOST_LOG_SEV(lg::get(), DEBUG)
        << ' ');

    // for each possible callee
    for (f_t SCalledProcN : Callees) { // still line 14
      // check if a special summary for the called procedure exists
      FlowFunctionPtrType SpecialSum =
          CachedFlowEdgeFunctions.getSummaryFlowFunction(N, SCalledProcN);
      // if a special summary is available, treat this as a normal flow
      // and use the summary flow and edge functions
      if (SpecialSum) {
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Found and process special summary");
        for (n_t ReturnSiteN : ReturnSiteNs) {
          container_type Res = computeSummaryFlowFunction(SpecialSum, D1, D2);
          INC_COUNTER("SpecialSummary-FF Application", 1,
                      PAMM_SEVERITY_LEVEL::Full);
          ADD_TO_HISTOGRAM("Data-flow facts", res.size(), 1,
                           PAMM_SEVERITY_LEVEL::Full);
          saveEdges(N, ReturnSiteN, D2, Res, false);
          for (d_t D3 : Res) {
            EdgeFunctionPtrType SumEdgFnE =
                CachedFlowEdgeFunctions.getSummaryEdgeFunction(N, D2,
                                                               ReturnSiteN, D3);
            INC_COUNTER("SpecialSummary-EF Queries", 1,
                        PAMM_SEVERITY_LEVEL::Full);
            LOG_IF_ENABLE(
                BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "Queried Summary Edge Function: " << SumEdgFnE->str();
                BOOST_LOG_SEV(lg::get(), DEBUG)
                << "Compose: " << SumEdgFnE->str() << " * " << F->str();
                BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
            propagate(D1, ReturnSiteN, D3, F->composeWith(SumEdgFnE), N, false);
          }
        }
      } else {
        // compute the call-flow function
        FlowFunctionPtrType Function =
            CachedFlowEdgeFunctions.getCallFlowFunction(N, SCalledProcN);
        INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        container_type Res = computeCallFlowFunction(Function, D1, D2);
        ADD_TO_HISTOGRAM("Data-flow facts", res.size(), 1,
                         PAMM_SEVERITY_LEVEL::Full);
        // for each callee's start point(s)
        std::set<n_t> StartPointsOf = ICF->getStartPointsOf(SCalledProcN);
        if (StartPointsOf.empty()) {
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                            << "Start points of '" +
                                   ICF->getFunctionName(SCalledProcN) +
                                   "' currently not available!";
                        BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        }
        // if startPointsOf is empty, the called function is a declaration
        for (n_t SP : StartPointsOf) {
          saveEdges(N, SP, D2, Res, true);
          // for each result node of the call-flow function
          for (d_t D3 : Res) {
            using TableCell =
                typename Table<n_t, d_t, EdgeFunctionPtrType>::Cell;
            // create initial self-loop
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Create initial self-loop with D: "
                          << IDEProblem.DtoString(D3));
            propagate(D3, SP, D3, EdgeIdentity<l_t>::getInstance(), N,
                      false); // line 15
            // register the fact that <sp,d3> has an incoming edge from <n,d2>
            // line 15.1 of Naeem/Lhotak/Rodriguez
            addIncoming(SP, D3, N, D2);
            // line 15.2, copy to avoid concurrent modification exceptions by
            // other threads
            // const std::set<TableCell> endSumm(endSummary(sP, d3));
            // std::cout << "ENDSUMM" << std::endl;
            // std::cout << "Size: " << endSumm.size() << std::endl;
            // std::cout << "sP: " << IDEProblem.NtoString(sP)
            //           << "\nd3: " << IDEProblem.DtoString(d3)
            //           << std::endl;
            // printEndSummaryTab();
            // still line 15.2 of Naeem/Lhotak/Rodriguez
            // for each already-queried exit value <eP,d4> reachable from
            // <sP,d3>, create new caller-side jump functions to the return
            // sites because we have observed a potentially new incoming
            // edge into <sP,d3>
            for (const TableCell Entry : endSummary(SP, D3)) {
              n_t EP = Entry.getRowKey();
              d_t D4 = Entry.getColumnKey();
              EdgeFunctionPtrType FCalleeSummary = Entry.getValue();
              // for each return site
              for (n_t RetSiteN : ReturnSiteNs) {
                // compute return-flow function
                FlowFunctionPtrType RetFunction =
                    CachedFlowEdgeFunctions.getRetFlowFunction(N, SCalledProcN,
                                                               EP, RetSiteN);
                INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
                const container_type ReturnedFacts = computeReturnFlowFunction(
                    RetFunction, D3, D4, N, Container{D2});
                ADD_TO_HISTOGRAM("Data-flow facts", returnedFacts.size(), 1,
                                 PAMM_SEVERITY_LEVEL::Full);
                saveEdges(EP, RetSiteN, D4, ReturnedFacts, true);
                // for each target value of the function
                for (d_t D5 : ReturnedFacts) {
                  // update the caller-side summary function
                  // get call edge function
                  EdgeFunctionPtrType F4 =
                      CachedFlowEdgeFunctions.getCallEdgeFunction(
                          N, D2, SCalledProcN, D3);
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                << "Queried Call Edge Function: " << F4->str());
                  // get return edge function
                  EdgeFunctionPtrType F5 =
                      CachedFlowEdgeFunctions.getReturnEdgeFunction(
                          N, SCalledProcN, EP, D4, RetSiteN, D5);
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                << "Queried Return Edge Function: "
                                << F5->str());
                  if (SolverConfig.emitESG()) {
                    for (auto SP : ICF->getStartPointsOf(SCalledProcN)) {
                      IntermediateEdgeFunctions[std::make_tuple(N, D2, SP, D3)]
                          .push_back(F4);
                    }
                    IntermediateEdgeFunctions[std::make_tuple(EP, D4, RetSiteN,
                                                              D5)]
                        .push_back(F5);
                  }
                  INC_COUNTER("EF Queries", 2, PAMM_SEVERITY_LEVEL::Full);
                  // compose call * calleeSummary * return edge functions
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                    << "Compose: " << F5->str() << " * "
                                    << FCalleeSummary->str() << " * "
                                    << F4->str();
                                BOOST_LOG_SEV(lg::get(), DEBUG)
                                << "         (return * calleeSummary * call)");
                  EdgeFunctionPtrType FPrime =
                      F4->composeWith(FCalleeSummary)->composeWith(F5);
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                    << "       = " << FPrime->str();
                                BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
                  d_t D5RestoredCtx = restoreContextOnReturnedFact(N, D2, D5);
                  // propagte the effects of the entire call
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                    << "Compose: " << FPrime->str() << " * "
                                    << F->str();
                                BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
                  propagate(D1, RetSiteN, D5RestoredCtx,
                            F->composeWith(FPrime), N, false);
                }
              }
            }
          }
        }
      }
    }
    // line 17-19 of Naeem/Lhotak/Rodriguez
    // process intra-procedural flows along call-to-return flow functions
    for (n_t ReturnSiteN : ReturnSiteNs) {
      FlowFunctionPtrType CallToReturnFlowFunction =
          CachedFlowEdgeFunctions.getCallToRetFlowFunction(N, ReturnSiteN,
                                                           Callees);
      INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
      container_type ReturnFacts =
          computeCallToReturnFlowFunction(CallToReturnFlowFunction, D1, D2);
      ADD_TO_HISTOGRAM("Data-flow facts", returnFacts.size(), 1,
                       PAMM_SEVERITY_LEVEL::Full);
      saveEdges(N, ReturnSiteN, D2, ReturnFacts, false);
      for (d_t D3 : ReturnFacts) {
        EdgeFunctionPtrType EdgeFnE =
            CachedFlowEdgeFunctions.getCallToRetEdgeFunction(N, D2, ReturnSiteN,
                                                             D3, Callees);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Queried Call-to-Return Edge Function: "
                      << EdgeFnE->str());
        if (SolverConfig.emitESG()) {
          IntermediateEdgeFunctions[std::make_tuple(N, D2, ReturnSiteN, D3)]
              .push_back(EdgeFnE);
        }
        INC_COUNTER("EF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        auto FPrime = F->composeWith(EdgeFnE);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Compose: " << EdgeFnE->str() << " * " << F->str()
                          << " = " << FPrime->str();
                      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        propagate(D1, ReturnSiteN, D3, FPrime, N, false);
      }
    }
  }

  /// Lines 33-37 of the algorithm.
  /// Simply propagate normal, intra-procedural flows.
  /// @param edge
  ///
  virtual void processNormalFlow(const PathEdge<n_t, d_t> Edge) { //NOLINT
    PAMM_GET_INSTANCE;
    INC_COUNTER("Process Normal", 1, PAMM_SEVERITY_LEVEL::Full);
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Process normal at target: "
                  << IDEProblem.NtoString(Edge.getTarget()));
    d_t D1 = Edge.factAtSource();
    n_t N = Edge.getTarget();
    d_t D2 = Edge.factAtTarget();
    EdgeFunctionPtrType F = jumpFunction(Edge);
    for (const auto Fn : ICF->getSuccsOf(N)) {
      FlowFunctionPtrType FlowFunction =
          CachedFlowEdgeFunctions.getNormalFlowFunction(N, Fn);
      INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
      const container_type Res =
          computeNormalFlowFunction(FlowFunction, D1, D2);
      ADD_TO_HISTOGRAM("Data-flow facts", res.size(), 1,
                       PAMM_SEVERITY_LEVEL::Full);
      saveEdges(N, Fn, D2, Res, false);
      for (d_t D3 : Res) {
        EdgeFunctionPtrType G =
            CachedFlowEdgeFunctions.getNormalEdgeFunction(N, D2, Fn, D3);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Queried Normal Edge Function: " << G->str());
        EdgeFunctionPtrType Fprime = F->composeWith(G);
        if (SolverConfig.emitESG()) {
          IntermediateEdgeFunctions[std::make_tuple(N, D2, Fn, D3)].push_back(
              G);
        }
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Compose: " << G->str() << " * " << F->str()
                          << " = " << Fprime->str();
                      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        INC_COUNTER("EF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        propagate(D1, Fn, D3, Fprime, nullptr, false);
      }
    }
  }

  void propagateValueAtStart(const std::pair<n_t, d_t> NAndD, n_t N) {
    PAMM_GET_INSTANCE;
    d_t D = NAndD.second;
    f_t P = ICF->getFunctionOf(N);
    for (const n_t C : ICF->getCallsFromWithin(P)) {
      auto LookupResults = JumpFn->forwardLookup(D, C);
      if (!LookupResults) {
        continue;
      }
      for (auto Entry : LookupResults->get()) {
        d_t DPrime = Entry.first;
        EdgeFunctionPtrType FPrime = Entry.second;
        n_t SP = N;
        l_t Value = val(SP, D);
        INC_COUNTER("Value Propagation", 1, PAMM_SEVERITY_LEVEL::Full);
        propagateValue(C, DPrime, FPrime->computeTarget(Value));
      }
    }
  }

  void propagateValueAtCall(const std::pair<n_t, d_t> NAndD, n_t N) {
    PAMM_GET_INSTANCE;
    d_t D = NAndD.second;
    for (const f_t Q : ICF->getCalleesOfCallAt(N)) {
      FlowFunctionPtrType CallFlowFunction =
          CachedFlowEdgeFunctions.getCallFlowFunction(N, Q);
      INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
      for (const d_t DPrime : CallFlowFunction->computeTargets(D)) {
        EdgeFunctionPtrType EdgeFn =
            CachedFlowEdgeFunctions.getCallEdgeFunction(N, D, Q, DPrime);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Queried Call Edge Function: " << EdgeFn->str());
        if (SolverConfig.emitESG()) {
          for (const auto SP : ICF->getStartPointsOf(Q)) {
            IntermediateEdgeFunctions[std::make_tuple(N, D, SP, DPrime)]
                .push_back(EdgeFn);
          }
        }
        INC_COUNTER("EF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        for (const n_t StartPoint : ICF->getStartPointsOf(Q)) {
          INC_COUNTER("Value Propagation", 1, PAMM_SEVERITY_LEVEL::Full);
          propagateValue(StartPoint, DPrime, EdgeFn->computeTarget(val(N, D)));
        }
      }
    }
  }

  void propagateValue(n_t NHashN, d_t NHashD, const l_t &L) {
    l_t ValNHash = val(NHashN, NHashD);
    l_t LPrime = joinValueAt(NHashN, NHashD, ValNHash, L);
    if (!(LPrime == ValNHash)) {
      setVal(NHashN, NHashD, std::move(LPrime));
      valuePropagationTask(std::pair<n_t, d_t>(NHashN, NHashD));
    }
  }

  l_t val(n_t NHashN, d_t NHashD) {
    if (Valtab.contains(NHashN, NHashD)) {
      return Valtab.get(NHashN, NHashD);
    }        // implicitly initialized to top; see line [1] of Fig. 7 in SRH96 paper
      return IDEProblem.topElement();
   
  }

  void setVal(n_t NHashN, d_t NHashD, l_t L) {
    LOG_IF_ENABLE([&]() {
      BOOST_LOG_SEV(lg::get(), DEBUG)
          << "Function : " << ICF->getFunctionOf(NHashN)->getName().str();
      BOOST_LOG_SEV(lg::get(), DEBUG)
          << "Inst.    : " << IDEProblem.NtoString(NHashN);
      BOOST_LOG_SEV(lg::get(), DEBUG)
          << "Fact     : " << IDEProblem.DtoString(NHashD);
      BOOST_LOG_SEV(lg::get(), DEBUG)
          << "Value    : " << IDEProblem.LtoString(L);
      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
    }());
    // TOP is the implicit default value which we do not need to store.
    // if (l == IDEProblem.topElement()) {
    // do not store top values
    // valtab.remove(nHashN, nHashD);
    // } else {
    Valtab.insert(NHashN, NHashD, std::move(L));
    // }
  }

  EdgeFunctionPtrType jumpFunction(const PathEdge<n_t, d_t> Edge) { //NOLINT
    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG) << " ";
        BOOST_LOG_SEV(lg::get(), DEBUG) << "JumpFunctions Forward-Lookup:";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "   Source D: " << IDEProblem.DtoString(Edge.factAtSource());
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "   Target N: " << IDEProblem.NtoString(Edge.getTarget());
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "   Target D: " << IDEProblem.DtoString(Edge.factAtTarget()));

    auto FwdLookupRes =
        JumpFn->forwardLookup(Edge.factAtSource(), Edge.getTarget());
    if (FwdLookupRes) {
      auto &Ref = FwdLookupRes->get();
      if (auto Find = std::find_if(Ref.begin(), Ref.end(),
                                   [Edge](const auto &Pair) {
                                     return Edge.factAtTarget() == Pair.first;
                                   });
          Find != Ref.end()) {
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "  => EdgeFn: " << Find->second->str();
                      BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
        return Find->second;
      }
    }
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "  => EdgeFn: " << AllTop->str();
                  BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
    // JumpFn initialized to all-top, see line [2] in SRH96 paper
    return AllTop;
  }

  void addEndSummary(n_t SP, d_t D1, n_t EP, d_t D2, EdgeFunctionPtrType F) {
    // note: at this point we don't need to join with a potential previous f
    // because f is a jump function, which is already properly joined
    // within propagate(..)
    Endsummarytab.get(SP, D1).insert(EP, D2, std::move(F));
  }

  // should be made a callable at some point
  void pathEdgeProcessingTask(const PathEdge<n_t, d_t> Edge) {
    PAMM_GET_INSTANCE;
    INC_COUNTER("JumpFn Construction", 1, PAMM_SEVERITY_LEVEL::Full);
    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "-------------------------------------------- " << PathEdgeCount
            << ". Path Edge --------------------------------------------";
        BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "Process " << PathEdgeCount << ". path edge:";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "< D source: " << IDEProblem.DtoString(Edge.factAtSource()) << " ;";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "  N target: " << IDEProblem.NtoString(Edge.getTarget()) << " ;";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "  D target: " << IDEProblem.DtoString(Edge.factAtTarget()) << " >";
        BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');

    if (!ICF->isCallSite(Edge.getTarget())) {
      if (ICF->isExitInst(Edge.getTarget())) {
        processExit(Edge);
      }
      if (!ICF->getSuccsOf(Edge.getTarget()).empty()) {
        processNormalFlow(Edge);
      }
    } else {
      processCall(Edge);
    }
  }

  // should be made a callable at some point
  void valuePropagationTask(const std::pair<n_t, d_t> NAndD) {
    n_t N = NAndD.first;
    // our initial seeds are not necessarily method-start points but here they
    // should be treated as such the same also for unbalanced return sites in
    // an unbalanced problem
    if (ICF->isStartPoint(N) || InitialSeeds.count(N) ||
        UnbalancedRetSites.count(N)) {
      propagateValueAtStart(NAndD, N);
    }
    if (ICF->isCallSite(N)) {
      propagateValueAtCall(NAndD, N);
    }
  }

  // should be made a callable at some point
  void valueComputationTask(const std::vector<n_t> &Values) {
    PAMM_GET_INSTANCE;
    for (n_t N : Values) {
      for (n_t SP : ICF->getStartPointsOf(ICF->getFunctionOf(N))) {
        using TableCell = typename Table<d_t, d_t, EdgeFunctionPtrType>::Cell;
        Table<d_t, d_t, EdgeFunctionPtrType> LookupByTarget;
        LookupByTarget = JumpFn->lookupByTarget(N);
        for (const TableCell &SourceValTargetValAndFunction :
             LookupByTarget.cellSet()) {
          d_t DPrime = SourceValTargetValAndFunction.getRowKey();
          d_t D = SourceValTargetValAndFunction.getColumnKey();
          EdgeFunctionPtrType FPrime = SourceValTargetValAndFunction.getValue();
          l_t TargetVal = val(SP, DPrime);
          setVal(N, D,
                 IDEProblem.join(val(N, D),
                                 FPrime->computeTarget(std::move(TargetVal))));
          INC_COUNTER("Value Computation", 1, PAMM_SEVERITY_LEVEL::Full);
        }
      }
    }
  }

  virtual void saveEdges(n_t SourceNode, n_t SinkStmt, d_t SourceVal,
                         const container_type &DestVals, bool InterP) {
    if (!SolverConfig.recordEdges()) {
      return;
    }
    Table<n_t, n_t, std::map<d_t, container_type>> &TgtMap =
        (InterP) ? ComputedInterPathEdges : ComputedIntraPathEdges;
    TgtMap.get(SourceNode, SinkStmt)[SourceVal].insert(DestVals.begin(),
                                                       DestVals.end());
  }

  /// Computes the final values for edge functions.
  void computeValues() {
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "Start computing values");
    // Phase II(i)
    std::map<n_t, std::set<d_t>> AllSeeds(InitialSeeds);
    for (n_t UnbalancedRetSite : UnbalancedRetSites) {
      if (AllSeeds[UnbalancedRetSite].empty()) {
        AllSeeds.emplace(UnbalancedRetSite, std::set<d_t>({ZeroValue}));
      }
    }
    // do processing
    for (const auto &Seed : AllSeeds) {
      n_t StartPoint = Seed.first;
      for (d_t Val : Seed.second) {
        // initialize the initial seeds with the top element as we have no
        // information at the beginning of the value computation problem
        setVal(StartPoint, Val, IDEProblem.topElement());
        std::pair<n_t, d_t> SuperGraphNode(StartPoint, Val);
        valuePropagationTask(SuperGraphNode);
      }
    }
    // Phase II(ii)
    // we create an array of all nodes and then dispatch fractions of this
    // array to multiple threads
    const std::set<n_t> AllNonCallStartNodes = ICF->allNonCallStartNodes();
    valueComputationTask(
        {AllNonCallStartNodes.begin(), AllNonCallStartNodes.end()});
  }

  /// Schedules the processing of initial seeds, initiating the analysis.
  /// Clients should only call this methods if performing synchronization on
  /// their own. Normally, solve() should be called instead.
  void submitInitialSeeds() {
    PAMM_GET_INSTANCE;
    for (const auto &[StartPoint, Facts] : InitialSeeds) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "Start point: " << IDEProblem.NtoString(StartPoint));
      for (const auto &Fact : Facts) {
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "\tFact: " << IDEProblem.DtoString(Fact);
                      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        if (!IDEProblem.isZeroValue(Fact)) {
          INC_COUNTER("Gen facts", 1, PAMM_SEVERITY_LEVEL::Core);
        }
        propagate(ZeroValue, StartPoint, Fact, EdgeIdentity<l_t>::getInstance(),
                  nullptr, false);
      }
      JumpFn->addFunction(ZeroValue, StartPoint, ZeroValue,
                          EdgeIdentity<l_t>::getInstance());
    }
  }

  /// Lines 21-32 of the algorithm.
  ///
  /// Stores callee-side summaries.
  /// Also, at the side of the caller, propagates intra-procedural flows to
  /// return sites using those newly computed summaries.
  ///
  /// @param edge an edge whose target node resembles a method exit
  ///
  virtual void processExit(const PathEdge<n_t, d_t> Edge) { //NOLINT
    PAMM_GET_INSTANCE;
    INC_COUNTER("Process Exit", 1, PAMM_SEVERITY_LEVEL::Full);
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Process exit at target: "
                  << IDEProblem.NtoString(Edge.getTarget()));
    n_t N = Edge.getTarget(); // an exit node; line 21...
    EdgeFunctionPtrType F = jumpFunction(Edge);
    f_t FunctionThatNeedsSummary = ICF->getFunctionOf(N);
    d_t D1 = Edge.factAtSource();
    d_t D2 = Edge.factAtTarget();
    // for each of the method's start points, determine incoming calls
    const std::set<n_t> StartPointsOf =
        ICF->getStartPointsOf(FunctionThatNeedsSummary);
    std::map<n_t, container_type> Inc;
    for (n_t SP : StartPointsOf) {
      // line 21.1 of Naeem/Lhotak/Rodriguez
      // register end-summary
      addEndSummary(SP, D1, N, D2, F);
      for (auto Entry : incoming(D1, SP)) {
        Inc[Entry.first] = Container{Entry.second};
      }
    }
    printEndSummaryTab();
    printIncomingTab();
    // for each incoming call edge already processed
    //(see processCall(..))
    for (auto Entry : Inc) {
      // line 22
      n_t C = Entry.first;
      // for each return site
      for (n_t RetSiteC : ICF->getReturnSitesOfCallAt(C)) {
        // compute return-flow function
        FlowFunctionPtrType RetFunction =
            CachedFlowEdgeFunctions.getRetFlowFunction(
                C, FunctionThatNeedsSummary, N, RetSiteC);
        INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        // for each incoming-call value
        for (d_t D4 : Entry.second) {
          const container_type Targets =
              computeReturnFlowFunction(RetFunction, D1, D2, C, Entry.second);
          ADD_TO_HISTOGRAM("Data-flow facts", targets.size(), 1,
                           PAMM_SEVERITY_LEVEL::Full);
          saveEdges(N, RetSiteC, D2, Targets, true);
          // for each target value at the return site
          // line 23
          for (d_t D5 : Targets) {
            // compute composed function
            // get call edge function
            EdgeFunctionPtrType F4 =
                CachedFlowEdgeFunctions.getCallEdgeFunction(
                    C, D4, ICF->getFunctionOf(N), D1);
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Queried Call Edge Function: " << F4->str());
            // get return edge function
            EdgeFunctionPtrType F5 =
                CachedFlowEdgeFunctions.getReturnEdgeFunction(
                    C, ICF->getFunctionOf(N), N, D2, RetSiteC, D5);
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Queried Return Edge Function: " << F5->str());
            if (SolverConfig.emitESG()) {
              for (auto SP : ICF->getStartPointsOf(ICF->getFunctionOf(N))) {
                IntermediateEdgeFunctions[std::make_tuple(C, D4, SP, D1)]
                    .push_back(F4);
              }
              IntermediateEdgeFunctions[std::make_tuple(N, D2, RetSiteC, D5)]
                  .push_back(F5);
            }
            INC_COUNTER("EF Queries", 2, PAMM_SEVERITY_LEVEL::Full);
            // compose call function * function * return function
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                              << "Compose: " << F5->str() << " * " << F->str()
                              << " * " << F4->str();
                          BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "         (return * function * call)");
            EdgeFunctionPtrType FPrime = F4->composeWith(F)->composeWith(F5);
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                              << "       = " << FPrime->str();
                          BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
            // for each jump function coming into the call, propagate to
            // return site using the composed function
            auto RevLookupResult = JumpFn->reverseLookup(C, D4);
            if (RevLookupResult) {
              for (auto ValAndFunc : RevLookupResult->get()) {
                EdgeFunctionPtrType F3 = ValAndFunc.second;
                if (!F3->equal_to(AllTop)) {
                  d_t D3 = ValAndFunc.first;
                  d_t D5RestoredCtx = restoreContextOnReturnedFact(C, D4, D5);
                  LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                                    << "Compose: " << FPrime->str() << " * "
                                    << F3->str();
                                BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
                  propagate(D3, RetSiteC, D5RestoredCtx,
                            F3->composeWitD5RestoredCtxfalse);
                }
              }
            }
          }
        }
      }
    }
    // handling for unbalanced problems where we return out of a method with a
    // fact for which we have no incoming flow.
    // note: we propagate that way only values that originate from ZERO, as
    // conditionally generated values should only
    // be propagated into callers that have an incoming edge for this
    // condition
    if (SolverConfig.followReturnsPastSeeds() && Inc.empty() &&
        IDEProblem.isZeroValue(D1)) {
      const std::set<n_t> Callers = ICF->getCallersOf(FunctionThatNeedsSummary);
      for (n_t C : Callers) {
        for (n_t RetSiteC : ICF->getReturnSitesOfCallAt(C)) {
          FlowFunctionPtrType RetFunction =
              CachedFlowEdgeFunctions.getRetFlowFunction(
                  C, FunctionThatNeedsSummary, N, RetSiteC);
          INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
          const container_type Targets = computeReturnFlowFunction(
              RetFunction, D1, D2, C, Container{ZeroValue});
          ADD_TO_HISTOGRAM("Data-flow facts", targets.size(), 1,
                           PAMM_SEVERITY_LEVEL::Full);
          saveEdges(N, RetSiteC, D2, Targets, true);
          for (d_t D5 : Targets) {
            EdgeFunctionPtrType F5 =
                CachedFlowEdgeFunctions.getReturnEdgeFunction(
                    C, ICF->getFunctionOf(N), N, D2, RetSiteC, D5);
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Queried Return Edge Function: " << F5->str());
            if (SolverConfig.emitESG()) {
              IntermediateEdgeFunctions[std::make_tuple(N, D2, RetSiteC, D5)]
                  .push_back(F5);
            }
            INC_COUNTER("EF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                              << "Compose: " << F5->str() << " * " << F->str();
                          BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
            propagteUnbalancedReturnFlow(RetSiteC, D5, F->composeWith(F5), C);
            // register for value processing (2nd IDE phase)
            UnbalancedRetSites.insert(RetSiteC);
          }
        }
      }
      // in cases where there are no callers, the return statement would
      // normally not be processed at all; this might be undesirable if
      // the flow function has a side effect such as registering a taint;
      // instead we thus call the return flow function will a null caller
      if (Callers.empty()) {
        FlowFunctionPtrType RetFunction =
            CachedFlowEdgeFunctions.getRetFlowFunction(
                nullptr, FunctionThatNeedsSummary, N, nullptr);
        INC_COUNTER("FF Queries", 1, PAMM_SEVERITY_LEVEL::Full);
        RetFunction->computeTargets(D2);
      }
    }
  }

  void propagteUnbalancedReturnFlow(n_t RetSiteC, d_t TargetVal,
                                    EdgeFunctionPtrType EdgeFunction,
                                    n_t RelatedCallSite) {
    propagate(ZeroValue, RetSiteC, TargetVal, std::move(EdgeFunction),
              RelatedCallSite, true);
  }

  /// This method will be called for each incoming edge and can be used to
  /// transfer knowledge from the calling edge to the returning edge, without
  /// affecting the summary edges at the callee.
  /// @param callSite
  ///
  /// @param d4
  ///            Fact stored with the incoming edge, i.e., present at the
  ///            caller side
  /// @param d5
  ///            Fact that originally should be propagated to the caller.
  /// @return Fact that will be propagated to the caller.
  ///
  d_t restoreContextOnReturnedFact(n_t  /*CallSite*/, d_t  /*D4*/, d_t D5) {
    // TODO support LinkedNode and JoinHandlingNode
    //		if (d5 instanceof LinkedNode) {
    //			((LinkedNode<D>) d5).setCallingContext(d4);
    //		}
    //		if(d5 instanceof JoinHandlingNode) {
    //			((JoinHandlingNode<D>)
    // d5).setCallingContext(d4);
    //		}
    return D5;
  }

  /// Computes the normal flow function for the given set of start and end
  /// abstractions-
  /// @param flowFunction The normal flow function to compute
  /// @param d1 The abstraction at the method's start node
  /// @param d2 The abstraction at the current node
  /// @return The set of abstractions at the successor node
  ///
  container_type
  computeNormalFlowFunction(const FlowFunctionPtrType &FlowFunction, d_t  /*D1*/,
                            d_t D2) {
    return FlowFunction->computeTargets(D2);
  }

  container_type
  computeSummaryFlowFunction(const FlowFunctionPtrType &SummaryFlowFunction,
                             d_t  /*D1*/, d_t D2) {
    return SummaryFlowFunction->computeTargets(D2);
  }

  /// Computes the call flow function for the given call-site abstraction
  /// @param callFlowFunction The call flow function to compute
  /// @param d1 The abstraction at the current method's start node.
  /// @param d2 The abstraction at the call site
  /// @return The set of caller-side abstractions at the callee's start node
  ///
  container_type
  computeCallFlowFunction(const FlowFunctionPtrType &CallFlowFunction, d_t  /*D1*/,
                          d_t D2) {
    return CallFlowFunction->computeTargets(D2);
  }

  /// Computes the call-to-return flow function for the given call-site
  /// abstraction
  /// @param callToReturnFlowFunction The call-to-return flow function to
  /// compute
  /// @param d1 The abstraction at the current method's start node.
  /// @param d2 The abstraction at the call site
  /// @return The set of caller-side abstractions at the return site
  ///
  container_type computeCallToReturnFlowFunction(
      const FlowFunctionPtrType &CallToReturnFlowFunction, d_t  /*D1*/, d_t D2) {
    return CallToReturnFlowFunction->computeTargets(D2);
  }

  /// Computes the return flow function for the given set of caller-side
  /// abstractions.
  /// @param retFunction The return flow function to compute
  /// @param d1 The abstraction at the beginning of the callee
  /// @param d2 The abstraction at the exit node in the callee
  /// @param callSite The call site
  /// @param callerSideDs The abstractions at the call site
  /// @return The set of caller-side abstractions at the return site
  ///
  container_type
  computeReturnFlowFunction(const FlowFunctionPtrType &RetFunction, d_t  /*D1*/,
                            d_t D2, n_t  /*CallSite*/,
                            const Container & /*callerSideDs*/) {
    return RetFunction->computeTargets(D2);
  }

  /// Propagates the flow further down the exploded super graph, merging any
  /// edge function that might already have been computed for targetVal at
  /// target.
  ///
  /// @param sourceVal the source value of the propagated summary edge
  /// @param target the target statement
  /// @param targetVal the target value at the target statement
  /// @param f the new edge function computed from (s0,sourceVal) to
  /// (target,targetVal)
  /// @param relatedCallSite for call and return flows the related call
  /// statement, nullptr otherwise (this value is not used within this
  /// implementation but may be useful for subclasses of IDESolver)
  /// @param isUnbalancedReturn true if this edge is propagating an
  /// unbalanced return (this value is not used within this implementation
  /// but may be useful for subclasses of {@link IDESolver})
  ///
  void
  propagate(d_t SourceVal, n_t Target, d_t TargetVal, //NOLINT
            const EdgeFunctionPtrType &F,
            /* deliberately exposed to clients */ n_t  /*RelatedCallSite*/,
            /* deliberately exposed to clients */ bool  /*IsUnbalancedReturn*/) {
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "Propagate flow";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Source value  : " << IDEProblem.DtoString(SourceVal);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Target        : " << IDEProblem.NtoString(Target);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Target value  : " << IDEProblem.DtoString(TargetVal);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Edge function : " << F.get()->str()
                  << " (result of previous compose)";
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');

    EdgeFunctionPtrType JumpFnE = [&]() {
      const auto RevLookupResult = JumpFn->reverseLookup(Target, TargetVal);
      if (RevLookupResult) {
        const auto &JumpFnContainer = RevLookupResult->get();
        const auto Find = std::find_if(
            JumpFnContainer.begin(), JumpFnContainer.end(),
            [SourceVal](auto &KVpair) { return KVpair.first == SourceVal; });
        if (Find != JumpFnContainer.end()) {
          return Find->second;
        }
      }
      // jump function is initialized to all-top if no entry was found
      return AllTop;
    }();
    EdgeFunctionPtrType FPrime = JumpFnE->joinWith(F);
    bool NewFunction = !(FPrime->equal_to(JumpFnE));

    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Join: " << JumpFnE->str() << " & " << F.get()->str()
                      << (JumpFnE->equal_to(F) ? " (EF's are equal)" : " ");
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "    = " << FPrime->str()
                  << (NewFunction ? " (new jump func)" : " ");
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
    if (NewFunction) {
      JumpFn->addFunction(SourceVal, Target, TargetVal, FPrime);
      const PathEdge<n_t, d_t> Edge(SourceVal, Target, TargetVal);
      PathEdgeCount++;
      pathEdgeProcessingTask(Edge);

      LOG_IF_ENABLE(if (!IDEProblem.isZeroValue(TargetVal)) {
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "EDGE: <F: " << Target->getFunction()->getName().str()
            << ", D: " << IDEProblem.DtoString(SourceVal) << '>';
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << " ---> <N: " << IDEProblem.NtoString(Target) << ',';
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "       D: " << IDEProblem.DtoString(TargetVal) << ',';
        BOOST_LOG_SEV(lg::get(), DEBUG) << "      EF: " << FPrime->str() << '>';
        BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
      });
    } else {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "PROPAGATE: No new function!");
    }
  }

  l_t joinValueAt(n_t  /*Unit*/, d_t  /*Fact*/, l_t Curr, l_t NewVal) {
    return IDEProblem.join(std::move(Curr), std::move(NewVal));
  }

  std::set<typename Table<n_t, d_t, EdgeFunctionPtrType>::Cell>
  endSummary(n_t SP, d_t D3) {
    if constexpr (PammCurrSevLevel >= PammSeverityLevel::Core) {
      auto Key = std::make_pair(SP, D3);
      auto FindNd = FSummaryReuse.find(Key);
      if (FindNd == FSummaryReuse.end()) {
        FSummaryReuse.emplace(Key, 0);
      } else {
        FSummaryReuse[Key] += 1;
      }
    }
    return Endsummarytab.get(SP, D3).cellSet();
  }

  std::map<n_t, container_type> incoming(d_t D1, n_t SP) {
    return Incomingtab.get(SP, D1);
  }

  void addIncoming(n_t SP, d_t D3, n_t N, d_t D2) {
    Incomingtab.get(SP, D3)[N].insert(D2);
  }

  void printIncomingTab() const { //NOLINT
#ifdef DYNAMIC_LOG
    if (boost::log::core::get()->get_logging_enabled()) {
      BOOST_LOG_SEV(lg::get(), DEBUG) << "Start of incomingtab entry";
      for (auto Cell : Incomingtab.cellSet()) {
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "sP: " << IDEProblem.NtoString(Cell.getRowKey());
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "d3: " << IDEProblem.DtoString(Cell.getColumnKey());
        for (auto Entry : Cell.getValue()) {
          BOOST_LOG_SEV(lg::get(), DEBUG)
              << "  n: " << IDEProblem.NtoString(Entry.first);
          for (auto Fact : Entry.second) {
            BOOST_LOG_SEV(lg::get(), DEBUG)
                << "  d2: " << IDEProblem.DtoString(Fact);
          }
        }
        BOOST_LOG_SEV(lg::get(), DEBUG) << "---------------";
      }
      BOOST_LOG_SEV(lg::get(), DEBUG) << "End of incomingtab entry";
      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
    }
#endif
  }

  void printEndSummaryTab() const {//NOLINT
#ifdef DYNAMIC_LOG 
    if (boost::log::core::get()->get_logging_enabled()) {
      BOOST_LOG_SEV(lg::get(), DEBUG) << "Start of endsummarytab entry";
      for (auto Cell : Endsummarytab.cellVec()) {
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "sP: " << IDEProblem.NtoString(Cell.getRowKey());
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "d1: " << IDEProblem.DtoString(Cell.getColumnKey());
        for (auto InnerCell : Cell.getValue().cellVec()) {
          BOOST_LOG_SEV(lg::get(), DEBUG)
              << "  eP: " << IDEProblem.NtoString(InnerCell.getRowKey());
          BOOST_LOG_SEV(lg::get(), DEBUG)
              << "  d2: " << IDEProblem.DtoString(InnerCell.getColumnKey());
          BOOST_LOG_SEV(lg::get(), DEBUG)
              << "  EF: " << InnerCell.getValue()->str();
          BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
        }
        BOOST_LOG_SEV(lg::get(), DEBUG) << "---------------";
        BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
      }
      BOOST_LOG_SEV(lg::get(), DEBUG) << "End of endsummarytab entry";
      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ';
    }
#endif
  }

  void printComputedPathEdges() {
    std::cout << "\n**********************************************************";
    std::cout << "\n*          Computed intra-procedural path egdes          *";
    std::cout
        << "\n**********************************************************\n";

    // Sort intra-procedural path edges
    auto Cells = ComputedIntraPathEdges.cellVec();
    StmtLess Stmtless(ICF);
    sort(Cells.begin(), Cells.end(), [&Stmtless](auto A, auto B) {
      return Stmtless(A.getRowKey(), B.getRowKey());
    });
    for (auto Cell : Cells) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      std::string N2Label = IDEProblem.NtoString(Edge.second);
      std::cout << "\nN1: " << IDEProblem.NtoString(Edge.first) << '\n'
                << "N2: " << N2Label << "\n----"
                << std::string(N2Label.size(), '-') << '\n';
      for (auto D1ToD2Set : Cell.getValue()) {
        auto D1Fact = D1ToD2Set.first;
        std::cout << "D1: " << IDEProblem.DtoString(D1Fact) << '\n';
        for (auto D2Fact : D1ToD2Set.second) {
          std::cout << "\tD2: " << IDEProblem.DtoString(D2Fact) << '\n';
        }
        std::cout << '\n';
      }
    }

    std::cout << "\n**********************************************************";
    std::cout << "\n*          Computed inter-procedural path edges          *";
    std::cout
        << "\n**********************************************************\n";

    // Sort intra-procedural path edges
    Cells = ComputedInterPathEdges.cellVec();
    sort(Cells.begin(), Cells.end(), [&Stmtless](auto A, auto B) {
      return Stmtless(A.getRowKey(), B.getRowKey());
    });
    for (auto Cell : Cells) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      std::string N2Label = IDEProblem.NtoString(Edge.second);
      std::cout << "\nN1: " << IDEProblem.NtoString(Edge.first) << '\n'
                << "N2: " << N2Label << "\n----"
                << std::string(N2Label.size(), '-') << '\n';
      for (auto D1ToD2Set : Cell.getValue()) {
        auto D1Fact = D1ToD2Set.first;
        std::cout << "D1: " << IDEProblem.DtoString(D1Fact) << '\n';
        for (auto D2Fact : D1ToD2Set.second) {
          std::cout << "\tD2: " << IDEProblem.DtoString(D2Fact) << '\n';
        }
        std::cout << '\n';
      }
    }
  }

  /// The invariant for computing the number of generated (#gen) and killed
  /// (#kill) facts:
  ///   (1) #Valid facts at the last statement <= #gen - #kill
  ///   (2) #gen >= #kill
  ///
  /// The total number of valid facts can be smaller than the difference of
  /// generated and killed facts, due to set semantics, i.e., a fact can be
  /// generated multiple times but appears only once.
  ///
  /// Zero value is not counted!
  ///
  /// @brief Computes and prints statistics of the analysis run, e.g. number of
  /// generated/killed facts, number of summary-reuses etc.
  ///
  void computeAndPrintStatistics() { //NOLINT
    PAMM_GET_INSTANCE;
    // Stores all valid facts at return site in caller context; return-site is
    // key
    std::unordered_map<n_t, std::set<d_t>> ValidInCallerContext;
    size_t NumGenFacts = 0;
    size_t NumKillFacts = 0;
    size_t NumIntraPathEdges = 0;
    size_t NumInterPathEdges = 0;
    // --- Intra-procedural Path Edges ---
    // d1 --> d2-Set
    // Case 1: d1 in d2-Set
    // Case 2: d1 not in d2-Set, i.e., d1 was killed. d2-Set could be empty.
    for (auto Cell : ComputedIntraPathEdges.cellSet()) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "N1: " << IDEProblem.NtoString(Edge.first);
                    BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "N2: " << IDEProblem.NtoString(Edge.second));
      for (auto &[D1, D2s] : Cell.getValue()) {
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "d1: " << IDEProblem.DtoString(D1));
        NumIntraPathEdges += D2s.size();
        // Case 1
        if (D2s.find(D1) != D2s.end()) {
          NumGenFacts += D2s.size() - 1;
        }
        // Case 2
        else {
          NumGenFacts += D2s.size();
          // We ignore the zero value
          if (!IDEProblem.isZeroValue(D1)) {
            NumKillFacts++;
          }
        }
        // Store all valid facts after call-to-return flow
        if (ICF->isCallSite(Edge.first)) {
          ValidInCallerContext[Edge.second].insert(D2s.begin(), D2s.end());
        }
        LOG_IF_ENABLE([this](const auto &D2s) {
          for (auto D2 : D2s) {
            BOOST_LOG_SEV(lg::get(), DEBUG)
                << "d2: " << IDEProblem.DtoString(D2);
          }
          BOOST_LOG_SEV(lg::get(), DEBUG) << "----";
        }(D2s));
      }
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
    }
    // Stores all pairs of (Startpoint, Fact) for which a summary was applied
    std::set<std::pair<n_t, d_t>> ProcessSummaryFacts;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "==============================================";
                  BOOST_LOG_SEV(lg::get(), DEBUG) << "INTER PATH EDGES");
    for (auto Cell : ComputedInterPathEdges.cellSet()) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "N1: " << IDEProblem.NtoString(Edge.first);
                    BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "N2: " << IDEProblem.NtoString(Edge.second));
      // --- Call-flow Path Edges ---
      // Case 1: d1 --> empty set
      //   Can be ignored, since killing a fact in the caller context will
      //   actually happen during  call-to-return.
      //
      // Case 2: d1 --> d2-Set
      //   Every fact d_i != ZeroValue in d2-set will be generated in the
      // callee context, thus counts as a new fact. Even if d1 is passed as it
      // is, it will count as a new fact. The reason for this is, that d1 can
      // be killed in the callee context, but still be valid in the caller
      // context.
      //
      // Special Case: Summary was applied for a particular call
      //   Process the summary's #gen and #kill.
      if (ICF->isCallSite(Edge.first)) {
        for (auto &[D1, D2s] : Cell.getValue()) {
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "d1: " << IDEProblem.DtoString(D1));
          NumInterPathEdges += D2s.size();
          for (auto D2 : D2s) {
            if (!IDEProblem.isZeroValue(D2)) {
              NumGenFacts++;
            }
            // Special case
            if (ProcessSummaryFacts.find(std::make_pair(Edge.second, D2)) !=
                ProcessSummaryFacts.end()) {
              std::multiset<d_t> SummaryDMultiSet =
                  Endsummarytab.get(Edge.second, D2).columnKeySet();
              // remove duplicates from multiset
              std::set<d_t> SummaryDSet(SummaryDMultiSet.begin(),
                                        SummaryDMultiSet.end());
              // Process summary just as an intra-procedural edge
              if (SummaryDSet.find(D2) != SummaryDSet.end()) {
                NumGenFacts += SummaryDSet.size() - 1;
              } else {
                NumGenFacts += SummaryDSet.size();
                // We ignore the zero value
                if (!IDEProblem.isZeroValue(D1)) {
                  NumKillFacts++;
                }
              }
            } else {
              ProcessSummaryFacts.emplace(Edge.second, D2);
            }
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "d2: " << IDEProblem.DtoString(D2));
          }
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "----");
        }
      }
      // --- Return-flow Path Edges ---
      // Since every fact passed to the callee was counted as a new fact, we
      // have to count every fact propagated to the caller as a kill to
      // satisfy our invariant. Obviously, every fact not propagated to the
      // caller will count as a kill. If an actual new fact is propagated to
      // the caller, we have to increase the number of generated facts by one.
      // Zero value does not count towards generated/killed facts.
      if (ICF->isExitInst(Cell.getRowKey())) {
        for (auto &[D1, D2s] : Cell.getValue()) {
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "d1: " << IDEProblem.DtoString(D1));
          NumInterPathEdges += D2s.size();
          auto CallerFacts = ValidInCallerContext[Edge.second];
          for (auto D2 : D2s) {
            // d2 not valid in caller context
            if (CallerFacts.find(D2) == CallerFacts.end()) {
              NumGenFacts++;
            }
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "d2: " << IDEProblem.DtoString(D2));
          }
          if (!IDEProblem.isZeroValue(D1)) {
            NumKillFacts++;
          }
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "----");
        }
      }
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
    }
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "SUMMARY REUSE");
    std::size_t TotalSummaryReuse = 0;
    for (auto Entry : FSummaryReuse) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "N1: " << IDEProblem.NtoString(Entry.first.first);
                    BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "D1: " << IDEProblem.DtoString(Entry.first.second);
                    BOOST_LOG_SEV(lg::get(), DEBUG)
                    << "#Reuse: " << Entry.second);
      TotalSummaryReuse += Entry.second;
    }
    INC_COUNTER("Gen facts", NumGenFacts, PAMM_SEVERITY_LEVEL::Core);
    INC_COUNTER("Kill facts", NumKillFacts, PAMM_SEVERITY_LEVEL::Core);
    INC_COUNTER("Summary-reuse", TotalSummaryReuse, PAMM_SEVERITY_LEVEL::Core);
    INC_COUNTER("Intra Path Edges", NumIntraPathEdges,
                PAMM_SEVERITY_LEVEL::Core);
    INC_COUNTER("Inter Path Edges", NumInterPathEdges,
                PAMM_SEVERITY_LEVEL::Core);

    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                      << "----------------------------------------------";
                  BOOST_LOG_SEV(lg::get(), INFO) << "=== Solver Statistics ===";
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "#Facts generated : " << GET_COUNTER("Gen facts");
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "#Facts killed    : " << GET_COUNTER("Kill facts");
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "#Summary-reuse   : " << GET_COUNTER("Summary-reuse");
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "#Intra Path Edges: " << GET_COUNTER("Intra Path Edges");
                  BOOST_LOG_SEV(lg::get(), INFO)
                  << "#Inter Path Edges: " << GET_COUNTER("Inter Path Edges"));
    if constexpr (PammCurrSevLevel >= PammSeverityLevel::Full) {
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
              << "Flow function query count: " << GET_COUNTER("FF Queries");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Edge function query count: " << GET_COUNTER("EF Queries");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Data-flow value propagation count: "
          << GET_COUNTER("Value Propagation");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Data-flow value computation count: "
          << GET_COUNTER("Value Computation");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Special flow function usage count: "
          << GET_COUNTER("SpecialSummary-FF Application");
          BOOST_LOG_SEV(lg::get(), INFO) << "Jump function construciton count: "
                                         << GET_COUNTER("JumpFn Construction");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Phase I duration: " << PRINT_TIMER("DFA Phase I");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Phase II duration: " << PRINT_TIMER("DFA Phase II");
          BOOST_LOG_SEV(lg::get(), INFO)
          << "----------------------------------------------");
      CachedFlowEdgeFunctions.print();
    }
  }

public:
  void enableESGAsDot() { SolverConfig.setEmitESG(); }

  void
  emitESGAsDot(std::ostream &OS = std::cout, //NOLINT
               std::string DotConfigDir = PhasarConfig::phasarDirectory()) {
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Emit Exploded super-graph (ESG) as DOT graph";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Process intra-procedural path egdes";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "=============================================");
    DOTGraph<d_t> G;
    DOTConfig::importDOTConfig(std::move(DotConfigDir));
    DOTFunctionSubGraph *FG = nullptr;

    // Sort intra-procedural path edges
    auto Cells = ComputedIntraPathEdges.cellVec();
    StmtLess Stmtless(ICF);
    sort(Cells.begin(), Cells.end(), [&Stmtless](auto A, auto B) {
      return Stmtless(A.getRowKey(), B.getRowKey());
    });
    for (auto Cell : Cells) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      std::string N1Label = IDEProblem.NtoString(Edge.first);
      std::string N2Label = IDEProblem.NtoString(Edge.second);
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "N1: " << N1Label;
                    BOOST_LOG_SEV(lg::get(), DEBUG) << "N2: " << N2Label);
      std::string N1StmtId = ICF->getStatementId(Edge.first);
      std::string N2StmtId = ICF->getStatementId(Edge.second);
      std::string FnName = ICF->getFunctionOf(Edge.first)->getName().str();
      // Get or create function subgraph
      if (!FG || FG->Id != FnName) {
        FG = &G.functions[FnName];
        FG->Id = FnName;
      }

      // Create control flow nodes
      DOTNode N1(FnName, N1Label, N1StmtId);
      DOTNode N2(FnName, N2Label, N2StmtId);
      // Add control flow node(s) to function subgraph
      FG->Stmts.insert(N1);
      if (ICF->isExitInst(Edge.second)) {
        FG->Stmts.insert(N2);
      }

      // Set control flow edge
      FG->IntraCfEdges.emplace(N1, N2);

      DOTFactSubGraph *D1Fsg = nullptr;
      unsigned D1FactId = 0;
      unsigned D2FactId = 0;
      for (auto D1ToD2Set : Cell.getValue()) {
        auto D1Fact = D1ToD2Set.first;
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "d1: " << IDEProblem.DtoString(D1Fact));

        DOTNode D1;
        if (IDEProblem.isZeroValue(D1Fact)) {
          D1 = {FnName, "Λ", N1StmtId, 0, false, true};
          D1FactId = 0;
        } else {
          // Get the fact-ID
          D1FactId = G.getFactID(D1Fact);
          std::string D1Label = IDEProblem.DtoString(D1Fact);

          // Get or create the fact subgraph
          D1Fsg = FG->getOrCreateFactSG(D1FactId, D1Label);

          // Insert D1 to fact subgraph
          D1 = {FnName, D1Label, N1StmtId, D1FactId, false, true};
          D1Fsg->Nodes.insert(std::make_pair(N1StmtId, D1));
        }

        DOTFactSubGraph *D2Fsg = nullptr;
        for (auto D2Fact : D1ToD2Set.second) {
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "d2: " << IDEProblem.DtoString(D2Fact));
          // We do not need to generate any intra-procedural nodes and edges
          // for the zero value since they will be auto-generated
          if (!IDEProblem.isZeroValue(D2Fact)) {
            // Get the fact-ID
            D2FactId = G.getFactID(D2Fact);
            std::string D2Label = IDEProblem.DtoString(D2Fact);
            DOTNode D2 = {FnName, D2Label, N2StmtId, D2FactId, false, true};
            std::string EFLabel;
            auto EFVec = IntermediateEdgeFunctions[std::make_tuple(
                Edge.first, D1Fact, Edge.second, D2Fact)];
            for (auto EF : EFVec) {
              EFLabel += EF->str() + ", ";
            }
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "EF LABEL: " << EFLabel);
            if (D1FactId == D2FactId && !IDEProblem.isZeroValue(D1Fact)) {
              assert(D1Fsg && "D1_FSG was nullptr but should be valid.");
              D1Fsg->Nodes.insert(std::make_pair(N2StmtId, D2));
              D1Fsg->Edges.emplace(D1, D2, true, EFLabel);
            } else {
              // Get or create the fact subgraph
              D2Fsg = FG->getOrCreateFactSG(D2FactId, D2Label);

              D2Fsg->Nodes.insert(std::make_pair(N2StmtId, D2));
              FG->CrossFactEdges.emplace(D1, D2, true, EFLabel);
            }
          }
        }
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "----------");
      }
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
    }

    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "=============================================";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "Process inter-procedural path edges";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "=============================================");
    Cells = ComputedInterPathEdges.cellVec();
    sort(Cells.begin(), Cells.end(), [&Stmtless](auto A, auto B) {
      return Stmtless(A.getRowKey(), B.getRowKey());
    });
    for (auto Cell : Cells) {
      auto Edge = std::make_pair(Cell.getRowKey(), Cell.getColumnKey());
      std::string N1Label = IDEProblem.NtoString(Edge.first);
      std::string N2Label = IDEProblem.NtoString(Edge.second);
      std::string FNameOfN1 = ICF->getFunctionOf(Edge.first)->getName().str();
      std::string FNameOfN2 = ICF->getFunctionOf(Edge.second)->getName().str();
      std::string N1StmtId = ICF->getStatementId(Edge.first);
      std::string N2StmtId = ICF->getStatementId(Edge.second);
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "N1: " << N1Label;
                    BOOST_LOG_SEV(lg::get(), DEBUG) << "N2: " << N2Label);

      // Add inter-procedural control flow edge
      DOTNode N1(FNameOfN1, N1Label, N1StmtId);
      DOTNode N2(FNameOfN2, N2Label, N2StmtId);

      // Handle recursion control flow as intra-procedural control flow
      // since those eges never leave the function subgraph
      FG = nullptr;
      if (FNameOfN1 == FNameOfN2) {
        // This function subgraph is guaranteed to exist
        FG = &G.functions[FNameOfN1];
        FG->IntraCfEdges.emplace(N1, N2);
      } else {
        // Check the case where the callee is a single statement function,
        // thus does not contain intra-procedural path edges. We have to
        // generate the function sub graph here!
        if (!G.functions.count(FNameOfN1)) {
          FG = &G.functions[FNameOfN1];
          FG->Id = FNameOfN1;
          FG->Stmts.insert(N1);
        } else if (!G.functions.count(FNameOfN2)) {
          FG = &G.functions[FNameOfN2];
          FG->Id = FNameOfN2;
          FG->Stmts.insert(N2);
        }
        G.interCFEdges.emplace(N1, N2);
      }

      // Create D1 and D2, if D1 == D2 == lambda then add Edge(D1, D2) to
      // interLambdaEges otherwise add Edge(D1, D2) to interFactEdges
      unsigned D1FactId = 0;
      unsigned D2FactId = 0;
      for (auto D1ToD2Set : Cell.getValue()) {
        auto D1Fact = D1ToD2Set.first;
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "d1: " << IDEProblem.DtoString(D1Fact));
        DOTNode D1;
        if (IDEProblem.isZeroValue(D1Fact)) {
          D1 = {FNameOfN1, "Λ", N1StmtId, 0, false, true};
        } else {
          // Get the fact-ID
          D1FactId = G.getFactID(D1Fact);
          std::string D1Label = IDEProblem.DtoString(D1Fact);
          D1 = {FNameOfN1, D1Label, N1StmtId, D1FactId, false, true};
          // FG should already exist even for single statement functions
          if (!G.containsFactSG(FNameOfN1, D1FactId)) {
            FG = &G.functions[FNameOfN1];
            auto *D1Fsg = FG->getOrCreateFactSG(D1FactId, D1Label);
            D1Fsg->Nodes.insert(std::make_pair(N1StmtId, D1));
          }
        }

        auto D2Set = D1ToD2Set.second;
        for (auto D2Fact : D2Set) {
          LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "d2: " << IDEProblem.DtoString(D2Fact));
          DOTNode D2;
          if (IDEProblem.isZeroValue(D2Fact)) {
            D2 = {FNameOfN2, "Λ", N2StmtId, 0, false, true};
          } else {
            // Get the fact-ID
            D2FactId = G.getFactID(D2Fact);
            std::string D2Label = IDEProblem.DtoString(D2Fact);
            D2 = {FNameOfN2, D2Label, N2StmtId, D2FactId, false, true};
            // FG should already exist even for single statement functions
            if (!G.containsFactSG(FNameOfN2, D2FactId)) {
              FG = &G.functions[FNameOfN2];
              auto *D2Fsg = FG->getOrCreateFactSG(D2FactId, D2Label);
              D2Fsg->Nodes.insert(std::make_pair(N2StmtId, D2));
            }
          }

          if (IDEProblem.isZeroValue(D1Fact) &&
              IDEProblem.isZeroValue(D2Fact)) {
            // Do not add lambda recursion edges as inter-procedural edges
            if (D1.FuncName != D2.FuncName) {
              G.interLambdaEdges.emplace(D1, D2, true, "AllBottom", "BOT");
            }
          } else {
            // std::string EFLabel = EF ? EF->str() : " ";
            std::string EFLabel;
            auto EFVec = IntermediateEdgeFunctions[std::make_tuple(
                Edge.first, D1Fact, Edge.second, D2Fact)];
            for (auto EF : EFVec) {
              LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                            << "Partial EF Label: " << EF->str());
              EFLabel.append(EF->str() + ", ");
            }
            LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "EF LABEL: " << EFLabel);
            G.interFactEdges.emplace(D1, D2, true, EFLabel);
          }
        }
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << "----------");
      }
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG) << " ");
    }
    OS << G;
  }

  /// @brief: Allows less-than comparison based on the statement ID.
  struct StmtLess {
    const i_t *ICF;
    StringIdLess StrIdLess;
    StmtLess(const i_t *ICF) : ICF(ICF), StrIdLess(StringIdLess()) {}
    bool operator()(n_t Lhs, n_t Rhs) {
      return StrIdLess(ICF->getStatementId(Lhs), ICF->getStatementId(Rhs));
    }
  };
};

template <typename Problem>
IDESolver(Problem &) -> IDESolver<typename Problem::ProblemAnalysisDomain,
                                  typename Problem::container_type>;

template <typename Problem>
using IDESolver_P = IDESolver<typename Problem::ProblemAnalysisDomain,
                              typename Problem::container_type>;

} // namespace psr

#endif
