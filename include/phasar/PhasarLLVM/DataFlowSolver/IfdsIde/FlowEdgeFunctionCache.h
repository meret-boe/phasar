/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

#ifndef PHASAR_PHASARLLVM_IFDSIDE_FLOWEDGEFUNCTIONCACHE_H_
#define PHASAR_PHASARLLVM_IFDSIDE_FLOWEDGEFUNCTIONCACHE_H_

#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFact.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/EdgeFunctions.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/IDETabulationProblem.h"
#include "phasar/Utils/EquivalenceClassMap.h"
#include "phasar/Utils/Logger.h"
#include "phasar/Utils/PAMMMacros.h"

#include "llvm/ADT/DenseMap.h"

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>

namespace psr {
template <typename KeyT> class DefaultMapKeyCompressor {
public:
  using KeyType = KeyT;
  using CompressedType = KeyT;

  [[nodiscard]] inline CompressedType getCompressedID(KeyT Key) { return Key; }
};

template <typename... Ts> class MapKeyCompressorCombinator : public Ts... {
public:
  using Ts::getCompressedID...;
};

class LLVMMapKeyCompressor {
public:
  using KeyType = const llvm::Value *;
  using CompressedType = uint32_t;

  [[nodiscard]] inline CompressedType getCompressedID(KeyType Key) {
    auto Search = Map.find(Key);
    if (Search == Map.end()) {
      return Map.insert(std::make_pair(Key, Map.size() + 1)).first->getSecond();
    }
    return Search->getSecond();
  }

private:
  llvm::DenseMap<KeyType, CompressedType> Map{};
};

/**
 * This class caches flow and edge functions to avoid their reconstruction.
 * When a flow or edge function must be applied to multiple times, a cached
 * version is used if existend, otherwise a new one is created and inserted
 * into the cache.
 */
template <typename AnalysisDomainTy,
          typename Container = std::set<typename AnalysisDomainTy::d_t>>
class FlowEdgeFunctionCache {
  using IDEProblemType = IDETabulationProblem<AnalysisDomainTy, Container>;
  using FlowFunctionPtrType = typename IDEProblemType::FlowFunctionPtrType;
  using EdgeFunctionPtrType = typename IDEProblemType::EdgeFunctionPtrType;

  using n_t = typename AnalysisDomainTy::n_t;
  using d_t = typename AnalysisDomainTy::d_t;
  using f_t = typename AnalysisDomainTy::f_t;
  using t_t = typename AnalysisDomainTy::t_t;

  using DTKeyCompressorType = std::conditional_t<
      std::is_base_of_v<llvm::Value, std::remove_pointer_t<d_t>>,
      LLVMMapKeyCompressor, DefaultMapKeyCompressor<d_t>>;
  using NTKeyCompressorType = std::conditional_t<
      std::is_base_of_v<llvm::Value, std::remove_pointer_t<n_t>>,
      LLVMMapKeyCompressor, DefaultMapKeyCompressor<n_t>>;

  using MapKeyCompressorType = std::conditional_t<
      std::is_same_v<NTKeyCompressorType, DTKeyCompressorType>,
      NTKeyCompressorType,
      MapKeyCompressorCombinator<NTKeyCompressorType, DTKeyCompressorType>>;

private:
  MapKeyCompressorType KeyCompressor;

  using EdgeFuncInstKey = uint64_t;
  using EdgeFuncNodeKey = std::conditional_t<
      std::is_base_of_v<llvm::Value, std::remove_pointer_t<d_t>>, uint64_t,
      std::pair<d_t, d_t>>;
  using InnerEdgeFunctionMapType =
      EquivalenceClassMap<EdgeFuncNodeKey, EdgeFunctionPtrType>;

  IDETabulationProblem<AnalysisDomainTy, Container> &Problem;
  // Auto add zero
  bool AutoAddZero;
  d_t ZeroValue;

  struct NormalEdgeFlowData {
    NormalEdgeFlowData(FlowFunctionPtrType Val)
        : FlowFuncPtr(std::move(Val)), EdgeFunctionMap{} {}
    NormalEdgeFlowData(InnerEdgeFunctionMapType Map)
        : FlowFuncPtr(nullptr), EdgeFunctionMap{std::move(Map)} {}

    FlowFunctionPtrType FlowFuncPtr;
    InnerEdgeFunctionMapType EdgeFunctionMap;
  };

  // Caches for the flow/edge functions
  std::map<EdgeFuncInstKey, NormalEdgeFlowData> NormalFunctionCache;

  // Caches for the flow functions
  std::map<std::tuple<n_t, f_t>, FlowFunctionPtrType> CallFlowFunctionCache;
  std::map<std::tuple<n_t, f_t, n_t, n_t>, FlowFunctionPtrType>
      ReturnFlowFunctionCache;
  std::map<std::tuple<n_t, n_t, std::set<f_t>>, FlowFunctionPtrType>
      CallToRetFlowFunctionCache;
  // Caches for the edge functions
  std::map<std::tuple<n_t, d_t, f_t, d_t>, EdgeFunctionPtrType>
      CallEdgeFunctionCache;
  std::map<std::tuple<n_t, f_t, n_t, d_t, n_t, d_t>, EdgeFunctionPtrType>
      ReturnEdgeFunctionCache;
  std::map<EdgeFuncInstKey, InnerEdgeFunctionMapType>
      CallToRetEdgeFunctionCache;
  std::map<std::tuple<n_t, d_t, n_t, d_t>, EdgeFunctionPtrType>
      SummaryEdgeFunctionCache;

public:
  // Ctor allows access to the IDEProblem in order to get access to flow and
  // edge function factory functions.
  FlowEdgeFunctionCache(
      IDETabulationProblem<AnalysisDomainTy, Container> &Problem)
      : Problem(Problem),
        AutoAddZero(Problem.getIFDSIDESolverConfig().autoAddZero()),
        ZeroValue(Problem.getZeroValue()) {
    PAMM_GET_INSTANCE;
    REG_COUNTER("Normal-FF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Normal-FF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the call flow functions
    REG_COUNTER("Call-FF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Call-FF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for return flow functions
    REG_COUNTER("Return-FF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Return-FF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the call to return flow functions
    REG_COUNTER("CallToRet-FF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("CallToRet-FF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the summary flow functions
    // REG_COUNTER("Summary-FF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    // REG_COUNTER("Summary-FF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the normal edge functions
    REG_COUNTER("Normal-EF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Normal-EF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the call edge functions
    REG_COUNTER("Call-EF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Call-EF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the return edge functions
    REG_COUNTER("Return-EF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Return-EF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the call to return edge functions
    REG_COUNTER("CallToRet-EF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("CallToRet-EF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
    // Counters for the summary edge functions
    REG_COUNTER("Summary-EF Construction", 0, PAMM_SEVERITY_LEVEL::Full);
    REG_COUNTER("Summary-EF Cache Hit", 0, PAMM_SEVERITY_LEVEL::Full);
  }

  ~FlowEdgeFunctionCache() = default;

  FlowEdgeFunctionCache(const FlowEdgeFunctionCache &FEFC) = default;
  FlowEdgeFunctionCache &operator=(const FlowEdgeFunctionCache &FEFC) = default;

  FlowEdgeFunctionCache(FlowEdgeFunctionCache &&FEFC) noexcept = default;
  FlowEdgeFunctionCache &
  operator=(FlowEdgeFunctionCache &&FEFC) noexcept = default;

  FlowFunctionPtrType getNormalFlowFunction(n_t Curr, n_t Succ) { //NOLINT
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Normal flow function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Curr Inst : " << Problem.NtoString(Curr);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Succ Inst : " << Problem.NtoString(Succ));
    auto Key = createEdgeFunctionInstKey(Curr, Succ);
    auto SearchNormalFlowFunction = NormalFunctionCache.find(Key);
    if (SearchNormalFlowFunction != NormalFunctionCache.end()) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      INC_COUNTER("Normal-FF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      if (SearchNormalFlowFunction->second.FlowFuncPtr != nullptr) {
        return SearchNormalFlowFunction->second.FlowFuncPtr;
      } 
      auto Ff =
            (AutoAddZero)
                ? std::make_shared<ZeroedFlowFunction<d_t, Container>>(
                      Problem.getNormalFlowFunction(Curr, Succ), ZeroValue)
                : Problem.getNormalFlowFunction(Curr, Succ);
        SearchNormalFlowFunction->second.FlowFuncPtr = Ff;
      return Ff;
     
    } 
      INC_COUNTER("Normal-FF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ff = (AutoAddZero)
                    ? std::make_shared<ZeroedFlowFunction<d_t, Container>>(
                          Problem.getNormalFlowFunction(Curr, Succ), ZeroValue)
                    : Problem.getNormalFlowFunction(Curr, Succ);
      NormalFunctionCache.insert(std::make_pair(Key, NormalEdgeFlowData(Ff)));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ff;
   
  }

  FlowFunctionPtrType getCallFlowFunction(n_t CallSite, f_t DestFun) { //NOLINT
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Call flow function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Call Stmt : " << Problem.NtoString(CallSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(F) Dest Fun : " << Problem.FtoString(DestFun));
    auto Key = std::tie(CallSite, DestFun);
    auto SearchCallFlowFunction = CallFlowFunctionCache.find(Key);
    if (SearchCallFlowFunction != CallFlowFunctionCache.end()) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      INC_COUNTER("Call-FF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      return SearchCallFlowFunction->second;
    } 
      INC_COUNTER("Call-FF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ff =
          (AutoAddZero)
              ? std::make_shared<ZeroedFlowFunction<d_t, Container>>(
                    Problem.getCallFlowFunction(CallSite, DestFun), ZeroValue)
              : Problem.getCallFlowFunction(CallSite, DestFun);
      CallFlowFunctionCache.insert(std::make_pair(Key, Ff));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ff;
   
  }

  FlowFunctionPtrType getRetFlowFunction(n_t CallSite, f_t CalleeFun, //NOLINT
                                         n_t ExitInst, n_t RetSite) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Return flow function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Call Site : " << Problem.NtoString(CallSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(F) Callee    : " << Problem.FtoString(CalleeFun);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Exit Stmt : " << Problem.NtoString(ExitInst);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Ret Site  : " << Problem.NtoString(RetSite));
    auto Key = std::tie(CallSite, CalleeFun, ExitInst, RetSite);
    auto SearchReturnFlowFunction = ReturnFlowFunctionCache.find(Key);
    if (SearchReturnFlowFunction != ReturnFlowFunctionCache.end()) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      INC_COUNTER("Return-FF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      return SearchReturnFlowFunction->second;
    } 
      INC_COUNTER("Return-FF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ff = (AutoAddZero)
                    ? std::make_shared<ZeroedFlowFunction<d_t, Container>>(
                          Problem.getRetFlowFunction(CallSite, CalleeFun,
                                                     ExitInst, RetSite),
                          ZeroValue)
                    : Problem.getRetFlowFunction(CallSite, CalleeFun, ExitInst,
                                                 RetSite);
      ReturnFlowFunctionCache.insert(std::make_pair(Key, Ff));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ff;
   
  }

  FlowFunctionPtrType getCallToRetFlowFunction(n_t CallSite, n_t RetSite, //NOLINT
                                               std::set<f_t> Callees) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "Call-to-Return flow function factory call";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(N) Call Site : " << Problem.NtoString(CallSite);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(N) Ret Site  : " << Problem.NtoString(RetSite);
        BOOST_LOG_SEV(lg::get(), DEBUG) << "(F) Callee's  : "; for (auto callee
                                                                    : Callees) {
          BOOST_LOG_SEV(lg::get(), DEBUG) << "  " << Problem.FtoString(callee);
        });
    auto Key = std::tie(CallSite, RetSite, Callees);
    auto SearchCallToRetFlowFunction = CallToRetFlowFunctionCache.find(Key);
    if (SearchCallToRetFlowFunction != CallToRetFlowFunctionCache.end()) {
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      INC_COUNTER("CallToRet-FF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      return SearchCallToRetFlowFunction->second;
    } 
      INC_COUNTER("CallToRet-FF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ff =
          (AutoAddZero)
              ? std::make_shared<ZeroedFlowFunction<d_t, Container>>(
                    Problem.getCallToRetFlowFunction(CallSite, RetSite,
                                                     Callees),
                    ZeroValue)
              : Problem.getCallToRetFlowFunction(CallSite, RetSite, Callees);
      CallToRetFlowFunctionCache.insert(std::make_pair(Key, Ff));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Flow function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ff;
   
  }

  FlowFunctionPtrType getSummaryFlowFunction(n_t CallSite, f_t DestFun) {
    // PAMM_GET_INSTANCE;
    // INC_COUNTER("Summary-FF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Summary flow function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Call Stmt : " << Problem.NtoString(CallSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(F) Dest Mthd : " << Problem.FtoString(DestFun);
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
    auto Ff = Problem.getSummaryFlowFunction(CallSite, DestFun);
    return Ff;
  }

  EdgeFunctionPtrType getNormalEdgeFunction(n_t Curr, d_t CurrNode, n_t Succ, //NOLINT
                                            d_t SuccNode) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Normal edge function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Curr Inst : " << Problem.NtoString(Curr);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Curr Node : " << Problem.DtoString(CurrNode);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Succ Inst : " << Problem.NtoString(Succ);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Succ Node : " << Problem.DtoString(SuccNode));

    EdgeFuncInstKey OuterMapKey = createEdgeFunctionInstKey(Curr, Succ);
    auto SearchInnerMap = NormalFunctionCache.find(OuterMapKey);
    if (SearchInnerMap != NormalFunctionCache.end()) {
      auto SearchEdgeFunc = SearchInnerMap->second.EdgeFunctionMap.find(
          createEdgeFunctionNodeKey(CurrNode, SuccNode));
      if (SearchEdgeFunc != SearchInnerMap->second.EdgeFunctionMap.end()) {
        INC_COUNTER("Normal-EF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Edge function fetched from cache";
                      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        return SearchEdgeFunc->second;
      }
      INC_COUNTER("Normal-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ef = Problem.getNormalEdgeFunction(Curr, CurrNode, Succ, SuccNode);

      SearchInnerMap->second.EdgeFunctionMap.insert(
          createEdgeFunctionNodeKey(CurrNode, SuccNode), Ef);

      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ef;
    }
    INC_COUNTER("Normal-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
    auto Ef = Problem.getNormalEdgeFunction(Curr, CurrNode, Succ, SuccNode);

    NormalFunctionCache.try_emplace(
        OuterMapKey, NormalEdgeFlowData(InnerEdgeFunctionMapType{std::make_pair(
                         createEdgeFunctionNodeKey(CurrNode, SuccNode), Ef)}));

    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Edge function constructed";
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
    return Ef;
  }

  EdgeFunctionPtrType getCallEdgeFunction(n_t CallSite, d_t SrcNode, //NOLINT
                                          f_t DestinationFunction,
                                          d_t DestNode) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG) << "Call edge function factory call";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(N) Call Stmt : " << Problem.NtoString(CallSite);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(D) Src Node  : " << Problem.DtoString(SrcNode);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(F) Dest Fun : " << Problem.FtoString(DestinationFunction);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(D) Dest Node : " << Problem.DtoString(DestNode));
    auto Key = std::tie(CallSite, SrcNode, DestinationFunction, DestNode);
    auto SearchCallEdgeFunction = CallEdgeFunctionCache.find(Key);
    if (SearchCallEdgeFunction != CallEdgeFunctionCache.end()) {
      INC_COUNTER("Call-EF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return SearchCallEdgeFunction->second;
    } 
      INC_COUNTER("Call-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ef = Problem.getCallEdgeFunction(CallSite, SrcNode,
                                            DestinationFunction, DestNode);
      CallEdgeFunctionCache.insert(std::make_pair(Key, Ef));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ef;
   
  }

  EdgeFunctionPtrType getReturnEdgeFunction(n_t CallSite, f_t CalleeFunction, //NOLINT
                                            n_t ExitInst, d_t ExitNode,
                                            n_t ReSite, d_t RetNode) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Return edge function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Call Site : " << Problem.NtoString(CallSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(F) Callee    : " << Problem.FtoString(CalleeFunction);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Exit Stmt : " << Problem.NtoString(ExitInst);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Exit Node : " << Problem.DtoString(ExitNode);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Ret Site  : " << Problem.NtoString(ReSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Ret Node  : " << Problem.DtoString(RetNode));
    auto Key =
        std::tie(CallSite, CalleeFunction, ExitInst, ExitNode, ReSite, RetNode);
    auto SearchReturnEdgeFunction = ReturnEdgeFunctionCache.find(Key);
    if (SearchReturnEdgeFunction != ReturnEdgeFunctionCache.end()) {
      INC_COUNTER("Return-EF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return SearchReturnEdgeFunction->second;
    } 
      INC_COUNTER("Return-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ef = Problem.getReturnEdgeFunction(
          CallSite, CalleeFunction, ExitInst, ExitNode, ReSite, RetNode);
      ReturnEdgeFunctionCache.insert(std::make_pair(Key, Ef));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ef;
   
  }

  EdgeFunctionPtrType getCallToRetEdgeFunction(n_t CallSite, d_t CallNode, //NOLINT
                                               n_t RetSite, d_t RetSiteNode,
                                               std::set<f_t> Callees) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(
        BOOST_LOG_SEV(lg::get(), DEBUG)
            << "Call-to-Return edge function factory call";
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(N) Call Site : " << Problem.NtoString(CallSite);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(D) Call Node : " << Problem.DtoString(CallNode);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(N) Ret Site  : " << Problem.NtoString(RetSite);
        BOOST_LOG_SEV(lg::get(), DEBUG)
        << "(D) Ret Node  : " << Problem.DtoString(RetSiteNode);
        BOOST_LOG_SEV(lg::get(), DEBUG) << "(F) Callee's  : "; for (auto callee
                                                                    : Callees) {
          BOOST_LOG_SEV(lg::get(), DEBUG) << "  " << Problem.FtoString(callee);
        });

    EdgeFuncInstKey OuterMapKey = createEdgeFunctionInstKey(CallSite, RetSite);
    auto SearchInnerMap = CallToRetEdgeFunctionCache.find(OuterMapKey);
    if (SearchInnerMap != CallToRetEdgeFunctionCache.end()) {
      auto SearchEdgeFunc = SearchInnerMap->second.find(
          createEdgeFunctionNodeKey(CallNode, RetSiteNode));
      if (SearchEdgeFunc != SearchInnerMap->second.end()) {
        INC_COUNTER("CTR-EF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
        LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                          << "Edge function fetched from cache";
                      BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
        return SearchEdgeFunc->second;
      }
      INC_COUNTER("CTR-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ef = Problem.getCallToRetEdgeFunction(CallSite, CallNode, RetSite,
                                                 RetSiteNode, Callees);

      SearchInnerMap->second.insert(
          createEdgeFunctionNodeKey(CallNode, RetSiteNode), Ef);

      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ef;
    }

    INC_COUNTER("CTR-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
    auto Ef = Problem.getCallToRetEdgeFunction(CallSite, CallNode, RetSite,
                                               RetSiteNode, Callees);

    CallToRetEdgeFunctionCache.emplace(
        OuterMapKey,
        InnerEdgeFunctionMapType{std::make_pair(
            createEdgeFunctionNodeKey(CallNode, RetSiteNode), Ef)});
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Edge function constructed";
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
    return Ef;
  }

  EdgeFunctionPtrType getSummaryEdgeFunction(n_t CallSite, d_t CallNode, //NOLINT
                                             n_t RetSite, d_t RetSiteNode) {
    PAMM_GET_INSTANCE;
    LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                      << "Summary edge function factory call";
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Call Site : " << Problem.NtoString(CallSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Call Node : " << Problem.DtoString(CallNode);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(N) Ret Site  : " << Problem.NtoString(RetSite);
                  BOOST_LOG_SEV(lg::get(), DEBUG)
                  << "(D) Ret Node  : " << Problem.DtoString(RetSiteNode);
                  BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
    auto Key = std::tie(CallSite, CallNode, RetSite, RetSiteNode);
    auto SearchSummaryEdgeFunction = SummaryEdgeFunctionCache.find(Key);
    if (SearchSummaryEdgeFunction != SummaryEdgeFunctionCache.end()) {
      INC_COUNTER("Summary-EF Cache Hit", 1, PAMM_SEVERITY_LEVEL::Full);
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function fetched from cache";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return SearchSummaryEdgeFunction->second;
    } 
      INC_COUNTER("Summary-EF Construction", 1, PAMM_SEVERITY_LEVEL::Full);
      auto Ef = Problem.getSummaryEdgeFunction(CallSite, CallNode, RetSite,
                                               RetSiteNode);
      SummaryEdgeFunctionCache.insert(std::make_pair(Key, Ef));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), DEBUG)
                        << "Edge function constructed";
                    BOOST_LOG_SEV(lg::get(), DEBUG) << ' ');
      return Ef;
   
  }

  void print() { //NOLINT
    if constexpr (PammCurrSevLevel >= PammSeverityLevel::Full) {
      PAMM_GET_INSTANCE;
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "=== Flow-Edge-Function Cache Statistics ===");
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Normal-flow function cache hits: "
                    << GET_COUNTER("Normal-FF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Normal-flow function constructions: "
                    << GET_COUNTER("Normal-FF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-flow function cache hits: "
                    << GET_COUNTER("Call-FF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-flow function constructions: "
                    << GET_COUNTER("Call-FF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Return-flow function cache hits: "
                    << GET_COUNTER("Return-FF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Return-flow function constructions: "
                    << GET_COUNTER("Return-FF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-to-Return-flow function cache hits: "
                    << GET_COUNTER("CallToRet-FF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-to-Return-flow function constructions: "
                    << GET_COUNTER("CallToRet-FF Construction"));
      // LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO) << "Summary-flow function
      // cache hits: "
      //                        << GET_COUNTER("Summary-FF Cache Hit"));
      // LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO) << "Summary-flow function
      // constructions: "
      //                         << GET_COUNTER("Summary-FF Construction"));
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Total flow function cache hits: "
          << GET_SUM_COUNT({"Normal-FF Cache Hit", "Call-FF Cache Hit",
                            "Return-FF Cache Hit", "CallToRet-FF Cache Hit"}));
      //"Summary-FF Cache Hit"});
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Total flow function constructions: "
          << GET_SUM_COUNT({"Normal-FF Construction", "Call-FF Construction",
                            "Return-FF Construction",
                            "CallToRet-FF Construction" /*,
                "Summary-FF Construction"*/}));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO) << ' ');
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Normal edge function cache hits: "
                    << GET_COUNTER("Normal-EF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Normal edge function constructions: "
                    << GET_COUNTER("Normal-EF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call edge function cache hits: "
                    << GET_COUNTER("Call-EF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call edge function constructions: "
                    << GET_COUNTER("Call-EF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Return edge function cache hits: "
                    << GET_COUNTER("Return-EF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Return edge function constructions: "
                    << GET_COUNTER("Return-EF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-to-Return edge function cache hits: "
                    << GET_COUNTER("CallToRet-EF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Call-to-Return edge function constructions: "
                    << GET_COUNTER("CallToRet-EF Construction"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Summary edge function cache hits: "
                    << GET_COUNTER("Summary-EF Cache Hit"));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "Summary edge function constructions: "
                    << GET_COUNTER("Summary-EF Construction"));
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Total edge function cache hits: "
          << GET_SUM_COUNT({"Normal-EF Cache Hit", "Call-EF Cache Hit",
                            "Return-EF Cache Hit", "CallToRet-EF Cache Hit",
                            "Summary-EF Cache Hit"}));
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Total edge function constructions: "
          << GET_SUM_COUNT({"Normal-EF Construction", "Call-EF Construction",
                            "Return-EF Construction",
                            "CallToRet-EF Construction",
                            "Summary-EF Construction"}));
      LOG_IF_ENABLE(BOOST_LOG_SEV(lg::get(), INFO)
                    << "----------------------------------------------");
    } else {
      LOG_IF_ENABLE(
          BOOST_LOG_SEV(lg::get(), INFO)
          << "Cache statistics only recorded on PAMM severity level: Full.");
    }
  }

private:
  inline EdgeFuncInstKey createEdgeFunctionInstKey(n_t N1, n_t N2) {
    uint64_t Val = 0;
    Val |= KeyCompressor.getCompressedID(N1);
    Val <<= 32;
    Val |= KeyCompressor.getCompressedID(N2);
    return Val;
  }

  inline EdgeFuncNodeKey createEdgeFunctionNodeKey(d_t D1, d_t D2) {
    if constexpr (std::is_base_of_v<llvm::Value, std::remove_pointer_t<d_t>>) {
      uint64_t Val = 0;
      Val |= KeyCompressor.getCompressedID(D1);
      Val <<= 32;
      Val |= KeyCompressor.getCompressedID(D2);
      return Val;
    } else {
      return std::make_pair(D1, D2);
    }
  }
};

} // namespace psr

#endif
