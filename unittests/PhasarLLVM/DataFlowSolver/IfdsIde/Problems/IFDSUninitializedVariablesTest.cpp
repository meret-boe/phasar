#include <memory>

#include "phasar/DB/ProjectIRDB.h"
#include "phasar/PhasarLLVM/ControlFlow/LLVMBasedICFG.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Problems/IFDSUninitializedVariables.h"
#include "phasar/PhasarLLVM/DataFlowSolver/IfdsIde/Solver/IFDSSolver.h"
#include "phasar/PhasarLLVM/Passes/ValueAnnotationPass.h"
#include "phasar/PhasarLLVM/Pointer/LLVMPointsToSet.h"
#include "phasar/PhasarLLVM/TypeHierarchy/LLVMTypeHierarchy.h"
#include "llvm/IR/Module.h"
#include "gtest/gtest.h"

#include "TestConfig.h"

using namespace std;
using namespace psr;

/* ============== TEST FIXTURE ============== */

class IFDSUninitializedVariablesTest : public ::testing::Test {
protected:

  using UninitCompactResult_t =
      std::tuple<std::string, std::size_t, std::string, int64_t>;


  const std::string PathToLlFiles =
      unittest::PathToLLTestFiles + "uninitialized_variables/";
  const std::set<std::string> EntryPoints = {"main"};

  unique_ptr<ProjectIRDB> IRDB;
  unique_ptr<LLVMTypeHierarchy> TH;
  unique_ptr<LLVMBasedICFG> ICFG;
  unique_ptr<LLVMPointsToInfo> PT;
  unique_ptr<IFDSUninitializedVariables> UninitProblem;

  IFDSUninitializedVariablesTest() = default;
  ~IFDSUninitializedVariablesTest() override = default;


  void SetUp() override {
    boost::log::core::get()->set_logging_enabled(false);
    ValueAnnotationPass::resetValueID();
  }

  void TearDown() override {}

  IFDSUninitializedVariables::uninit_results_t
  doAnalysis(const std::string &LlvmFilePath, bool PrintDump = false) {

    auto IR_Files = {PathToLlFiles + LlvmFilePath};
    IRDB = std::make_unique<ProjectIRDB>(IR_Files, IRDBOptions::WPA);
    LLVMTypeHierarchy TH(*IRDB);
    LLVMPointsToSet PT(*IRDB);
    LLVMBasedICFG ICFG(*IRDB, CallGraphAnalysisType::OTF, {"main"}, &TH, &PT,
                       Soundness::Soundy, /*IncludeGlobals*/ true);

    IFDSUninitializedVariables UninitProblem(IRDB.get(), &TH, &ICFG, &PT, {"main"});

    IFDSSolver Solver(UninitProblem);
    Solver.solve();
    if (PrintDump) {
      IRDB->print();
      ICFG.print();
      Solver.dumpResults();
    }

    return UninitProblem.getUnitializedResults(Solver.getSolverResults());
  }




  static void compareResults(IFDSUninitializedVariables::uninit_results_t &Results,
                             std::set<UninitCompactResult_t> &GroundTruth) {
    std::set<UninitCompactResult_t> RelevantResults;
    for (auto G : GroundTruth) {
      std::string FName = std::get<0>(G);
      unsigned Line = std::get<1>(G);
      if (Results.find(FName) != Results.end()) {
        if (auto It = Results[FName].find(Line); It != Results[FName].end()) {
          for (const auto &VarToVal : It->second.variableToValue) {
            RelevantResults.emplace(FName, Line, VarToVal.first,
                                    VarToVal.second);
          }
        }
      }
    }
    EXPECT_EQ(RelevantResults, GroundTruth);
  }
}; // Test Fixture

TEST_F(IFDSUninitializedVariablesTest, UninitTest_01_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("all_uninit_cpp_dbg.ll");
  // all_uninit.cpp does not contain undef-uses
  std::set<UninitCompactResult_t> GroundTruth;
  compareResults(Results, GroundTruth);
}

TEST_F(IFDSUninitializedVariablesTest, UninitTest_02_SHOULD_LEAK) {
  auto Results = doAnalysis("binop_uninit_cpp_dbg.ll");

  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "i", "?");
  // binop_uninit uses uninitialized variable i in 'int j = i + 10;'
  GroundTruth.emplace("main", 3, "j", "i + 10");
  compareResults(Results, GroundTruth);
}
TEST_F(IFDSUninitializedVariablesTest, UninitTest_03_SHOULD_LEAK) {
  auto Results = doAnalysis("callnoret_c_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "a", "?");
  // callnoret uses uninitialized variable a in 'return a + 10;' of addTen(int)
  GroundTruth.emplace("main", 4, "a", "a + 10");
  GroundTruth.emplace("main", 5, "d", "?");

  compareResults(Results,GroundTruth);
}

TEST_F(IFDSUninitializedVariablesTest, UninitTest_04_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("ctor_default_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // ctor.cpp does not contain undef-uses
  compareResults(Results, GroundTruth);
}

TEST_F(IFDSUninitializedVariablesTest, UninitTest_05_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("struct_member_init_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // struct_member_init.cpp does not contain undef-uses
  compareResults(Results, GroundTruth);
}
TEST_F(IFDSUninitializedVariablesTest, UninitTest_06_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("ctor_default_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // struct_member_uninit.cpp does not contain undef-uses
  compareResults(Results, GroundTruth);
}
/****************************************************************************************
 * fails due to field-insensitivity + struct ignorance + clang compiler hacks
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_07_SHOULD_LEAK) {
  Initialize({pathToLLFiles + "struct_member_uninit2_cpp_dbg.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve();
  // struct_member_uninit2.cpp contains a use of the uninitialized field _x.b
  map<int, set<string>> GroundTruth;
  // %5 = load i16, i16* %4; %4 is the uninitialized struct-member _x.b
  GroundTruth[4] = {"3"};



  compareResults(GroundTruth);
}
*****************************************************************************************/
TEST_F(IFDSUninitializedVariablesTest, UninitTest_08_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("global_variable_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // global_variable.cpp does not contain undef-uses
  compareResults(Results, GroundTruth);
}
/****************************************************************************************
 * failssince @i is uninitialized in the c++ code, but initialized in the
 * LLVM-IR
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_09_SHOULD_LEAK) {
  Initialize({pathToLLFiles + "global_variable_cpp_dbg.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve();
  // global_variable.cpp does not contain undef-uses
  map<int, set<string>> GroundTruth;
  // load i32, i32* @i
  GroundTruth[5] = {"0"};
  compareResults(GroundTruth);
}
*****************************************************************************************/
TEST_F(IFDSUninitializedVariablesTest, UninitTest_10_SHOULD_LEAK) {

  auto Results = doAnalysis("return_uninit_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "a", "?");
  compareResults(Results, GroundTruth);
}

TEST_F(IFDSUninitializedVariablesTest, UninitTest_11_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("sanitizer_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // all undef-uses are sanitized;
  // However, the uninitialized variable j is read, which causes the analysis to
  // report an undef-use
  GroundTruth.emplace("main", 4, "i", "?");
  GroundTruth.emplace("main", 4, "j", "?");
  compareResults(Results, GroundTruth);
}
//---------------------------------------------------------------------
// Not relevant any more; Test case covered by UninitTest_11
//---------------------------------------------------------------------
/* TEST_F(IFDSUninitializedVariablesTest, UninitTest_12_SHOULD_LEAK) {

  Initialize({pathToLLFiles + "sanitizer_uninit_cpp_dbg.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, true); Solver.solve();
  // The sanitized value is not used always => the phi-node is "tainted"
  map<int, set<string>> GroundTruth;
  GroundTruth[6] = {"2"};
  GroundTruth[13] = {"2"};
  compareResults(GroundTruth);
}
*/
TEST_F(IFDSUninitializedVariablesTest, UninitTest_13_SHOULD_NOT_LEAK) {
  auto Results = doAnalysis("ctor_default_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  // The undef-uses do not affect the program behaviour, but are of course still
  // found and reported
  GroundTruth.emplace("main", 3, "t", "?");
  compareResults(Results, GroundTruth);
}
TEST_F(IFDSUninitializedVariablesTest, UninitTest_14_SHOULD_LEAK) {
  auto Results = doAnalysis("uninit_c_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 3, "a", "?");
  GroundTruth.emplace("main", 4, "b", "?");
  GroundTruth.emplace("main", 7, "e", "a*b");
  compareResults(Results, GroundTruth);
}
/****************************************************************************************
 * Fails probably due to field-insensitivity
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_15_SHOULD_LEAK) {
  Initialize({pathToLLFiles + "dyn_mem_cpp_dbg.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve(); map<int, set<string>>
GroundTruth;
  // TODO remove GT[14] and GT[13]
  GroundTruth[14] = {"3"};
  GroundTruth[13] = {"2"};
  GroundTruth[15] = {"13", "14"};

  GroundTruth[35] = {"4"};
  GroundTruth[38] = {"35"};

  GroundTruth[28] = {"2"};
  GroundTruth[29] = {"3"};
  GroundTruth[30] = {"28", "29"};

  GroundTruth[33] = {"30"};

  // Analysis detects false positive at %12:

  // store i32* %3, i32** %6, align 8, !dbg !28
  // %12 = load i32*, i32** %6, align 8, !dbg !29


  compareResults(GroundTruth);
}
*****************************************************************************************/
TEST_F(IFDSUninitializedVariablesTest, UninitTest_16_SHOULD_LEAK) {
  auto Results = doAnalysis("growing_example_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "i", "?");
  GroundTruth.emplace("main", 3, "j", "?");
  GroundTruth.emplace("main", 4, "k", "?");
  GroundTruth.emplace("main", 3, "k", "?+12");
  compareResults(Results, GroundTruth);
}

/****************************************************************************************
 * Fails due to struct ignorance; general problem with field sensitivity: when
 * all structs would be treated as uninitialized per default, the analysis would
 * not be able to detect correct constructor calls
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_17_SHOULD_LEAK) {

  Initialize({pathToLLFiles + "struct_test_cpp.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve();

  map<int, set<string>> GroundTruth;
  // printf should leak both parameters => fails

  GroundTruth[8] = {"5", "7"};
  compareResults(GroundTruth);
}
*****************************************************************************************/
/****************************************************************************************
 * Fails, since the analysis is not able to detect memcpy calls
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_18_SHOULD_NOT_LEAK) {

  Initialize({pathToLLFiles + "array_init_cpp.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve();

  map<int, set<string>> GroundTruth;
  //

  compareResults(GroundTruth);
}
*****************************************************************************************/
/****************************************************************************************
 * fails due to missing alias information (and missing field/array element
 *information)
 *
*****************************************************************************************
TEST_F(IFDSUninitializedVariablesTest, UninitTest_19_SHOULD_NOT_LEAK) {

  Initialize({pathToLLFiles + "array_init_simple_cpp.ll"});
  IFDSSolver<IFDSUninitializedVariables::n_t,
IFDSUninitializedVariables::d_t,IFDSUninitializedVariables::f_t,IFDSUninitializedVariables::t_t,IFDSUninitializedVariables::v_t,IFDSUninitializedVariables::i_t>
Solver(*UninitProblem, false); Solver.solve();

  map<int, set<string>> GroundTruth;
  compareResults(GroundTruth);
}
*****************************************************************************************/
TEST_F(IFDSUninitializedVariablesTest, UninitTest_20_SHOULD_LEAK) {
  auto Results = doAnalysis("recursion_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "i", "?");
  GroundTruth.emplace("main", 3, "j", "i");
  compareResults(Results, GroundTruth);
}
TEST_F(IFDSUninitializedVariablesTest, UninitTest_21_SHOULD_LEAK) {
  auto Results = doAnalysis("virtual_call_cpp_dbg.ll");
  std::set<UninitCompactResult_t> GroundTruth;
  GroundTruth.emplace("main", 2, "i", "?");
  GroundTruth.emplace("main", 3, "baz", "?");
  compareResults(Results, GroundTruth);
}
int main(int Argc, char **Argv) {
  ::testing::InitGoogleTest(&Argc, Argv);
  return RUN_ALL_TESTS();
}
