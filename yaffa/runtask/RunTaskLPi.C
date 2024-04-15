
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskLambdaPion.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMultSelectionTask.h"
#include "AliPhysicsSelectionTask.h"
#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

void RunLocalAnalysis(bool isMC = false, std::string suffix = "0", int nFiles = 1) {
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");

    // create an input handler
    AliAODInputHandler *inputH = new AliAODInputHandler();
    mgr->SetInputEventHandler(inputH);

    // Physics Selection task
    AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(
        gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d,kTRUE)", isMC)));

    // mult task
    AliMultSelectionTask *taskMult = reinterpret_cast<AliMultSelectionTask *>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));

    // // pid response task
    // AliAnalysisTaskPIDResponse *PIDTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
    //     Form(".x %s(%d, kTRUE, kTRUE, \"1\")",
    //              gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isMC));

        // PID mathis task
    AliAnalysisTaskPIDResponse *PIDTask =
        reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(
            Form(".x %s(%d, kTRUE, kTRUE, \"1\")",
                 gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isMC)));

    // femto task
    AliAnalysisTaskLambdaPion *task =
        reinterpret_cast<AliAnalysisTaskLambdaPion *>(gInterpreter->ProcessLine(
            Form(".x %s(%d,\"kHM\",128,0,\"0\",false,false,AliAnalysisTaskLambdaPion::PCSettings::NoPC,false,\"%s\")", // old
            // Form(".x %s(%d,\"kHM\",128,false,AliAnalysisTaskLambdaPion::PCSettings::NoPC,false,\"%s\")", // new
                gSystem->ExpandPathName("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/AddTaskFemtoLambdaPion.C"), isMC, suffix.data())));

    task->SetPairCleaner(AliAnalysisTaskLambdaPion::PCSettings::NewPC);
    task->SetExcludeDausOf({3224, 3114});

    if (!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(1);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 250);

    // if you want to run locally, we need to define some input
    TChain *chain = new TChain("aodTree");
    if (isMC) {
        // chain->Add("~/an/LPi/run_local/sim/2021/LHC21g2a/274442/AOD/001/AliAOD.root");

        std::vector<const char*> files = {
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/043/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/007/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/017/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/004/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/002/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/041/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/045/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/047/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/022/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/050/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/016/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/026/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/006/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/010/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/030/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/028/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/001/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/027/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/032/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/048/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/019/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/034/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/003/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/049/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/036/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/031/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/040/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/023/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/008/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/005/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/018/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/044/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/033/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/011/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/021/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/046/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/042/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/015/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/012/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/024/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/029/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/009/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/037/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/035/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/039/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/038/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/020/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/013/AliAOD.root",
            "/data/grid/sim/2021/LHC21b3a/294925/AOD/025/AliAOD.root",
        };

        for (size_t iName = 0; iName < nFiles; iName++) chain->Add(files[iName]);

    } else {
        chain->Add("/home/daniel/an/LPi/run_local/data/LHC18b/000285064/pass2/AOD/001/AliAOD.root");
    }

    mgr->StartAnalysis("local", chain);
}
