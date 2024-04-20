
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

void RunTaskLPi(bool isMC, std::string inFileName, std::string suffix = "0", int nFiles = 1) {
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

    // PID mathis task
    AliAnalysisTaskPIDResponse *PIDTask =
        reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(
            Form(".x %s(%d, kTRUE, kTRUE, \"1\")",
                 gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isMC)));

    // femto task
    AliAnalysisTaskLambdaPion *task =
        reinterpret_cast<AliAnalysisTaskLambdaPion *>(gInterpreter->ProcessLine(
            Form(".x %s(%d,\"kHM\",128,0,\"0\",false,false,AliAnalysisTaskLambdaPion::PCSettings::NoPC,false,2,\"%s\")", // old
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
    chain->Add(inFileName.data());

    mgr->StartAnalysis("local", chain);
}
