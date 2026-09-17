#pragma once
#include "G4UserRunAction.hh"
#include "G4Threading.hh"
#include <fstream>

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    RunAction()  = default;
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*)   override;

    static std::ofstream* GetFile()  { return sFile; }
    static G4Mutex&       GetMutex() { return sMutex; }

    // Prescaled dump of every fired primary (not just the triggered ones), so the
    // sampled beam (pitch-angle distribution, gyrophase, energy, source position)
    // can be plotted and the PAD sampling checked. Separate file and mutex so it
    // never contends with the trigger CSV.
    static std::ofstream* GetBeamFile()  { return sBeamFile; }
    static G4Mutex&       GetBeamMutex() { return sBeamMutex; }

    // 2026-09-15: wall-thickness investigation. Dumps every primary that deposits energy
    // in the tube collimator's wall (V3CollimatorTubePV, TrackerHit planeID 7), regardless
    // of whether the full coincidence trigger fires -- independent file/mutex so it never
    // contends with the trigger CSV. Only meaningful when Det::kCollStyle == kTube.
    static std::ofstream* GetWallFile()  { return sWallFile; }
    static G4Mutex&       GetWallMutex() { return sWallMutex; }

private:
    static std::ofstream* sFile;
    static G4Mutex        sMutex;
    static std::ofstream* sBeamFile;
    static G4Mutex        sBeamMutex;
    static std::ofstream* sWallFile;
    static G4Mutex        sWallMutex;
};
