#include "G4RunManagerFactory.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <string>
#include <cstdlib>
#include <chrono>

int main(int argc, char** argv)
{
    // Mode is picked from the macro filename:
    //   imaging_all.mac  -> shadow run then open-field run, one execution
    //   imaging_open.mac -> open-field reference (tracker only, B field on)
    //   imaging.mac      -> shadow run (uranium + tracker, B field on)
    //   vis_imaging.mac  -> Qt GUI
    //   collimator_export.mac    -> batch-only VRML export of the collimator concept geometry
    //   (no macro)       -> Qt GUI, double-click default
    // All modes but the collimator concept use the AP9 spectrum beam.
    const std::string macroFile = (argc > 1) ? std::string(argv[1]) : std::string();
    const bool kCombined    = macroFile.find("imaging_all")  != std::string::npos;
    const bool kImagingOpen = macroFile.find("imaging_open") != std::string::npos;
    const bool kCollimator  = macroFile.find("collimator") != std::string::npos;
    const bool kVisMode     = !kCombined && !kImagingOpen && !kCollimator &&
                              (macroFile.empty() || macroFile.find("vis") != std::string::npos);

    // Double-click (empty macro) opens the full standoff world with no uranium so fired
    // protons are visible crossing the open standoff. Combined mode starts in shadow mode
    // and switches to open after the first BeamOn. Collimator macros (no "vis" substring) fall
    // through to the plain batch-execute branch below, so no interactive Qt session
    // blocks the export.
    const bool kVisOpen = kVisMode && macroFile.empty();
    DetectorConstruction::SetMode(kCollimator ? DetectorConstruction::Mode::kCollimatorConcept
        : (kImagingOpen || kVisOpen)
            ? DetectorConstruction::Mode::kImagingOpen
            : DetectorConstruction::Mode::kImaging);

    // sVisMode fires the same production beam but tells RunAction to skip CSV output,
    // so vis runs never clobber production data.
    PrimaryGeneratorAction::sVisMode = kVisMode;

    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);

    auto* mtRunManager = dynamic_cast<G4MTRunManager*>(runManager);
    if (mtRunManager)
        mtRunManager->SetNumberOfThreads(18);

    runManager->SetUserInitialization(new DetectorConstruction());

    G4PhysListFactory physFactory;
    G4VModularPhysicsList* phys = physFactory.GetReferencePhysList("FTFP_BERT");
    runManager->SetUserInitialization(phys);

    runManager->SetUserInitialization(new ActionInitialization());

    // Geant4's default seed is fixed, so without this every invocation replays the same
    // event sequence and separate chunked runs would just be duplicates. SEED env var
    // overrides for reproducible per-chunk seeds; otherwise seed off the clock.
    {
        const char* s = std::getenv("SEED");
        const long seed = s ? std::atol(s)
            : static_cast<long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        G4Random::setTheSeed(seed);
        G4cout << "[MAIN] RNG seed = " << seed << G4endl;
    }

    runManager->Initialize();

    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    auto* uiManager = G4UImanager::GetUIpointer();

    if (kCombined) {
        // Shadow run (uranium present).
        uiManager->ApplyCommand("/control/execute " + G4String(argv[1]));
        // Open-field run (no uranium): rebuild geometry then re-execute same mac.
        DetectorConstruction::SetMode(DetectorConstruction::Mode::kImagingOpen);
        runManager->ReinitializeGeometry();
        uiManager->ApplyCommand("/control/execute " + G4String(argv[1]));
    } else if (kVisMode) {
        auto* ui = new G4UIExecutive(argc, argv);
        if (!macroFile.empty())
            uiManager->ApplyCommand("/control/execute " + G4String(argv[1]));
        else
            uiManager->ApplyCommand("/control/execute vis_world.mac");   // double-click default
        ui->SessionStart();
        delete ui;
    } else {
        uiManager->ApplyCommand("/control/execute " + G4String(argv[1]));
    }

    delete visManager;
    delete runManager;
    return 0;
}
