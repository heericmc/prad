#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Threading.hh"
#include "G4Run.hh"
#include <fstream>

std::ofstream* RunAction::sFile      = nullptr;
G4Mutex        RunAction::sMutex     = G4MUTEX_INITIALIZER;
std::ofstream* RunAction::sBeamFile  = nullptr;
G4Mutex        RunAction::sBeamMutex = G4MUTEX_INITIALIZER;
std::ofstream* RunAction::sWallFile  = nullptr;
G4Mutex        RunAction::sWallMutex = G4MUTEX_INITIALIZER;

static const char* kWallHeader =
    "event_id,entry_x_mm,entry_y_mm,entry_z_mm,entry_ke_MeV,angle_from_z_deg,wall_edep_MeV,punched_through\n";

// Prescaled sample of every fired primary (post cone + cos-incidence acceptance), i.e.
// the beam as launched, independent of whether it reached the detector. This is the
// distribution to plot for PAD / gyrophase / source-position checks; the trigger CSV's
// birth columns are the same quantities but for the biased triggered subset.
static const char* kBeamHeader =
    "src_x_mm,src_y_mm,src_z_mm,ke_MeV,ux,uy,uz,pitch_deg,gyro_deg\n";

static const char* kHeader =
    "event_id,"
    "ssd0_x_mm,ssd0_y_mm,ssd0_edep_keV,ssd0_time_ns,"
    "ssd1_x_mm,ssd1_y_mm,ssd1_edep_keV,ssd1_ke_MeV,ssd1_time_ns,"
    "ssd2_x_mm,ssd2_y_mm,ssd2_edep_keV,ssd2_ke_MeV,ssd2_time_ns,"
    "poly_edep_MeV,scint_hit,"
    "geo_theta_x_deg,geo_theta_y_deg,geo_theta_deg,geo_phi_deg,"
    "mom_theta_deg,mom_phi_deg,"
    "mom0_theta_deg,mom0_phi_deg,"
    // Birth (primary-vertex) state, appended at the end so column-name lookups elsewhere
    // keep working unchanged.
    "src_x_mm,src_y_mm,src_z_mm,birth_ke_MeV,birth_pitch_deg,birth_gyro_deg\n";

// X-Y strip hodoscope variant (Det::kUseStripHodoscope): poly_edep_MeV/scint_hit are
// replaced by the fired strip index + deposited energy in each crossed layer.
static const char* kHeaderStripHodo =
    "event_id,"
    "ssd0_x_mm,ssd0_y_mm,ssd0_edep_keV,ssd0_time_ns,"
    "ssd1_x_mm,ssd1_y_mm,ssd1_edep_keV,ssd1_ke_MeV,ssd1_time_ns,"
    "ssd2_x_mm,ssd2_y_mm,ssd2_edep_keV,ssd2_ke_MeV,ssd2_time_ns,"
    "stripx_id,stripx_edep_MeV,stripy_id,stripy_edep_MeV,"
    "geo_theta_x_deg,geo_theta_y_deg,geo_theta_deg,geo_phi_deg,"
    "mom_theta_deg,mom_phi_deg,"
    "mom0_theta_deg,mom0_phi_deg,"
    "src_x_mm,src_y_mm,src_z_mm,birth_ke_MeV,birth_pitch_deg,birth_gyro_deg\n";

void RunAction::BeginOfRunAction(const G4Run*)
{
    if (!G4Threading::IsMasterThread()) return;
    if (PrimaryGeneratorAction::sVisMode) return;   // vis runs never touch production CSVs

    using Mode = DetectorConstruction::Mode;
    const bool shadow = (DetectorConstruction::GetMode() == Mode::kImaging);
    sFile = new std::ofstream(shadow ? "imaging_shadow.csv" : "imaging_open.csv");
    if (sFile) *sFile << (DetectorConstruction::kUseStripHodoscope ? kHeaderStripHodo : kHeader);

    // BEAM_DUMP_PRESCALE = 0 disables the beam dump; otherwise 1 primary in N is written.
    if (PrimaryGeneratorAction::BeamDumpPrescale() > 0) {
        sBeamFile = new std::ofstream(shadow ? "beam_sample_shadow.csv"
                                             : "beam_sample_open.csv");
        if (sBeamFile) *sBeamFile << kBeamHeader;
    }

    // Wall-thickness investigation (2026-09-15): only meaningful for the tube collimator,
    // but harmless (empty, header-only) otherwise -- no volume is ever named
    // V3CollimatorTubePV in the lattice build, so TrackerSD's planeID 7 branch never fires.
    sWallFile = new std::ofstream(shadow ? "wall_hits_shadow.csv" : "wall_hits_open.csv");
    if (sWallFile) *sWallFile << kWallHeader;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    if (!G4Threading::IsMasterThread()) return;
    if (PrimaryGeneratorAction::sVisMode) return;   // vis runs write no CSV/stats

    if (sFile) {
        sFile->close();
        delete sFile;
        sFile = nullptr;
    }
    if (sBeamFile) {
        sBeamFile->close();
        delete sBeamFile;
        sBeamFile = nullptr;
    }
    if (sWallFile) {
        sWallFile->close();
        delete sWallFile;
        sWallFile = nullptr;
    }

    using Mode = DetectorConstruction::Mode;
    switch (DetectorConstruction::GetMode()) {
        case Mode::kImaging: {
            std::ofstream stats("imaging_shadow_stats.csv");
            stats << "n_fired\n" << run->GetNumberOfEvent() << '\n';
            break;
        }
        case Mode::kImagingOpen: {
            std::ofstream stats("imaging_open_stats.csv");
            stats << "n_fired\n" << run->GetNumberOfEvent() << '\n';
            break;
        }
    }
}
