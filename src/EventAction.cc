#include "EventAction.hh"
#include "RunAction.hh"
#include "TrackerHit.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

#include <cmath>

void EventAction::EndOfEventAction(const G4Event* event)
{
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID("TrackerSD/TrackerHC");

    auto* hce = event->GetHCofThisEvent();
    if (!hce) return;

    auto* hc = dynamic_cast<TrackerHitsCollection*>(hce->GetHC(fHCID));
    if (!hc) return;

    TrackerHit* h0 = nullptr;   // tracking layer 0  planeID 0
    TrackerHit* h1 = nullptr;   // tracking layer 1  planeID 1
    TrackerHit* h2 = nullptr;   // tracking layer 2  planeID 2
    TrackerHit* hs = nullptr;   // scintillator      planeID 3  (unused if kUseStripHodoscope)
    TrackerHit* hp = nullptr;   // poly filter        planeID 4  (unused if kUseStripHodoscope)
    TrackerHit* hStripX = nullptr;   // strip hodoscope, x-layer  planeID 5
    TrackerHit* hStripY = nullptr;   // strip hodoscope, y-layer  planeID 6
    TrackerHit* hWall   = nullptr;   // v3 tube collimator wall   planeID 7 (investigation only)

    for (G4int i = 0; i < hc->entries(); ++i) {
        TrackerHit* h = (*hc)[i];
        switch (h->GetPlaneID()) {
            case 0: if (!h0) h0 = h; break;
            case 1: if (!h1) h1 = h; break;
            case 2: if (!h2) h2 = h; break;
            case 3: if (!hs) hs = h; break;
            case 4: if (!hp) hp = h; break;
            case 5: if (!hStripX) hStripX = h; break;
            case 6: if (!hStripY) hStripY = h; break;
            case 7: if (!hWall) hWall = h; break;
        }
    }

    // Wall-thickness investigation (2026-09-15): dump every primary that touches the tube
    // collimator's wall, REGARDLESS of whether it goes on to trigger -- independent of the
    // coincidence gate below, so this also captures protons the wall successfully stops.
    // "punched_through" = also reached all 3 Si layers (the tracking coincidence core,
    // independent of the trigger-stage requirement) -- a simple proxy for "the wall didn't
    // stop it enough to matter for this measurement."
    if (hWall) {
        G4AutoLock wlock(&RunAction::GetWallMutex());
        auto* wout = RunAction::GetWallFile();
        if (wout && wout->is_open()) {
            const G4ThreeVector& wp = hWall->GetPosition();
            const G4ThreeVector& wd = hWall->GetMomDir();
            const G4double angleFromZ = std::acos(std::min(1.0, std::max(-1.0, wd.z()))) / deg;
            const G4int punched = DetectorConstruction::kTwoLayerViz
                ? (h0 && h1) ? 1 : 0
                : (h0 && h1 && h2) ? 1 : 0;
            *wout << event->GetEventID()          << ','
                  << wp.x()/mm                     << ','
                  << wp.y()/mm                     << ','
                  << wp.z()/mm                     << ','
                  << hWall->GetKineticEnergy()/MeV << ','
                  << angleFromZ                    << ','
                  << hWall->GetEdep()/MeV          << ','
                  << punched
                  << '\n';
        }
    }

    // Quadruple coincidence: all 3 Si tracking layers plus the trigger stage, whichever
    // it is -- the scintillator (which requires the proton to clear the 25 cm poly
    // filter, E >~ 200 MeV) or, for the strip hodoscope variant, both crossed strip
    // layers (no hard energy threshold -- see StripHodo namespace comment).
    if (DetectorConstruction::kTwoLayerViz) {
        if (!h0 || !h1) return;   // no layer 2 exists in this geometry -- see the header
    } else {
        if (!h0 || !h1 || !h2) return;
    }
    if (DetectorConstruction::kUseStripHodoscope) {
        if (!hStripX || !hStripY) return;
    } else {
        if (!hs) return;
    }

    // geo_theta comes from layers 0 -> 1 -- this is what defines the back-projection to
    // the uranium plane. Layer 2 and the scintillator are coincidence-only.
    const G4ThreeVector d = h1->GetPosition() - h0->GetPosition();
    const G4double geo_tx  = std::atan2(d.x(), d.z()) / deg;
    const G4double geo_ty  = std::atan2(d.y(), d.z()) / deg;
    const G4double geo_th  = std::atan2(std::sqrt(d.x()*d.x() + d.y()*d.y()), d.z()) / deg;
    G4double geo_phi = std::atan2(d.y(), d.x()) / deg;
    if (geo_phi < 0.) geo_phi += 360.;

    // mom_theta is the momentum entering layer-1 silicon, i.e. downstream of most of the
    // material the 0->1 chord is sensitive to, so a scatter kick in that span tilts both
    // the chord and this direction and partly cancels out. The residual (geo - mom) reads
    // low and is not a clean sigma_MCS -- kept only as a cross-check.
    const G4ThreeVector& md = h1->GetMomDir();
    const G4double mom_th  = std::acos(std::min(1.0, md.z())) / deg;
    G4double mom_phi = std::atan2(md.y(), md.x()) / deg;
    if (mom_phi < 0.) mom_phi += 360.;

    // mom0_theta is the primary's launch direction straight off the event vertex, before
    // any tracker material, so (geo_theta_0->1 - mom0) is the true back-projection error
    // with no shared scattering to cancel out.
    G4ThreeVector md0(0., 0., 1.);
    // Birth state of the primary: launch point, launch energy, and launch direction in
    // the sampler's (pitch, gyrophase) coordinates. src_* is the true source position, so
    // (back-projected x_src - src_x_mm) is the back-projection error with no model in the way.
    G4double src_x = 0., src_y = 0., src_z = 0., birth_ke = 0.;
    if (const auto* pv = event->GetPrimaryVertex(0)) {
        src_x = pv->GetPosition().x() / mm;
        src_y = pv->GetPosition().y() / mm;
        src_z = pv->GetPosition().z() / mm;
        if (const auto* pp = pv->GetPrimary(0)) {
            md0      = pp->GetMomentum().unit();
            birth_ke = pp->GetKineticEnergy() / MeV;
        }
    }
    const G4double mom0_th  = std::acos(std::min(1.0, md0.z())) / deg;
    G4double mom0_phi = std::atan2(md0.y(), md0.x()) / deg;
    if (mom0_phi < 0.) mom0_phi += 360.;

    // B is along +y: pitch = angle from +y, gyrophase measured in the x-z plane. Same
    // convention as PrimaryGeneratorAction::DumpBeamSample, so the triggered subset can
    // be histogrammed directly against the full fired beam.
    const G4double birth_pitch = std::acos(std::max(-1.0, std::min(1.0, md0.y()))) / deg;
    const G4double birth_gyro  = std::atan2(md0.x(), md0.z()) / deg;

    G4AutoLock lock(&RunAction::GetMutex());
    auto* out = RunAction::GetFile();
    if (!out || !out->is_open()) return;

    *out << event->GetEventID()        << ','
         << h0->GetPosition().x()/mm   << ','
         << h0->GetPosition().y()/mm   << ','
         << h0->GetEdep()/keV          << ','
         << h0->GetTime()/ns           << ','
         << h1->GetPosition().x()/mm   << ','
         << h1->GetPosition().y()/mm   << ','
         << h1->GetEdep()/keV          << ','
         << h1->GetKineticEnergy()/MeV << ','
         << h1->GetTime()/ns           << ',';

    if (DetectorConstruction::kTwoLayerViz) {
        // No layer 2 in this geometry -- five empty fields keep the column count (and
        // every downstream column's position) identical to the 3-layer CSV format, so
        // imaging_analysis.m's fixed-column reads don't need a special case.
        *out << ",,,,,";
    } else {
        *out << h2->GetPosition().x()/mm   << ','
             << h2->GetPosition().y()/mm   << ','
             << h2->GetEdep()/keV          << ','
             << h2->GetKineticEnergy()/MeV << ','
             << h2->GetTime()/ns           << ',';
    }

    if (DetectorConstruction::kUseStripHodoscope) {
        // Strip hodoscope trigger stage: which strip fired + how much energy it deposited,
        // in each crossed layer. imaging_analysis.m turns (stripx_id, stripy_id) into the
        // 1 cm^2 hit cell and (stripx_edep_MeV + stripy_edep_MeV) into a proton energy
        // estimate via dE/dx.
        *out << hStripX->GetStripID()   << ','
             << hStripX->GetEdep()/MeV  << ','
             << hStripY->GetStripID()   << ','
             << hStripY->GetEdep()/MeV  << ',';
    } else {
        const G4double poly_edep_MeV = hp ? hp->GetEdep() / MeV : 0.;
        const G4int    scint_hit     = hs ? 1 : 0;
        *out << poly_edep_MeV << ',' << scint_hit << ',';
    }

    *out << geo_tx  << ',' << geo_ty   << ','
         << geo_th  << ',' << geo_phi  << ','
         << mom_th  << ',' << mom_phi  << ','
         << mom0_th << ',' << mom0_phi << ','
         << src_x   << ',' << src_y    << ',' << src_z << ','
         << birth_ke << ',' << birth_pitch << ',' << birth_gyro
         << '\n';
}
