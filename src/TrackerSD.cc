#include "TrackerSD.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

TrackerSD::TrackerSD(G4String name, G4String hcName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(hcName);
}

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

    hce->AddHitsCollection(fHCID, fHitsCollection);
    fTrackHitMap.clear();
}

G4bool TrackerSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // Only the primary proton's own trajectory counts toward the coincidence trigger --
    // we're measuring whether THIS proton (direct, or scattered off the uranium into the
    // FOV) continues on to the tracker, not whether some spallation product it generated
    // along the way happens to land on silicon. Hadronic secondaries (kept alive by
    // StackingAction) are transported for correct energy bookkeeping but never hit here.
    const G4Track* track = step->GetTrack();
    if (track->GetTrackID() != 1) return false;

    // Determine which plane from the physical volume name
    const G4String pvName =
        step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
    G4int planeID = -1;
    if      (pvName == "SSD0PV")  planeID = 0;   // tracking layer 0
    else if (pvName == "SSD1PV")  planeID = 1;   // tracking layer 1
    else if (pvName == "SSD2PV")  planeID = 2;   // tracking layer 2
    else if (pvName == "ScintPV") planeID = 3;
    else if (pvName == "PolyPV")  planeID = 4;
    else if (pvName == "StripXPV") planeID = 5;   // strip hodoscope, x-measuring layer
    else if (pvName == "StripYPV") planeID = 6;   // strip hodoscope, y-measuring layer
    else if (pvName == "V3CollimatorTubePV") planeID = 7;   // wall-thickness investigation, 2026-09-15
    else return false;

    const G4int    trackID = track->GetTrackID();
    const G4double edep    = step->GetTotalEnergyDeposit();

    auto key = std::make_pair(trackID, planeID);
    auto it  = fTrackHitMap.find(key);

    if (it == fTrackHitMap.end()) {
        // First step of this track in this plane: record entry kinematics
        auto* hit = new TrackerHit();
        hit->SetPlaneID(planeID);
        hit->SetTrackID(trackID);
        hit->SetParticleName(track->GetParticleDefinition()->GetParticleName());

        const G4StepPoint* pre = step->GetPreStepPoint();
        hit->SetPosition(pre->GetPosition());
        hit->SetMomDir(pre->GetMomentumDirection());
        hit->SetKineticEnergy(pre->GetKineticEnergy());
        hit->SetTime(pre->GetGlobalTime());
        hit->AddEdep(edep);
        if (planeID == 5 || planeID == 6)   // strip hodoscope: which strip fired
            hit->SetStripID(pre->GetTouchable()->GetCopyNumber());

        const G4int idx = fHitsCollection->entries();
        fHitsCollection->insert(hit);
        fTrackHitMap[key] = idx;
    } else {
        // Accumulate energy deposit from subsequent steps inside this plane
        (*fHitsCollection)[it->second]->AddEdep(edep);
    }

    return true;
}
