#pragma once
#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"
#include <map>
#include <utility>

class G4Step;
class G4HCofThisEvent;

class TrackerSD : public G4VSensitiveDetector
{
public:
    TrackerSD(G4String name, G4String hcName);
    ~TrackerSD() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

private:
    TrackerHitsCollection* fHitsCollection = nullptr;
    G4int fHCID = -1;

    // (trackID, planeID) -> index in hit collection
    std::map<std::pair<G4int,G4int>, G4int> fTrackHitMap;
};
