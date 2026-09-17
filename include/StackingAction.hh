#pragma once
#include "G4UserStackingAction.hh"
#include "globals.hh"

class G4Track;

// Kills secondary electrons, positrons, and photons at creation to save CPU. None of
// them can penetrate the 25 cm polyethylene filter to reach the scintillator, so they
// can't create a coincidence trigger, and killing them doesn't touch the proton's
// continuous dE/dx (which stops protons in the uranium and sets the ~200 MeV poly
// cutoff). Protons, neutrons, and other hadrons are kept alive for correct energy
// bookkeeping, though TrackerSD::ProcessHits only ever records the primary track
// (ID 1), so a secondary hadron reaching the tracker can't register a hit anyway.
class StackingAction : public G4UserStackingAction
{
public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
};
