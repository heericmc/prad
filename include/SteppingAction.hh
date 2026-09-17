#pragma once
#include "G4UserSteppingAction.hh"

// SteppingAction is retained for compilation but is no longer registered.
// Hit collection is handled by TrackerSD / EventAction.
class SteppingAction : public G4UserSteppingAction
{
public:
    void UserSteppingAction(const G4Step*) override {}
};
