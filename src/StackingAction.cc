#include "StackingAction.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include <cstdlib>
#include <string>

namespace {
    // Enabled by default; set env KILL_EM=0 to disable (for A/B validation).
    const bool kKillEM = [] {
        const char* e = std::getenv("KILL_EM");
        return (e == nullptr) || (std::string(e) != "0");
    }();
}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
    // Only ever affects secondaries (primaries are protons); kill e-, e+, gamma.
    if (kKillEM && track->GetParentID() > 0) {
        const G4ParticleDefinition* pd = track->GetParticleDefinition();
        if (pd == G4Gamma::Definition()    ||
            pd == G4Electron::Definition() ||
            pd == G4Positron::Definition())
            return fKill;
    }
    return fUrgent;
}
