#pragma once
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"

class TrackerHit : public G4VHit
{
public:
    TrackerHit()  = default;
    ~TrackerHit() override = default;

    void* operator new(size_t);
    void  operator delete(void*);

    void SetPlaneID       (G4int id)          { fPlaneID = id; }
    void SetTrackID       (G4int id)          { fTrackID = id; }
    void SetParticleName  (G4String n)        { fParticleName = n; }
    void SetPosition      (G4ThreeVector p)   { fPos = p; }
    void SetMomDir        (G4ThreeVector d)   { fMomDir = d; }
    void SetKineticEnergy (G4double e)        { fKineticEnergy = e; }
    void SetTime          (G4double t)        { fTime = t; }
    void AddEdep          (G4double e)        { fEdep += e; }
    void SetStripID       (G4int id)          { fStripID = id; }   // strip hodoscope only (planeID 5/6)

    G4int         GetPlaneID()       const { return fPlaneID; }
    G4int         GetTrackID()       const { return fTrackID; }
    G4String      GetParticleName()  const { return fParticleName; }
    G4ThreeVector GetPosition()      const { return fPos; }
    G4ThreeVector GetMomDir()        const { return fMomDir; }
    G4double      GetKineticEnergy() const { return fKineticEnergy; }
    G4double      GetTime()          const { return fTime; }
    G4double      GetEdep()          const { return fEdep; }
    G4int         GetStripID()       const { return fStripID; }

private:
    G4int         fPlaneID       = -1;
    G4int         fStripID       = -1;
    G4int         fTrackID       = -1;
    G4String      fParticleName;
    G4ThreeVector fPos;
    G4ThreeVector fMomDir;
    G4double      fKineticEnergy = 0.;
    G4double      fTime          = 0.;
    G4double      fEdep          = 0.;
};

using TrackerHitsCollection = G4THitsCollection<TrackerHit>;

extern G4ThreadLocal G4Allocator<TrackerHit>* gTrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
    if (!gTrackerHitAllocator)
        gTrackerHitAllocator = new G4Allocator<TrackerHit>;
    return gTrackerHitAllocator->MallocSingle();
}

inline void TrackerHit::operator delete(void* hit)
{
    gTrackerHitAllocator->FreeSingle(static_cast<TrackerHit*>(hit));
}
