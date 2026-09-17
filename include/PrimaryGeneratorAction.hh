#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    // Set by MAIN.cc in vis mode; only changes whether RunAction writes CSV output.
    static bool sVisMode;

    // 1 fired primary in N is written to beam_sample_{shadow,open}.csv so the sampled
    // beam can be plotted (pitch angle, gyrophase, energy, source position). Read once
    // from env BEAM_DUMP_PRESCALE; 0 disables. Default 100000 gives ~20k rows per 2e9
    // events, plenty for the distributions and negligible for I/O.
    static long BeamDumpPrescale();

private:
    void DumpBeamSample(G4double xs, G4double ys, G4double zs,
                        G4double E_MeV, const G4ThreeVector& dir) const;

    void     LoadSpectrum(const G4String& filename);
    void     LoadDirectionalSpectrum(const G4String& filename, G4double pitchDeg);
    void     LoadPAD(const G4String& filename);        // build the pitch-angle sampling CDF
    G4double SampleEnergy() const;
    G4double SamplePitch() const;                      // draw a pitch angle [deg] from the PAD
    void     GenerateInvertedBeam(G4Event* event);     // legacy aiming beam (PROD_BEAM=vr)

    // 2026-09-15: wall-thickness investigation (env WALL_TEST=1). Fires the TRUE,
    // un-redirected AP9 PAD directly at the v3 tube collimator's outer (+x) surface, from
    // a source plane positioned just outside it spanning the collimator's full z-length --
    // no WouldHitUranium redirect, no uranium, no tight-world-truncation artifact (the
    // normal imaging beam structurally cannot produce genuinely wide-angle hits on the
    // wall; see CLAUDE.md). Validates the graded (front-thinned) wall design in
    // DetectorConstruction.cc's CollTube namespace against real background, not a guess.
    void GenerateWallTestBeam(G4Event* event);

    // Collimator stand-in: no shielding geometry is built, this just imposes the narrow
    // angular acceptance a real boresight collimator would give. Samples direction
    // uniformly in solid angle within halfAngleRad of +z -- a small-FOV instrument's
    // pitch is pinned near 90 deg (dir.y() ~ 0), so this bypasses the full PAD rather
    // than rejection-sampling it (which would be astronomically inefficient at ~1 deg).
    G4ThreeVector SampleNarrowConeDir(G4double halfAngleRad) const;

    // Straight-line ray-vs-uranium-box pre-check (bending over the ~5 mm source-to-slab
    // gap is negligible). A sampling-efficiency/classification heuristic only -- Geant4
    // tracks the exact helical trajectory and physics regardless of this result.
    bool WouldHitUranium(const G4ThreeVector& pos, const G4ThreeVector& dir) const;

    G4ParticleGun*        fGun = nullptr;
    std::vector<G4double> fEnergies;
    std::vector<G4double> fCDF;
    G4double              fMeanDxBend_mm = 0.;  // flux-weighted mean bending at SSD0 [mm]
    G4double              fZdet          = 0.;  // SSD0 z in G4 world [G4 units]
    G4double              fZsrc          = 0.;  // source plane z [G4 units]
    G4double              fSrcHX         = 0.;  // source half-width [G4 units]
    G4double              fUraniumZ      = 0.;  // uranium slab centre z [G4 units]
    G4double              fUraniumHalfX  = 0.;  // uranium slab half-extents [G4 units]
    G4double              fUraniumHalfY  = 0.;
    G4double              fUraniumHalfZ  = 0.;
    G4double              fPitchLOS_deg  = 0.;  // line-of-sight pitch angle used for J(alpha,E) [deg]
    std::vector<G4double> fPitchGrid;           // fine pitch-angle grid [deg] for PAD sampling
    std::vector<G4double> fPitchCDF;            // cumulative PAD over fPitchGrid, normalized
};
