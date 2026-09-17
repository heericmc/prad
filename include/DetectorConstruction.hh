#pragma once
#include "G4VUserDetectorConstruction.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    enum class Mode { kImaging, kImagingOpen, kCollimatorConcept };

    // 2026-09-14: replaces the 25 cm PE range filter + single scintillator slab with an
    // X-Y strip hodoscope (two crossed layers of 1 cm x 1 cm scintillator strips) --
    // see the StripHodo namespace in DetectorConstruction.cc for the geometry, and
    // TrackerSD/EventAction for the per-strip readout (stripx_id/stripy_id +
    // stripx_edep_MeV/stripy_edep_MeV columns replace poly_edep_MeV/scint_hit in the CSV
    // when this is true). A compile-time toggle, not a runtime env var, like every other
    // geometry variant in this codebase -- exposed here (not just in the .cc) so
    // EventAction.cc and RunAction.cc can branch their CSV logic on it too.
    static constexpr bool kUseStripHodoscope = true;

    // 2026-09-15: exactly 2 Si tracking layers (not 3) in this dedicated repo -- see
    // Det::kTwoLayerViz's comment in DetectorConstruction.cc. Exposed here (not just in
    // the .cc's Det:: namespace) so EventAction.cc can drop the layer-2 coincidence
    // requirement and CSV columns, which would otherwise make every event fail the
    // coincidence gate (there is no planeID-2 sensitive detector in this geometry).
    static constexpr bool kTwoLayerViz = true;

    DetectorConstruction()  = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    static void SetMode(Mode m) { sMode = m; }
    static Mode GetMode()       { return sMode; }

    static G4double sDetOffsetX_mm;   // flux-weighted x offset of detector centre [mm]
    static G4double GetDetOffsetX_mm() { return sDetOffsetX_mm; }

    static G4double sStandoff_mm;     // uranium-SSD0 standoff distance [mm]
    static G4double GetStandoff_mm()  { return sStandoff_mm; }

    static G4double sSrcHX_mm;        // source illumination half-width [mm], auto-scaled with standoff
    static G4double GetSrcHX_mm()     { return sSrcHX_mm; }

    static G4double sActiveHalfX_mm;  // tracker active half-width x [mm] (Det::kActiveHalfX)
    static G4double sActiveHalfY_mm;  // tracker active half-height y [mm] (Det::kActiveHalfY)
    static G4double GetActiveHalfX_mm() { return sActiveHalfX_mm; }
    static G4double GetActiveHalfY_mm() { return sActiveHalfY_mm; }

    // Uranium slab half-extents [mm] (UraniumDet::kHalfX/Y/Z), exposed so the primary
    // generator can test whether a candidate trajectory would hit it. z-centre is
    // derivable from the standoff alone (-0.5*kStandoff) so isn't duplicated here.
    static G4double sUraniumHalfX_mm;
    static G4double sUraniumHalfY_mm;
    static G4double sUraniumHalfZ_mm;
    static G4double GetUraniumHalfX_mm() { return sUraniumHalfX_mm; }
    static G4double GetUraniumHalfY_mm() { return sUraniumHalfY_mm; }
    static G4double GetUraniumHalfZ_mm() { return sUraniumHalfZ_mm; }

    // 2026-09-15: v3 tube collimator's full z-span + outer half-width [mm], exposed so
    // PrimaryGeneratorAction's wall-test beam (WALL_TEST env var) can position a source
    // plane just outside it, to validate the graded (front-thinned) wall design against
    // the true un-redirected wide-angle flux. Only meaningful when Det::kCollStyle==kTube.
    static G4double sCollZFront_mm;
    static G4double sCollZBack_mm;
    static G4double sCollOuterHalfXY_mm;
    static G4double GetCollZFront_mm()      { return sCollZFront_mm; }
    static G4double GetCollZBack_mm()       { return sCollZBack_mm; }
    static G4double GetCollOuterHalfXY_mm() { return sCollOuterHalfXY_mm; }

private:
    static Mode sMode;

    G4LogicalVolume* fSSD0_LV  = nullptr;   // tracking layer 0 (Si sensor)
    G4LogicalVolume* fSSD1_LV  = nullptr;   // tracking layer 1
    G4LogicalVolume* fSSD2_LV  = nullptr;   // tracking layer 2
    G4LogicalVolume* fPoly_LV  = nullptr;   // unused when kUseStripHodoscope
    G4LogicalVolume* fScint_LV = nullptr;   // unused when kUseStripHodoscope
    G4LogicalVolume* fStripX_LV = nullptr;  // measures x: strips segmented along x, each spans full y
    G4LogicalVolume* fStripY_LV = nullptr;  // measures y: strips segmented along y, each spans full x
    // 2026-09-15: multi-zone taper (mass-optimized, see CollTube::kZones), one LV per
    // zone, all sensitive -- see EventAction.cc wall_hits.csv.
    std::vector<G4LogicalVolume*> fCollTubeWallZones_LV;
};
