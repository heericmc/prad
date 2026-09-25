// Geometry for a passive standoff nuclear-material detection sensor (orbital treaty
// verification), NOT a weapon/implosion simulation — see README.md. The "uranium slab"
// below is a generic dense stand-in target whose shadow in the ambient AP9 trapped-
// proton flux the tracker is meant to detect; no weapon-internals or explosive physics
// is modeled here or anywhere else in this repo.

#include "DetectorConstruction.hh"
#include "TrackerSD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <string>

DetectorConstruction::Mode    DetectorConstruction::sMode          = DetectorConstruction::Mode::kImaging;
G4double                      DetectorConstruction::sDetOffsetX_mm = 0.0;
G4double                      DetectorConstruction::sStandoff_mm   = 0.0;
G4double                      DetectorConstruction::sSrcHX_mm      = 0.0;
G4double                      DetectorConstruction::sActiveHalfX_mm = 0.0;
G4double                      DetectorConstruction::sActiveHalfY_mm = 0.0;
G4double                      DetectorConstruction::sUraniumHalfX_mm = 0.0;
G4double                      DetectorConstruction::sUraniumHalfY_mm = 0.0;
G4double                      DetectorConstruction::sUraniumHalfZ_mm = 0.0;
G4double                      DetectorConstruction::sCollZFront_mm = 0.0;
G4double                      DetectorConstruction::sCollZBack_mm  = 0.0;
G4double                      DetectorConstruction::sCollOuterHalfXY_mm = 0.0;

// Uranium slab dimensions
namespace UraniumDet {
    constexpr G4double kHalfX = 20.0 * cm;   // 40 cm wide
    constexpr G4double kHalfY = 20.0 * cm;   // 40 cm tall
    constexpr G4double kHalfZ =  2.5 * cm;   // 5 cm thick
}

// Tracker layer dimensions (3-layer ALTAI MAPS design).
// Each layer, beam-facing to downstream: kapton | copper | Si sensor | K13D2U cold
// plate. Chips are 15x30 mm, 50 um Si; a stave is a 2x5 chip array = 30x150 mm of
// silicon sitting on the leg-free centre of a wider 34x150 mm K13D2U cold plate --
// the extra 2 mm on each side is a 6.6 mm leg. Staves tile in x with the cold plates
// touching, so adjacent legs meet and leave a 4 mm silicon-to-silicon gap per stave.
// In y the rows are separated by a 10 mm gap of bare (continuous) cold plate, so the
// y pitch is 160 mm = 150 mm silicon + 10 mm gap.
namespace Det {
    constexpr G4double kSensorHalfX = 1.5     * cm;    // 30 mm silicon (2 chips x 15 mm)
    constexpr G4double kColdHalfX   = 1.7     * cm;    // 34 mm cold plate = sensor + 2 mm leg/side
    constexpr G4double kLegHalfX    = 0.1     * cm;    // 2 mm leg width (one per cold-plate edge)
    constexpr G4double kStavePitchX = 3.4     * cm;    // = cold-plate width (staves touch in x)
    constexpr G4int    kNStaveX     = 31;              // staves in x (odd -> one centred at 0)
    constexpr G4double kActiveHalfX = kNStaveX * kStavePitchX / 2.;  // 527 mm cold-plate half-extent
    constexpr G4double kSensorHalfY =  7.5    * cm;    // 150 mm silicon per stave row (5 chips x 30 mm)
    constexpr G4double kRowGapY     =  1.0    * cm;    // 10 mm carbon-fibre gap between stave rows
    constexpr G4double kStavePitchY = 2. * kSensorHalfY + kRowGapY;  // 160 mm row pitch
    constexpr G4int    kNStaveY     = 7;               // stave rows in y (odd -> one centred at 0)
    constexpr G4double kActiveHalfY = kNStaveY * kStavePitchY / 2.;  // 560 mm cold-plate half-extent
    // 31 x 7 staves = 9765 cm^2 of silicon. Rows are not continuous in y -- each pair is
    // separated by kRowGapY of bare cold plate. Fill factor 0.882 (x) x 0.938 (y) = 0.827.
    // For the smaller 5-stave build: kNStaveX = 5, kNStaveY = 1, kRowGapY = 0.
    constexpr G4double kSiHalfZ     =  0.0025 * cm;    // 50 um silicon sensor
    constexpr G4double kKaptonHalfZ =  0.00375* cm;    // 75 um kapton (top of sensor)
    constexpr G4double kCuHalfZ     =  0.0009 * cm;    // 18 um copper
    constexpr G4double kColdHalfZ   =  0.02   * cm;    // 400 um K13D2U cold plate
    constexpr G4double kColdDensity_gcm3 = 1.8;        // K13D2U density (logged for archive_run.ps1)
    constexpr G4double kLegHalfZ    =  0.33   * cm;    // 6.6 mm leg height (along beam)
    constexpr G4double kLayerGapZ   =  2.0    * cm;    // 20 mm layer-to-layer spacing
    // 2026-09-15: exactly 2 Si tracking layers (not 3), kTwoLayerGapZ apart (not
    // kLayerGapZ). Originally added purely for a one-off VRML export/render (see the
    // main proton_tracker_prototype repo, where this flag stays permanently FALSE --
    // every physics number/production run there is based on the real 3-layer/20 mm
    // design). THIS repo (proton_tracker_2layer10cm) is a fully independent copy of
    // src/include/mac_files/CMakeLists.txt + root data files, created specifically to
    // run a real 2-layer/10 cm + 10-deg tube-taper collimator simulation campaign
    // (100/200/400 m) without sharing any mutable source state with the main repo's
    // live production campaign. kTwoLayerViz is permanently TRUE here for exactly that
    // reason -- do not copy this change back into the main repo.
    constexpr bool     kTwoLayerViz   = DetectorConstruction::kTwoLayerViz;  // see the
                                        // header's copy of this flag for why it's exposed
                                        // there too (EventAction.cc needs it)
    constexpr G4double kTwoLayerGapZ  = 10.0   * cm;   // 100 mm, matching the old
                                                        // standalone CollimatorDet
                                                        // concept's spacing rationale
    constexpr G4double kPoly_HalfZ  = 12.5    * cm;    // 25 cm PE range filter
    constexpr G4double kScint_HalfZ =  1.0    * cm;    // 2 cm scintillator
    // v3 puts the range filter between layers 1 and 2 instead of behind all three.
    // Set false for the v2 order. Layers 0-1 keep their 20 mm spacing either way, so
    // geo_theta and the angular budget are unaffected.
    constexpr bool     kFilterBeforeLayer2 = false;
    // "v3" redefined 2026-09-10: a real HDPE lattice collimator (CollV3 namespace, below)
    // mounted immediately upstream of layer 0 -- physical shielding, not just the
    // PrimaryGeneratorAction source-term emulation. Independent of kFilterBeforeLayer2
    // (that flag never had an archived run taken with it, so no data collides with this
    // redefinition). GeometryTag below tags this "_v3coll" when true.
    constexpr bool     kHasCollimator = true;
    // Which collimator geometry to build when kHasCollimator is true. kLattice = CollV3's
    // periodic square-hole HDPE (or Al) lattice, ~1 m deep, uniform per-position angular
    // acceptance. kTube = CollTube's single large square HDPE tube around the tracker's
    // whole perimeter (2026-09-14), ~45 m deep, ON-AXIS acceptance only -- see the CollTube
    // namespace comment for why it can't match the lattice's uniform per-pixel guarantee.
    enum class CollStyle { kLattice, kTube };
    constexpr CollStyle kCollStyle = CollStyle::kTube;
}

// v3 collimator: HDPE square-hole lattice immediately upstream of layer 0, mounted on the
// same tilted/offset stack (shares detRot/detOffset_x in Construct()) so its bore points
// along the mean bent trajectory -- real shielding physics (absorption, wall-scattering),
// not the PrimaryGeneratorAction source-term emulation this project used before. Design
// point: literally collimator.md's "Fixed at 1 deg" table, 49% open area row (wall =
// 25cm*sin(1deg) = 4.36 mm, hole solved for 49% open area = 10.18 mm, depth ~107 cm).
// 2026-09-11 CORRECTION: this was first built "re-solved at 2 deg" instead of using the
// 1 deg row directly -- wrong; the collimator itself is the 1 deg design, PGA's
// COLLIMATOR_DEG=2 is a separate, deliberately-wider SAMPLING margin so the source term
// doesn't clip real physics near the true 1 deg boundary (see PrimaryGeneratorAction.cc
// and imaging_analysis.m's COLLIMATOR_GATE_DEG, which must use this 1 deg design angle,
// not the 2 deg sampling one, for its +3sigma gate). Every v3 archive taken before this
// fix (100 m - 1600 m) used the wrong 2 deg collimator and needs re-running. Lattice sized
// to fully cover the tracker's rectangular active area (+-527 x +-560 mm) -- kNHolesX/Y
// differ because the footprint isn't square, unlike the CollimatorDet visualization rig.
namespace CollV3 {
    // Material comparison, 2026-09-14: same lattice geometry (hole/wall/depth/open-area,
    // all purely geometric per the design formula below -- none of it depends on material),
    // swapped HDPE -> aluminum to see how a weaker-stopping-power wall changes wall-clip
    // scattering/leakage. Same treatment as the earlier U-vs-Al uranium-slab comparison:
    // dimensions unchanged, only material differs. GeometryTag/log both tag "Al" so this
    // can't silently land in the existing HDPE run_*_v3_coll_1deg/ archives.
    constexpr bool     kAlCollimator = true;  // true = G4_Al, false = G4_POLYETHYLENE (HDPE)
    constexpr G4double kDesignDeg = 1.0;      // acceptance half-angle design point (collimator.md)
    constexpr G4double kOpenArea  = 0.49;     // fraction, matches collimator.md's point design
    const G4double kAngleRad   = kDesignDeg * deg;
    const G4double kWallThick  = 25.0 * cm * std::sin(kAngleRad);
    const G4double kHoleWidth  = kWallThick * std::sqrt(kOpenArea) / (1.0 - std::sqrt(kOpenArea));
    const G4double kHoleHalfXY = kHoleWidth / 2.0;
    const G4double kPitch      = kHoleWidth + kWallThick;
    const G4double kCollHalfZ  = (std::sqrt(2.0) * kHoleWidth + kWallThick)
                                / (2.0 * std::tan(kAngleRad));
    const G4int    kNHolesX    = 2 * G4int(std::ceil(Det::kActiveHalfX / kPitch)) + 1;  // odd, centred
    const G4int    kNHolesY    = 2 * G4int(std::ceil(Det::kActiveHalfY / kPitch)) + 1;
    // Same "+wall/2" perimeter-rim trick as CollimatorDet (see its comment): N*pitch/2
    // alone leaves only a half-thickness rim at the outer edge.
    const G4double kCollHalfX  = kNHolesX * kPitch / 2.0 + kWallThick / 2.0;
    const G4double kCollHalfY  = kNHolesY * kPitch / 2.0 + kWallThick / 2.0;
    constexpr G4double kCollGap = 5.0 * mm;    // clearance, collimator to layer 0
}

// v3 "tube" collimator, 2026-09-14: a single large square HDPE tube wrapping the
// PERIMETER of the tracker's whole active area, replacing CollV3's periodic lattice of
// many small holes with one big bore. Selected via Det::kCollStyle = kTube.
//
// Wall thickness -- literal "must stop protons below 200 MeV" requirement: reuses the
// SAME verified 25 cm HDPE range already established for Det::kPoly_HalfZ's downstream
// range filter (measured threshold ~195 MeV: a proton needs roughly >=200 MeV to punch
// through 25 cm of G4_POLYETHYLENE at normal incidence). This is deliberately NOT
// CollV3::kWallThick's 25cm*sin(angle) form -- that formula sizes a thin lattice SEPTUM
// against a proton grazing nearly ALONG its length at the design angle, so the septum's
// own face-on thickness can be much less than 25 cm and still present a 25 cm path length
// to a grazing ray. This tube's wall can instead be struck face-on by a background proton
// from any direction (there's no adjacent hole forcing near-axial incidence), so its own
// transverse thickness must BE the full 25 cm stopping range.
//
// Length -- a single tube cannot reproduce the lattice's UNIFORM per-position acceptance:
// with no internal wall structure, nothing clips an off-centre ray before it reaches the
// far wall, so kDesignDeg here is only the ON-AXIS (boresight-centred) field of view, not
// a per-pixel guarantee like CollV3. Defined by the worst-case ray from the CENTRE of the
// exit face to the farthest CORNER of the entrance aperture (diagonal sqrt(2)*kBoreHalfXY,
// since the bore is a 2D square, not a 1D slit):
//     tan(kDesignDeg) = sqrt(2) * kBoreHalfXY / length
// A hit away from centre sees a WIDER effective angle (up to ~2x at the extreme opposite
// corner, where both entry and exit can be at opposite corners) -- an inherent cost of
// dropping the lattice's internal walls, not a tunable parameter. At the 1 deg design
// point and a bore sized to the full +-527 x +-560 mm tracker footprint, this comes out to
// tens of metres -- "really long" is the expected, correct answer, not a bug: it is the
// quantitative price of collimating a large-area detector to a fine angle with one big
// aperture instead of a fine lattice.
namespace CollTube {
    // 2026-09-14: widened 1 -> 10 deg for a short-bore comparison sweep (100-400 m).
    // Length scales as 1/tan(angle), so this shrinks the bore from ~45.4 m to ~4.5 m --
    // NB PrimaryGeneratorAction's COLLIMATOR_DEG (source-sampling cone) must be widened
    // to match (>=10 deg), or the open/shadow beams would keep sampling only the old
    // narrow 1 deg cone and never populate most of this wider bore's true acceptance.
    //
    // 2026-09-16: narrowed 10 -> 5 deg to bring the whole-detector pile-up rate down to
    // the maps-pileup-collimation.md Sec.9 target (~15 hits/us at 5 deg, vs ~62 at 10 deg,
    // vs a 10 hits/us budget) -- see that document's independently-AP9-CSV-derived table.
    // Length scales as 1/tan(angle) so this roughly DOUBLES the bore vs the 10 deg design
    // (~4.6 m -> ~9.2 m); NB PrimaryGeneratorAction's COLLIMATOR_DEG source-sampling cone
    // must be narrowed to match for a real production run, or the beam keeps sampling the
    // old wider cone and wastes most of its statistics outside this bore's true acceptance.
    //
    // 2026-09-16 (same day, later): re-derived from a principled bound instead of a round
    // number. The two-plane track-confusion probability pn(alpha) = X_unshielded *
    // f(alpha) * tau * pi * alpha_rad^2 * d^2 (maps-pileup-collimation.md Sec.2/3) was
    // independently solved (three separate numerical methods -- fzero, a 19001-pt grid
    // scan, and manual bisection -- all agreeing to 4 decimal places) for the angle where
    // pn crosses 1e-3, the ceiling below which contamination stays under this analysis's
    // statistical noise floor (~1/sqrt(N) for realistic trigger counts) regardless of
    // exactly how large N ends up being: alpha = 8.01 deg, pn(5 deg) = 1.5e-4 by
    // comparison. Widening 5 -> 8 deg raises the whole-detector rate (Sec.9 of that
    // document) from ~15.7 to ~34 hits/us -- still using the same idealized flat-aperture
    // estimate documented there as a conservative upper bound, not a measurement.
    // SUPERSEDED 2026-09-17 -- see the next note: the pn(alpha)=8.01deg solve above used
    // d = 2 cm (the ORIGINAL 3-layer/20mm project's layer pitch), not this repo's actual
    // d = Det::kTwoLayerGapZ = 10 cm. Left in place only as a record of the (invalid)
    // reasoning that originally picked 8 deg.
    //
    // 2026-09-17: re-derived correctly with d = 10 cm and, more importantly, with the
    // ceiling now binding on the REAL 3-point confirmation this design actually enforces
    // -- both Si planes AND a position-matched hit in the StripHodo trigger stage (not
    // just "some strip fired": Analysis must check the fired stripx_id/stripy_id cell is
    // consistent with the extrapolated plane0->plane1 track; see maps-pileup-collimation.md
    // Sec.3a for the CSV columns this needs). Two independent probabilities multiply,
    // since the tracker-plane and hodoscope accidental hits come from separate,
    // uncorrelated background protons:
    //   pn_2plane(alpha) = X_unshielded * f(alpha) * tau * pi * alpha_rad^2 * d^2   (d=10cm)
    //   p_strip(alpha)   = X_unshielded * f(alpha) * tau * A_cell                   (A_cell=1cm^2)
    //   pn_3point(alpha) = pn_2plane(alpha) * p_strip(alpha)
    // Solving pn_3point(alpha) = 1e-3 (three independent methods again -- direct
    // bisection, a narrower-bracket bisection, and a closed-form alpha^4-anchored fit --
    // agreeing to 0.01 deg) gives alpha = 11.72 deg. At the OLD 8 deg point, pn_3point was
    // already comfortably under the ceiling (1.0e-4, 10x margin) -- this widening spends
    // that margin on a shorter, lighter-per-unit-length (though not lighter overall, see
    // the ZoneSpec comment below) bore instead of leaving it unused. Bare two-plane
    // confusion at 11.72 deg is 0.115 (11.5%) -- far over the two-plane-only ceiling on
    // its own -- so the strip hodoscope's position match is now LOAD-BEARING, not optional
    // insurance: this design's confusion-ceiling claim is only true if the analysis
    // actually enforces that position match, not merely a "hodoscope fired" flag.
    // Length scales as 1/tan(angle), so this shrinks the front-section bore from ~5.66 m
    // (8 deg) to ~3.83 m (11.72 deg); the ZoneSpec taper below was re-segmented to match
    // (see that comment) -- **NOT YET RE-VALIDATED with a fresh WALL_TEST run** (the
    // zero-sub-200-MeV-leakage claim for the previous 8 deg/5-zone taper doesn't
    // automatically carry over to new zone boundaries; re-run WALL_TEST=1 after building
    // this before trusting it for a production campaign).
    constexpr G4double kDesignDeg  = 11.7;
    // 2026-09-14: HDPE -> aluminum, to use less material. NIST PSTAR: 12.4 cm of Al stops
    // protons <=200 MeV (vs 25 cm of HDPE for the same cutoff) -- Al is ~2x denser (2.70 vs
    // ~0.94-0.97 g/cm^3) so it wins on LINEAR range despite a somewhat worse per-gram
    // stopping power (HDPE's hydrogen content boosts electron density/gram). Flat/uniform
    // thickness, not angle-tapered: a tapered (sin(angle)-scaled, ~2.15 cm) wall was
    // considered but rejected for now -- it only helps protons grazing near the design-
    // angle boundary (same physics as CollV3's septum), giving zero benefit against a
    // genuinely wide-angle hit (sin(80 deg) ~ 1), and this tube (unlike the lattice) has no
    // periodic structure to geometrically rule those out. Revisit only after measuring the
    // actual wide-angle leakage rate in Geant4.
    constexpr bool     kAlWall     = true;
    constexpr G4double kWallThick  = kAlWall ? 124.0 * mm : 250.0 * mm;   // full (100%) thickness
    // 2026-09-15: multi-zone mass-optimized taper -- see the ZoneSpec table where the tube
    // is actually built (Construct(), Det::CollStyle::kTube branch) for the validated
    // per-zone thickness fractions and the wall-test-beam investigation behind them.
    // kWallThick above remains the reference (100%) value every zone's fraction scales.
    const G4double kAngleRad    = kDesignDeg * deg;
    // Bore: square, sized to fully contain the tracker's rectangular active area
    // (+-527 x +-560 mm) -- use the larger half-extent, +2 mm clearance. The margin
    // matters now that the bore also runs alongside the tracker itself (side shielding,
    // 2026-09-14): the cold-plate sheets are placed at EXACTLY +-560 mm in y (Det::
    // kActiveHalfY), so a bore sized to that exact number would touch them with zero
    // clearance -- the same failure mode tiltClear exists elsewhere in this file to avoid.
    const G4double kBoreHalfXY  = std::max(Det::kActiveHalfX, Det::kActiveHalfY) + 2.0 * mm;
    const G4double kOuterHalfXY = kBoreHalfXY + kWallThick;   // full-thickness outer (world sizing etc.)
    // kCollHalfZ is the ACCEPTANCE-ANGLE-derived half-depth only (front section, upstream
    // of layer 0) -- the actual built structure is longer: Construct() extends it downstream
    // past the whole detector stack (layer 0 through the trigger stage) as perimeter side
    // shielding, 2026-09-14, so background can't sneak in past the sides of the tracker/
    // hodoscope once it's already inside the tube. Same wall thickness/material throughout
    // (no angle benefit applies to side shielding -- see kAlWall comment above).
    const G4double kCollHalfZ   = (std::sqrt(2.0) * kBoreHalfXY) / (2.0 * std::tan(kAngleRad));
    constexpr G4double kCollGap = 5.0 * mm;        // clearance, collimator to layer 0
    constexpr G4double kSideMargin = 10.0 * mm;    // clearance past the last detector volume
}

// X-Y strip hodoscope, 2026-09-14: replaces the 25 cm PE range filter + single 2 cm
// scintillator slab with two crossed layers of 1 cm x 1 cm x 1 cm plastic scintillator
// strips, selected via DetectorConstruction::kUseStripHodoscope. Purpose: coarse (1 cm^2) hit position
// PLUS an energy-deposit (dE/dx) measurement per event, instead of the old poly filter's
// binary >200 MeV pass/fail cut. Removing the poly filter is deliberate -- sub-200 MeV
// protons now reach the strips too (some may even range out and stop inside the 1 cm
// plastic, Bragg-peaking), so the deposited energy itself becomes the energy proxy the
// post-analysis inverts, rather than a hard threshold.
//
// Layer naming is by what each layer MEASURES, not the strips' own long-axis direction:
// "StripX" strips are segmented ALONG x (1 cm pitch in x) and each spans the FULL y
// active extent, so which strip fired gives x-position. "StripY" is the transpose (
// segmented along y, spans full x), giving y-position. Together the two crossed layers
// give a 2D (x,y) hit cell -- see EventAction.cc / imaging_analysis.m for how the strip
// indices map back to a physical 1 cm^2 spot.
//
// No lattice/tube-style design-angle derivation here: strip width and thickness are
// fixed (1 cm each) design inputs, not solved from an acceptance angle.
namespace StripHodo {
    constexpr G4double kPitch   = 1.0 * cm;   // strip width AND thickness
    constexpr G4double kHalfZ   = kPitch / 2.0;
    // number of strips needed to span the tracker's full rectangular active area
    const G4int kNStripsX = G4int(std::ceil(2.0 * Det::kActiveHalfX / kPitch));  // measures x
    const G4int kNStripsY = G4int(std::ceil(2.0 * Det::kActiveHalfY / kPitch));  // measures y
    constexpr G4double kGap = 5.0 * mm;   // clearance: layer-2 to StripX, StripX to StripY
}

// Collimator concept: a stand-alone visualization geometry, no uranium, no beam
// physics, no field. Two 1 m^2 silicon planes 10 cm apart (centre to centre), with a
// square-hole HDPE collimator immediately upstream of the first plane.
//
// 10 cm separation: sigma_pos = sqrt(2)*(pixel pitch/sqrt(12))/L falls to ~10% of
// sigma_MCS (1.305 mrad, negligible in quadrature, <1% effect on sigma_tot) by
// L ~ 94 mm. 100 mm rounds that up with a little headroom, well short of the ~1 m a
// long-lever-arm design would need, since the pixel term is already subdominant to MCS
// at this spacing and a much longer arm would only buy mass/volume/acceptance cost for
// a shrinking return (see CLAUDE.md).
//
// Collimator: 49%-open-area point design from collimator.md's "Fixed at 1 deg" table
// (angle pinned at 1 deg rather than solved from the rate budget, so wall thickness is
// just 25 cm x sin(1 deg) regardless of open area; only hole width and depth vary).
// Table cells round to 4 mm wall / 10 mm hole / 107 cm depth; this uses the unrounded
// values (wall = 25cm*sin(1deg) = 4.363 mm, hole solved for exactly 49% open area =
// 10.181 mm) so the logged open area actually reads 49%, not the ~51% four-and-ten
// whole-mm rounding would give -- depth still comes out to the table's 107 cm (1074.6 mm
// exactly). Modelled as a solid HDPE block with the holes cut out via a single
// G4SubtractionSolid against a G4MultiUnion of the hole boxes (see below) rather than a
// grid of vacuum daughter volumes, since a daughter's invisible boundary wouldn't show
// up as a real perforation in an exported view.
namespace CollimatorDet {
    constexpr G4double kPlaneHalfXY = 500.0  * mm;   // 1 m x 1 m plane
    constexpr G4double kPlaneHalfZ  = 0.0025 * cm;   // 50 um Si (matches Det::kSiHalfZ)
    constexpr G4double kPlaneSep    = 100.0  * mm;   // 10 cm, plane-centre to plane-centre

    constexpr G4double kHoleHalfXY  = 5.09   * mm;   // 10.18 mm hole width / 2
    constexpr G4double kWallThick   = 4.36   * mm;   // 25 cm x sin(1 deg)
    constexpr G4double kPitch       = 2. * kHoleHalfXY + kWallThick;   // 14.54 mm
    constexpr G4int    kNHoles      = 69;            // 69 x 14.54 mm = 1003.5 mm footprint
    // N*pitch/2 alone puts the block's outer face at the midpoint of what would be the
    // next wall, leaving only a half-thickness rim (wall/2) between the outermost hole
    // and the edge -- visibly thinner than every interior wall. +wall/2 grows the block
    // just enough to give the perimeter a full wall thickness, matching the interior.
    constexpr G4double kCollHalfXY  = kNHoles * kPitch / 2. + kWallThick / 2.;  // 504.0 mm
    constexpr G4double kCollHalfZ   = 537.3  * mm;   // 1074.6 mm depth / 2 (~107 cm)
    constexpr G4double kCollGap     = 5.0    * mm;   // clearance, collimator to plane 0

    // Filter + scintillator downstream of plane 1, same spec as Det:: (25 cm PE range
    // filter, 2 cm plastic scintillator) -- the trigger/coincidence stack the analysis
    // needs behind the two tracking planes, reusing Det::kPoly_HalfZ / kScint_HalfZ.
    constexpr G4double kFilterGap   = 5.0    * mm;   // clearance, plane 1 to filter
    constexpr G4double kScintGap    = 5.0    * mm;   // clearance, filter to scintillator
}

// Standoff imaging world. kStandoff is the single variable to change for a different
// standoff distance -- everything else (uranium z, SSD positions, world size) derives
// from it. B = 0.1 G in +y (equatorial dipole at L ~ 1.4: B_eq = 0.31/L^3 ~ 0.11 G).
namespace ImagingDet {
    constexpr G4double kStandoff  = 400.0 * m;   // change for a different standoff
    // Floor on the source illumination half-width. fSrcHX auto-scales with standoff in
    // Construct() below so the back-projected open field stays flat across the analysis
    // cell, but never drops under this. PrimaryGeneratorAction reads the scaled value
    // via GetSrcHX_mm() and the world is sized from it, so the two can't drift apart.
    constexpr G4double kSrcHXFloor = 400.0 * mm;

    constexpr G4double kUraniumZ  = -0.5 * kStandoff;
    constexpr G4double kSSD0_Z    =  0.5 * kStandoff;
    constexpr G4double kWorldHZ   =  0.5 * kStandoff + 2.5 * m;
    // kWorldHX/HY is computed at runtime from the max bending at E_min = 200 MeV.

    constexpr G4double kBField    =   0.1 * gauss;  // 1e-5 T
}

// Flux-weighted mean lateral bending and detector alignment angle. Reads
// diff_flux_AP9_i316.csv, computes sum(J/r_g)/sum(J) (the flux-weighted mean inverse
// gyroradius), and returns the corresponding displacement and tilt angle.
// Returns { dx_mm (negative, toward -x), theta_rad (positive) }.
namespace {
std::pair<G4double,G4double>
ComputeAlignmentBending(G4double B_T, G4double d_m)
{
    const G4double mp = 938.272;  // proton rest mass [MeV/c^2]

    G4double sumJ = 0., sumJ_inv_rg = 0.;

    std::ifstream f("diff_flux_AP9_i316.csv");
    if (f.is_open()) {
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream ss(line);
            G4double eMeV, j;
            if (!(ss >> eMeV >> j) || eMeV < 200. || j <= 0.) continue;
            const G4double p   = std::sqrt(eMeV * (eMeV + 2. * mp));  // MeV/c
            const G4double r_g = p / (299.8 * B_T);                    // m
            sumJ       += j;
            sumJ_inv_rg += j / r_g;
        }
    }

    G4double dx_mm, theta_rad;
    if (sumJ > 0.) {
        const G4double mean_inv_rg = sumJ_inv_rg / sumJ;          // 1/m
        dx_mm    = -0.5 * d_m * d_m * mean_inv_rg * 1000.;        // mm (negative)
        theta_rad =       d_m         * mean_inv_rg;               // rad (positive)
    } else {
        // Fallback: 300 MeV representative proton near the AP9 peak at L = 1.4
        const G4double p0   = std::sqrt(300. * (300. + 2. * mp));
        const G4double r_g0 = p0 / (299.8 * B_T);
        dx_mm    = -0.5 * d_m * d_m / r_g0 * 1000.;
        theta_rad =       d_m          / r_g0;
        G4cout << "[Det] AP9 CSV not found; using 300 MeV fallback for alignment\n";
    }

    G4cout << "[Det] Geomagnetic alignment:"
           << "  B = " << B_T * 1e4 << " G"
           << "  d = " << d_m << " m"
           << "  dx = " << dx_mm << " mm"
           << "  theta = " << theta_rad * 1e3 << " mrad"
           << G4endl;
    return { dx_mm, theta_rad };
}
} // namespace

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();
    auto* vac  = nist->FindOrBuildMaterial("G4_Galactic");
    auto* Si   = nist->FindOrBuildMaterial("G4_Si");
    auto* poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    auto* scintMat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    // Collimator concept: self-contained, no uranium/beam/field machinery below.
    if (sMode == Mode::kCollimatorConcept) {
        const G4double zPlane0  = -CollimatorDet::kPlaneSep / 2.;
        const G4double zPlane1  =  CollimatorDet::kPlaneSep / 2.;
        const G4double zColl    = zPlane0 - CollimatorDet::kCollGap - CollimatorDet::kCollHalfZ;
        const G4double zFilter  = zPlane1 + CollimatorDet::kPlaneHalfZ + CollimatorDet::kFilterGap
                                 + Det::kPoly_HalfZ;
        const G4double zScint   = zFilter + Det::kPoly_HalfZ + CollimatorDet::kScintGap
                                 + Det::kScint_HalfZ;
        const G4double margin   = 50.0 * mm;
        const G4double worldHXY = CollimatorDet::kCollHalfXY + margin;
        const G4double worldHZ  = std::max(std::abs(zColl) + CollimatorDet::kCollHalfZ,
                                           zScint + Det::kScint_HalfZ) + margin;

        auto* worldLV = new G4LogicalVolume(
            new G4Box("World", worldHXY, worldHXY, worldHZ), vac, "WorldLV");
        auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                          nullptr, false, 0, true);
        worldLV->SetVisAttributes(new G4VisAttributes(false));

        // Two 1 m^2 silicon planes.
        auto* siVis = new G4VisAttributes(G4Colour(0.2, 0.5, 1.0, 0.85));
        siVis->SetForceSolid(true);
        auto* planeSolid = new G4Box("Plane", CollimatorDet::kPlaneHalfXY, CollimatorDet::kPlaneHalfXY,
                                     CollimatorDet::kPlaneHalfZ);
        auto* plane0LV = new G4LogicalVolume(planeSolid, Si, "Plane0LV");
        auto* plane1LV = new G4LogicalVolume(planeSolid, Si, "Plane1LV");
        plane0LV->SetVisAttributes(siVis);
        plane1LV->SetVisAttributes(siVis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPlane0), plane0LV, "Plane0PV",
                          worldLV, false, 0, true);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPlane1), plane1LV, "Plane1PV",
                          worldLV, false, 0, true);

        // HDPE collimator: solid block with a 13x13 lattice of full-depth holes actually
        // cut out via Boolean subtraction, so the exported surface shows real
        // perforations (a vacuum-filled daughter volume would stay invisible and the
        // block would render as an unbroken solid).
        auto* hdpe = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
        auto* hdpeVis = new G4VisAttributes(G4Colour(0.9, 0.85, 0.7, 1.0));
        hdpeVis->SetForceSolid(true);

        auto* holeShape = new G4Box("HoleShape", CollimatorDet::kHoleHalfXY, CollimatorDet::kHoleHalfXY,
                                    CollimatorDet::kCollHalfZ + 1.0 * mm);   // overshoot for a clean cut
        auto* holesUnion = new G4MultiUnion("Holes");
        for (G4int i = 0; i < CollimatorDet::kNHoles; ++i) {
            const G4double xc = (i - (CollimatorDet::kNHoles - 1) / 2.) * CollimatorDet::kPitch;
            for (G4int j = 0; j < CollimatorDet::kNHoles; ++j) {
                const G4double yc = (j - (CollimatorDet::kNHoles - 1) / 2.) * CollimatorDet::kPitch;
                G4Transform3D tr(G4RotationMatrix(), G4ThreeVector(xc, yc, 0));
                holesUnion->AddNode(*holeShape, tr);
            }
        }
        holesUnion->Voxelize();

        auto* collBlock = new G4Box("CollimatorBlock", CollimatorDet::kCollHalfXY, CollimatorDet::kCollHalfXY,
                                    CollimatorDet::kCollHalfZ);
        auto* collSolid = new G4SubtractionSolid("Collimator", collBlock, holesUnion);

        auto* collLV = new G4LogicalVolume(collSolid, hdpe, "CollimatorLV");
        collLV->SetVisAttributes(hdpeVis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zColl), collLV, "CollimatorPV",
                          worldLV, false, 0, true);

        // Trigger/coincidence stack behind plane 1: 25 cm PE range filter (same spec as
        // Det::kPoly_HalfZ, threshold ~195-200 MeV) then a 2 cm plastic scintillator
        // (Det::kScint_HalfZ), same materials as the main imaging tracker.
        auto* filterLV = new G4LogicalVolume(
            new G4Box("Filter", CollimatorDet::kPlaneHalfXY, CollimatorDet::kPlaneHalfXY, Det::kPoly_HalfZ),
            poly, "FilterLV");
        auto* filterVis = new G4VisAttributes(G4Colour(0.5, 0.8, 1.0, 0.4));
        filterVis->SetForceSolid(true);
        filterLV->SetVisAttributes(filterVis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zFilter), filterLV, "FilterPV",
                          worldLV, false, 0, true);

        auto* scintLV = new G4LogicalVolume(
            new G4Box("Scint", CollimatorDet::kPlaneHalfXY, CollimatorDet::kPlaneHalfXY, Det::kScint_HalfZ),
            scintMat, "ScintLV");
        auto* scintVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.2, 0.9));
        scintVis->SetForceSolid(true);
        scintLV->SetVisAttributes(scintVis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zScint), scintLV, "ScintPV",
                          worldLV, false, 0, true);

        G4cout << "[Det] Collimator concept geometry: 2 planes 1 m^2, " << CollimatorDet::kPlaneSep / mm
               << " mm apart; collimator " << CollimatorDet::kNHoles << "x" << CollimatorDet::kNHoles
               << " holes, " << 2.*CollimatorDet::kHoleHalfXY / mm << " mm hole / "
               << CollimatorDet::kWallThick / mm << " mm wall (open area "
               << 100. * std::pow(2.*CollimatorDet::kHoleHalfXY / CollimatorDet::kPitch, 2) << "%), "
               << 2.*CollimatorDet::kCollHalfZ / mm << " mm deep; filter+scint at z = "
               << zFilter / mm << ", " << zScint / mm << " mm" << G4endl;

        return worldPV;
    }

    auto* copper = nist->FindOrBuildMaterial("G4_Cu");
    // Kapton polyimide (C22H10N2O5) built explicitly at 1.413 g/cm3.
    auto* kapton = new G4Material("Kapton", 1.413 * g/cm3, 4);
    kapton->AddElement(nist->FindOrBuildElement("O"),  5);
    kapton->AddElement(nist->FindOrBuildElement("C"), 22);
    kapton->AddElement(nist->FindOrBuildElement("N"),  2);
    kapton->AddElement(nist->FindOrBuildElement("H"), 10);
    // K13D2U pitch-based carbon-fibre cold plate, modelled as carbon at 1.8 g/cm3
    // (the composite radiation length is what matters for MCS).
    // 400 um at 1.8 g/cm3 = 0.072 g/cm2 = 1.686e-3 X0.
    auto* k13d2u = new G4Material("K13D2U", Det::kColdDensity_gcm3 * g/cm3, 1);
    k13d2u->AddElement(nist->FindOrBuildElement("C"), 1.0);

    // Expose standoff and compute world x/y from max bending.
    sStandoff_mm    = ImagingDet::kStandoff / mm;
    sActiveHalfX_mm = Det::kActiveHalfX / mm;   // exposed so the VR beam tracks the real size
    sActiveHalfY_mm = Det::kActiveHalfY / mm;
    sUraniumHalfX_mm = UraniumDet::kHalfX / mm; // exposed so PGA can ray-test the uranium box
    sUraniumHalfY_mm = UraniumDet::kHalfY / mm; // (unconditional on mode -- fixed dimensions,
    sUraniumHalfZ_mm = UraniumDet::kHalfZ / mm; //  only the volume placement itself is mode-gated)

    // Max bending at E_min = 200 MeV over the full drift (standoff + 30 mm source offset).
    const G4double d_mm    = ImagingDet::kStandoff / mm + 30.0;
    const G4double p_min   = std::sqrt(200. * (200. + 2. * 938.272));  // MeV/c
    const G4double r_g_min = p_min / (299.792 * ImagingDet::kBField / tesla);  // m
    const G4double maxBend_mm = (d_mm * d_mm) / (2000.0 * r_g_min);

    // Compute detector alignment from flux-weighted mean bending before sizing the
    // world, since the offset dx grows as d^2 and the world has to know about it.
    const auto [dx_mm, theta_rad] =
        ComputeAlignmentBending(ImagingDet::kBField / tesla,
                                ImagingDet::kStandoff / m);
    sDetOffsetX_mm = dx_mm;                      // expose to PrimaryGeneratorAction
    const G4double detOffset_x = dx_mm * mm;

    // Auto-scale the source half-width so the back-projected open field stays flat
    // across the analysis cell. The reconstructed source position is blurred by
    // sigma_src = sigma_tot * d; if the +-fSrcHX patch doesn't extend well past the
    // +-CELL analysis cell, the edge rolls off inside the cell and t_5sigma is biased
    // high. So grow the patch with standoff:
    //     fSrcHX = max(kSrcHXFloor, CELL + K*sigma_src),   CELL = max(U_HX, 1.75*sigma_src)
    // matching imaging_analysis.m's CELL_MM rule. K = SRC_MARGIN below keeps the loss
    // well under 1%; cost is that events scale as fSrcHX^2.
    {
        const G4double d_km      = ImagingDet::kStandoff / m / 1000.0;
        // Highland at x/X0 = 3.73e-3 per layer (K13D2U at 1.8 g/cm3), flux-weighted
        // over the detected AP9 spectrum.
        const G4double sigMCS    = 1.305, sigPos = 0.612, sigChrom = 0.297 * d_km;  // mrad
        const G4double sigTot    = std::sqrt(sigMCS*sigMCS + sigPos*sigPos + sigChrom*sigChrom);
        const G4double sigSrc_mm = sigTot * 1e-3 * (ImagingDet::kStandoff / mm);    // mm
        const G4double cell_mm   = std::max(200.0, 1.75 * sigSrc_mm);               // U_HX = 200 mm
        constexpr G4double SRC_MARGIN = 2.5;   // sigma_src past the cell
        sSrcHX_mm = std::max(ImagingDet::kSrcHXFloor / mm, cell_mm + SRC_MARGIN * sigSrc_mm);
    }
    const G4double srcHX = sSrcHX_mm * mm;

    // A tight world is a loss-free speed-up: the forward-hemisphere beam fires wide-angle
    // protons that would otherwise traverse a huge world before exiting and can't reach
    // the detector anyway. It must still enclose all three of:
    //   (a) the detector at its bending offset,  |dx| + kActiveHalfX  (kActiveHalfY in y),
    //   (b) the source plane,                    sSrcHX_mm (auto-scaled above),
    //   (c) the bend excursion of a proton launched at the source edge, srcHX + maxBend,
    // otherwise primaries can spawn outside the world and vanish.
    // (d) the collimator's own outer half-extent, when present -- CollV3's lattice (~504
    //     mm) never bound this in practice (smaller than kActiveHalfY), but CollTube's
    //     single tube (~810 mm) does.
    const G4double collOuterHalfXY = !Det::kHasCollimator ? 0.0
        : (Det::kCollStyle == Det::CollStyle::kTube ? CollTube::kOuterHalfXY
                                                      : std::max(CollV3::kCollHalfX, CollV3::kCollHalfY));
    const G4double worldHXY   = std::max({ std::abs(detOffset_x) + Det::kActiveHalfX,
                                           Det::kActiveHalfY,
                                           srcHX + maxBend_mm * mm,
                                           std::abs(detOffset_x) + collOuterHalfXY })
                              + 150.0 * mm;
    // CollTube is ~45 m long regardless of standoff (its length depends only on the
    // design angle + bore size, not kStandoff), so the default kWorldHZ (standoff/2 +
    // 2.5 m) can be too short at small standoffs and let the tube's upstream face poke
    // out of the world. Grow worldHZ to guarantee it always fits.
    G4double worldHZ = ImagingDet::kWorldHZ;
    if (Det::kHasCollimator && Det::kCollStyle == Det::CollStyle::kTube) {
        worldHZ = std::max(worldHZ, 2.0 * CollTube::kCollHalfZ + 1.0 * m);
    }

    G4cout << "[Det] Standoff = " << ImagingDet::kStandoff / m << " m"
           << "  max bend (200 MeV) = " << maxBend_mm << " mm"
           << "  source HX = " << sSrcHX_mm << " mm"
           << "  world HX/HY = " << worldHXY / m << " m"
           << "  world HZ = " << worldHZ / m << " m" << G4endl;
    G4cout << "[Det] Tracker = " << Det::kNStaveX << " x " << Det::kNStaveY
           << " staves, envelope +/-" << Det::kActiveHalfX / mm << " x +/-"
           << Det::kActiveHalfY / mm << " mm, silicon "
           << Det::kNStaveX * Det::kNStaveY * (2*Det::kSensorHalfX) * (2*Det::kSensorHalfY) / cm2
           << " cm2, y row gap " << Det::kRowGapY / mm << " mm" << G4endl;
    // Geometry tag for archive_run.ps1's cross-check, so a run taken with a different
    // tracker build doesn't silently land in the same run_<N>m/ folder as another.
    G4cout << "[Det] GeometryTag = " << Det::kNStaveX << "x" << Det::kNStaveY
           << "_coldplate" << Det::kColdDensity_gcm3 << "gcm3"
           << (DetectorConstruction::kUseStripHodoscope ? "_stripHodo"
               : (G4String("_filter") + (Det::kFilterBeforeLayer2 ? "BeforeL2" : "AfterL2")))
           << (!Det::kHasCollimator ? ""
               : (Det::kCollStyle == Det::CollStyle::kTube ? "_v3collTube"
                  : (CollV3::kAlCollimator ? "_v3collAl" : "_v3coll")))
           << (Det::kTwoLayerViz ? "_TWOLAYERVIZ_NOT_FOR_PRODUCTION" : "")
           << G4endl;

    // rotateY(-theta) makes local +z -> (-sin theta, 0, cos theta) = proton direction
    auto* detRot = new G4RotationMatrix();
    detRot->rotateY(-theta_rad);

    // Active half-widths (cold plate + poly + scint span the full active area).
    const G4double hx = Det::kActiveHalfX;
    const G4double hy = Det::kActiveHalfY;

    // World volume
    auto* worldLV = new G4LogicalVolume(
        new G4Box("World",
                  worldHXY,
                  worldHXY,
                  worldHZ),
        vac, "WorldLV");
    auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV",
                                      nullptr, false, 0, true);
    worldLV->SetVisAttributes(new G4VisAttributes(false));

    // Uranium slab (kImaging only; absent in kImagingOpen for the open-field reference)
    if (sMode == Mode::kImaging) {
        auto* uMat = nist->FindOrBuildMaterial("G4_U");
        auto* uLV  = new G4LogicalVolume(
            new G4Box("UraniumCase",
                      UraniumDet::kHalfX,
                      UraniumDet::kHalfY,
                      UraniumDet::kHalfZ),
            uMat, "UraniumCaseLV");
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ImagingDet::kUraniumZ),
                          uLV, "UraniumCasePV", worldLV, false, 0, true);
        auto* uVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.1, 0.9));
        uVis->SetForceSolid(true);
        uLV->SetVisAttributes(uVis);
    } else if (sMode == Mode::kImagingOpen) {
        // 2026-09-24: vacuum placeholder at the uranium position in the open run. Physically
        // inert (G4_Galactic, same as the world), but it gives the world's smart voxels a
        // daughter at the far (-z) end. Without it every daughter sits at the detector end
        // of the ~standoff-long world, and a proton crossing the empty stretch paid for
        // navigating hundreds of tracker volumes on every step: SAMPLING=direct open runs
        // were ~20x slower per event than the shadow runs (400 m: 17.5 vs 0.8 ms/event).
        auto* phLV = new G4LogicalVolume(
            new G4Box("UraniumPlaceholder", UraniumDet::kHalfX, UraniumDet::kHalfY, UraniumDet::kHalfZ),
            nist->FindOrBuildMaterial("G4_Galactic"), "UraniumPlaceholderLV");
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, ImagingDet::kUraniumZ),
                          phLV, "UraniumPlaceholderPV", worldLV, false, 0, true);
        phLV->SetVisAttributes(new G4VisAttributes(false));
    }

    // Tracker: 3 MAPS layers (ALTAI staves), 20 mm apart.
    // Per layer, beam-facing to downstream: kapton 75 um | copper 18 um | Si 50 um |
    // K13D2U cold plate 400 um. Each stave's 30 mm sensor stack sits on the leg-free
    // centre of a 34 mm cold plate, the outer 2 mm on each side being a 6.6 mm leg.
    // Staves tile at the 34 mm pitch so the cold plates touch (continuous sheet) while
    // adjacent legs meet, leaving a 4 mm bare-silicon gap every stave.
    auto* siVis   = new G4VisAttributes(G4Colour(0.2, 0.5, 1.0, 0.85)); siVis->SetForceSolid(true);
    auto* kapVis  = new G4VisAttributes(G4Colour(0.85, 0.65, 0.2, 0.5)); kapVis->SetForceSolid(true);
    auto* cuVis   = new G4VisAttributes(G4Colour(0.9, 0.45, 0.2, 0.6)); cuVis->SetForceSolid(true);
    auto* coldVis = new G4VisAttributes(G4Colour(0.35, 0.35, 0.35, 0.6)); coldVis->SetForceSolid(true);
    auto* legVis  = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.85)); legVis->SetForceSolid(true);

    const G4double xLeg = Det::kColdHalfX - Det::kLegHalfX;  // 16 mm: leg centre from stave centre

    // Tilt clearance. detRot is applied per volume and G4PVPlacement rotates each solid
    // about its own centre. The cold plate is one full-area sheet centred at detOffset_x
    // while the silicon is kNStaveX separate columns each centred at its own xc, so the
    // z-gap between a column's downstream face and the cold plate's upstream face is
    // uc*sin(theta) (uc = column centre in local x) -- nominally the faces touch, so the
    // tilt drives that negative on one side of the axis and the navigator resolves the
    // overlap in favour of the cold plate, silently shadowing that silicon column.
    // Clearance has to cover the outermost leg (kActiveHalfX + kColdHalfX from the axis),
    // since the legs meet the same full-area sheet from downstream. Gaps are G4_Galactic
    // so this is physics-inert, and kLayerGapZ still sets the exact Si-to-Si spacing.
    const G4double tiltClear =
        (Det::kActiveHalfX + Det::kColdHalfX) * std::sin(theta_rad) + 0.05 * mm;

    // Stack ordering along the beam.
    // v3 (kFilterBeforeLayer2 = true):  L0 | L1 | 25 cm HDPE | L2 | scintillator
    // v2 (false):                       L0 | L1 | L2 | 25 cm HDPE | scintillator
    // Layers 0 and 1 stay kLayerGapZ apart either way, so geo_theta (the L0->L1 chord)
    // and the tracking resolution are unchanged. What v3 changes is that layer 2 now
    // sits downstream of the filter, so its hit requires the proton to have punched
    // through 25 cm of PE and is displaced by filter scattering -- coincidence-only
    // information that never enters the back-projection.
    const G4double kSensorFront = 2.*Det::kKaptonHalfZ + 2.*Det::kCuHalfZ;   // upstream of Si
    const G4double kSensorBack  = 2.*tiltClear + 2.*Det::kColdHalfZ
                                + 2.*Det::kLegHalfZ;                        // downstream of Si
    G4double zSiPlane[3] = {0., 0., 0.}, poly_z = 0., scint_z = 0., stripX_z = 0., stripY_z = 0.;
    zSiPlane[0] = ImagingDet::kSSD0_Z;
    zSiPlane[1] = zSiPlane[0] + (Det::kTwoLayerViz ? Det::kTwoLayerGapZ : Det::kLayerGapZ);
    if (Det::kTwoLayerViz) {
        // Exactly 2 Si planes -- no layer 2 at all (zSiPlane[2] stays 0, unused: the
        // build loop below only runs L=0,1). Hodoscope sits directly behind layer 1.
        const G4double l1Bot = zSiPlane[1] + Det::kSiHalfZ + kSensorBack;
        stripX_z = l1Bot + tiltClear + StripHodo::kGap + StripHodo::kHalfZ;
        stripY_z = stripX_z + StripHodo::kHalfZ + StripHodo::kGap + StripHodo::kHalfZ;
    } else if (DetectorConstruction::kUseStripHodoscope) {
        // No poly filter in this variant (see StripHodo namespace comment), so
        // kFilterBeforeLayer2 doesn't apply -- always the simple v2-style layer spacing.
        zSiPlane[2] = zSiPlane[1] + Det::kLayerGapZ;
        const G4double l2Bot = zSiPlane[2] + Det::kSiHalfZ + kSensorBack;
        stripX_z = l2Bot + tiltClear + StripHodo::kGap + StripHodo::kHalfZ;
        stripY_z = stripX_z + StripHodo::kHalfZ + StripHodo::kGap + StripHodo::kHalfZ;
    } else if (Det::kFilterBeforeLayer2) {
        const G4double l1Bot = zSiPlane[1] + Det::kSiHalfZ + kSensorBack;
        poly_z      = l1Bot + tiltClear + Det::kPoly_HalfZ;
        zSiPlane[2] = poly_z + Det::kPoly_HalfZ + tiltClear + kSensorFront + Det::kSiHalfZ;
        scint_z     = zSiPlane[2] + Det::kSiHalfZ + kSensorBack + tiltClear + Det::kScint_HalfZ;
    } else {
        zSiPlane[2] = zSiPlane[1] + Det::kLayerGapZ;
        const G4double l2Bot = zSiPlane[2] + Det::kSiHalfZ + kSensorBack;
        poly_z      = l2Bot + tiltClear + Det::kPoly_HalfZ;
        scint_z     = poly_z + Det::kPoly_HalfZ + tiltClear + Det::kScint_HalfZ;
    }

    // v3 collimator: immediately upstream of layer 0's kapton front face, sharing the
    // stack's tilt/offset. Own tilt clearance (same mechanism as tiltClear above, but
    // using the collimator's own -- larger -- half-extent, since it's that volume's own
    // face-tilt-across-its-width that matters here). Two styles, see Det::kCollStyle:
    // kLattice = CollV3's periodic square-hole lattice (below); kTube = CollTube's single
    // large square tube around the tracker's whole perimeter (further below).
    if (Det::kHasCollimator && Det::kCollStyle == Det::CollStyle::kLattice) {
        const G4double collClear = CollV3::kCollHalfX * std::sin(theta_rad) + 0.05 * mm;
        const G4double zCollDownstream = zSiPlane[0] - Det::kSiHalfZ - kSensorFront
                                        - CollV3::kCollGap - collClear;
        const G4double zColl = zCollDownstream - CollV3::kCollHalfZ;

        auto* collMat = nist->FindOrBuildMaterial(CollV3::kAlCollimator ? "G4_Al" : "G4_POLYETHYLENE");
        auto* hdpeVis = new G4VisAttributes(CollV3::kAlCollimator ? G4Colour(0.75, 0.78, 0.80, 1.0)
                                                                    : G4Colour(0.9, 0.85, 0.7, 1.0));
        hdpeVis->SetForceSolid(true);

        auto* holeShape = new G4Box("V3HoleShape", CollV3::kHoleHalfXY, CollV3::kHoleHalfXY,
                                    CollV3::kCollHalfZ + 1.0 * mm);   // overshoot for a clean cut
        auto* holesUnion = new G4MultiUnion("V3Holes");
        for (G4int i = 0; i < CollV3::kNHolesX; ++i) {
            const G4double xc = (i - (CollV3::kNHolesX - 1) / 2.) * CollV3::kPitch;
            for (G4int j = 0; j < CollV3::kNHolesY; ++j) {
                const G4double yc = (j - (CollV3::kNHolesY - 1) / 2.) * CollV3::kPitch;
                G4Transform3D tr(G4RotationMatrix(), G4ThreeVector(xc, yc, 0));
                holesUnion->AddNode(*holeShape, tr);
            }
        }
        holesUnion->Voxelize();

        auto* collBlock = new G4Box("V3CollimatorBlock", CollV3::kCollHalfX, CollV3::kCollHalfY,
                                    CollV3::kCollHalfZ);
        auto* collSolid = new G4SubtractionSolid("V3Collimator", collBlock, holesUnion);
        auto* collLV = new G4LogicalVolume(collSolid, collMat, "V3CollimatorLV");
        collLV->SetVisAttributes(hdpeVis);
        new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, 0, zColl), collLV,
                          "V3CollimatorPV", worldLV, false, 0, true);

        G4cout << "[Det] v3 collimator: " << (CollV3::kAlCollimator ? "Al" : "HDPE") << ", "
               << CollV3::kNHolesX << "x" << CollV3::kNHolesY
               << " holes, " << CollV3::kHoleWidth / mm << " mm hole / "
               << CollV3::kWallThick / mm << " mm wall (open area "
               << 100. * CollV3::kOpenArea << "%), " << 2.*CollV3::kCollHalfZ / mm
               << " mm deep, design angle " << CollV3::kDesignDeg << " deg, centre z = "
               << zColl / mm << " mm" << G4endl;
    }
    else if (Det::kHasCollimator && Det::kCollStyle == Det::CollStyle::kTube) {
        const G4double collClear = CollTube::kOuterHalfXY * std::sin(theta_rad) + 0.05 * mm;
        const G4double zCollDownstream = zSiPlane[0] - Det::kSiHalfZ - kSensorFront
                                        - CollTube::kCollGap - collClear;
        // Front face is fixed by the acceptance-angle depth (kCollHalfZ), same as before.
        // Back face is now extended past the WHOLE detector stack (perimeter side
        // shielding, 2026-09-14) instead of stopping at layer 0 -- background can no
        // longer sneak in past the sides of the tracker/trigger stage once it's inside
        // the bore. detStackBack is whichever trigger stage is active (strip hodoscope
        // or poly+scint); both were already computed above alongside zSiPlane[].
        const G4double zFront = zCollDownstream - 2.0 * CollTube::kCollHalfZ;
        const G4double detStackBack = DetectorConstruction::kUseStripHodoscope
            ? (stripY_z + StripHodo::kHalfZ)
            : (scint_z + Det::kScint_HalfZ);
        const G4double zBack  = detStackBack + CollTube::kSideMargin;
        const G4double tubeHalfZ = (zBack - zFront) / 2.0;
        const G4double zColl     = (zFront + zBack) / 2.0;

        auto* collMat = nist->FindOrBuildMaterial(CollTube::kAlWall ? "G4_Al" : "G4_POLYETHYLENE");

        // Mass-optimized taper, 2026-09-15: validated via the WALL_TEST beam
        // (PrimaryGeneratorAction.cc fires the TRUE un-redirected AP9 PAD directly at the
        // wall -- see the dated investigation). Punch-through rate at a FIXED thickness
        // falls off ~exponentially with distance from the detector (measured ~700 mm
        // decay length: 12.3% within 300 mm of the detector, down to 0.00% beyond ~4.3 m,
        // at a uniform 25%-thickness diagnostic build) -- a real self-collimation effect,
        // not angle-grazing: a punch-through far from the detector still has to cross the
        // remaining bore length without clipping another zone to ever reach the tracker.
        // distFromBack_mm are CUMULATIVE boundaries from the detector-facing end (zBack);
        // thickFrac is that zone's wall thickness as a fraction of the full kWallThick.
        //
        // RE-SEGMENTED 2026-09-17 (first pass) to uniform 1000 mm zones, to match the
        // widened kDesignDeg = 11.7 deg design -- SUPERSEDED same day, see the next note:
        // holding the near-detector (100%) zone at full thickness for a full 1000 mm
        // instead of the original design's 500 mm cost ~836 kg on its own, which very
        // nearly cancelled the mass saved by the shorter bore (net +121 kg vs the 8 deg
        // design, despite being 1.82 m shorter) -- not all length costs the same mass,
        // since the near-detector zone's cross-section is far larger than the far zones'.
        //
        // RE-DESIGNED 2026-09-17 (second pass, same day): segmented far more finely, in
        // 200 mm steps, with the thickness fraction interpolated LINEARLY between the
        // same validated anchor points from the original 8 deg/500-1500-2500-3500 mm
        // taper (0mm->100%, 500mm->75%, 1500mm->50%, 2500mm->32%, 3500mm->20%, held flat
        // beyond 3500 mm -- no data past that point) instead of holding each old bin's
        // value constant across a whole new 1000 mm zone. Each 200 mm zone uses the
        // fraction at its NEAR edge (closest to the detector, highest flux) for the same
        // conservative-rebinning reason as the first pass. Result: 3.07 t (20 zones) vs
        // 4.18 t (4 zones/1000mm) -- 27% lighter, and 24% lighter than the original
        // 8 deg/5-zone design (4.06 t) despite the wider FOV. Zone count and lengths are
        // DERIVED from totalLen (itself derived from kDesignDeg), not hardcoded, so this
        // stays correct automatically if the design angle changes again.
        // **NOT YET RE-VALIDATED**: still a first-pass interpolation of the previously
        // measured anchor points, not a fresh WALL_TEST measurement at these finer
        // boundaries -- re-run WALL_TEST=1 after building to confirm zero sub-200 MeV
        // leakage in every zone before trusting it for a production run.
        struct ThickAnchor { G4double dist_mm; G4double frac; };
        const ThickAnchor kAnchors[] = {
            {    0.0 * mm, 1.00 },
            {  500.0 * mm, 0.75 },
            { 1500.0 * mm, 0.50 },
            { 2500.0 * mm, 0.32 },
            { 3500.0 * mm, 0.20 },   // held flat beyond this point -- no data past here
        };
        const G4int nAnchors = sizeof(kAnchors) / sizeof(kAnchors[0]);
        auto thickFracAt = [&](G4double distFromDetector_mm) -> G4double {
            if (distFromDetector_mm <= kAnchors[0].dist_mm) return kAnchors[0].frac;
            for (G4int i = 0; i < nAnchors - 1; ++i) {
                if (distFromDetector_mm >= kAnchors[i].dist_mm &&
                    distFromDetector_mm <= kAnchors[i + 1].dist_mm) {
                    const G4double t = (distFromDetector_mm - kAnchors[i].dist_mm)
                                      / (kAnchors[i + 1].dist_mm - kAnchors[i].dist_mm);
                    return kAnchors[i].frac + t * (kAnchors[i + 1].frac - kAnchors[i].frac);
                }
            }
            return kAnchors[nAnchors - 1].frac;   // held flat beyond the last anchor
        };

        const G4double totalLen = zBack - zFront;
        constexpr G4double kZoneStep = 200.0 * mm;
        const G4int nZones = G4int(std::ceil(totalLen / kZoneStep));
        constexpr G4double kZoneEps = 0.01 * mm;   // tiny shrink per zone so joins don't overlap

        fCollTubeWallZones_LV.clear();
        G4double prevD = 0.0;
        G4double reportMaxOuter = 0.0;
        for (G4int iz = 0; iz < nZones; ++iz) {
            const G4double thisD = std::min(prevD + kZoneStep, totalLen);
            if (thisD <= prevD) break;
            const G4double zoneLen   = thisD - prevD;
            const G4double zoneHalfZ = zoneLen / 2.0 - kZoneEps;
            const G4double zCentre   = zBack - (prevD + thisD) / 2.0;
            const G4double thickFrac = thickFracAt(prevD);   // near-edge (conservative) value
            const G4double outerHXY  = CollTube::kBoreHalfXY + thickFrac * CollTube::kWallThick;
            reportMaxOuter = std::max(reportMaxOuter, outerHXY);

            const G4String tag = "TubeZone" + std::to_string(iz);
            auto* outer = new G4Box(tag + "Outer", outerHXY, outerHXY, zoneHalfZ);
            auto* bore  = new G4Box(tag + "Bore", CollTube::kBoreHalfXY, CollTube::kBoreHalfXY,
                                    zoneHalfZ + 1.0 * mm);
            auto* solid = new G4SubtractionSolid(tag, outer, bore);
            auto* lv = new G4LogicalVolume(solid, collMat, tag + "LV");
            auto* vis = new G4VisAttributes(CollTube::kAlWall
                ? G4Colour(0.55 + 0.30*thickFrac, 0.60, 0.65, 1.0)
                : G4Colour(0.65 + 0.30*thickFrac, 0.65, 0.45, 1.0));
            vis->SetForceSolid(true);
            lv->SetVisAttributes(vis);
            // Every zone shares the SAME PV name (TrackerSD's planeID lookup keys on name,
            // not copy number) so the wall-hit investigation sees the whole taper as one.
            new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, 0, zCentre), lv,
                              "V3CollimatorTubePV", worldLV, false, iz, true);
            fCollTubeWallZones_LV.push_back(lv);

            G4cout << "[Det]   tube zone " << iz << ": " << thickFrac*100 << "% ("
                   << thickFrac*CollTube::kWallThick/mm << " mm), "
                   << zoneLen/mm << " mm long, " << prevD/mm << "-" << thisD/mm
                   << " mm from detector" << G4endl;
            prevD = thisD;
        }

        // Expose the collimator's full z-span + outer half-width (largest zone, i.e. the
        // full-thickness one) for the wall-test beam (WALL_TEST env var) to position its
        // source plane safely outside every zone.
        sCollZFront_mm = zFront / mm;
        sCollZBack_mm  = zBack / mm;
        sCollOuterHalfXY_mm = reportMaxOuter / mm;

        G4cout << "[Det] v3 tube collimator: " << (CollTube::kAlWall ? "Al" : "HDPE")
               << ", bore +/-" << CollTube::kBoreHalfXY / mm << " mm, " << nZones
               << "-zone taper, total " << 2.*tubeHalfZ / m << " m, on-axis design angle "
               << CollTube::kDesignDeg << " deg, centre z = " << zColl / mm << " mm" << G4endl;
    }

    G4LogicalVolume* siLV[3] = { nullptr, nullptr, nullptr };
    const G4int nSiLayers = Det::kTwoLayerViz ? 2 : 3;
    for (G4int L = 0; L < nSiLayers; ++L) {
        const G4double zSi  = zSiPlane[L];                                 // Si plane
        const G4double zCu  = zSi - Det::kSiHalfZ - Det::kCuHalfZ;         // copper (upstream)
        const G4double zKap = zCu - Det::kCuHalfZ - Det::kKaptonHalfZ;     // kapton (upstream)
        // + tiltClear on each face-to-face join with the full-area cold-plate sheet.
        const G4double zCol = zSi + Det::kSiHalfZ + tiltClear + Det::kColdHalfZ;
        const G4double zLeg = zCol + Det::kColdHalfZ + tiltClear + Det::kLegHalfZ;
        const G4String s = std::to_string(L);

        // Continuous 400 um K13D2U cold-plate sheet (touching 34 mm cells -> full area).
        auto* colLV = new G4LogicalVolume(new G4Box("Cold"+s, hx, hy, Det::kColdHalfZ), k13d2u, "ColdLV"+s);
        colLV->SetVisAttributes(coldVis);
        new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, 0, zCol), colLV, "Cold"+s+"PV", worldLV, false, 0, true);

        // Per-stave sensor stack (30 x 150 mm) plus two 2 mm legs at the cold-plate edges.
        // The sensor sub-slabs tile in y at the 160 mm pitch, leaving kRowGapY (10 mm) of
        // bare cold plate between rows. The legs run the full cold-plate height and so
        // bridge the inter-row gaps.
        auto* kapLV = new G4LogicalVolume(new G4Box("Kapton"+s, Det::kSensorHalfX, Det::kSensorHalfY, Det::kKaptonHalfZ), kapton, "KaptonLV"+s);
        kapLV->SetVisAttributes(kapVis);
        auto* cuLV  = new G4LogicalVolume(new G4Box("Cu"+s, Det::kSensorHalfX, Det::kSensorHalfY, Det::kCuHalfZ), copper, "CuLV"+s);
        cuLV->SetVisAttributes(cuVis);
        siLV[L] = new G4LogicalVolume(new G4Box("SSD"+s, Det::kSensorHalfX, Det::kSensorHalfY, Det::kSiHalfZ), Si, "SSD"+s+"LV");
        siLV[L]->SetVisAttributes(siVis);
        auto* legLV = new G4LogicalVolume(new G4Box("Leg"+s, Det::kLegHalfX, hy, Det::kLegHalfZ), k13d2u, "LegLV"+s);
        legLV->SetVisAttributes(legVis);

        for (G4int i = 0; i < Det::kNStaveX; ++i) {
            const G4double xc = (i - (Det::kNStaveX - 1) / 2.) * Det::kStavePitchX;  // stave centre
            new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc + xLeg, 0, zLeg), legLV, "Leg"+s+"PV", worldLV, false, 2*i,   false);
            new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc - xLeg, 0, zLeg), legLV, "Leg"+s+"PV", worldLV, false, 2*i+1, false);
            for (G4int j = 0; j < Det::kNStaveY; ++j) {
                const G4double yc = (j - (Det::kNStaveY - 1) / 2.) * Det::kStavePitchY;  // row centre
                const G4int    cp = i * Det::kNStaveY + j;                               // unique copy no.
                new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc, yc, zKap), kapLV,   "Kapton"+s+"PV", worldLV, false, cp, false);
                new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc, yc, zCu ), cuLV,    "Cu"+s+"PV",     worldLV, false, cp, false);
                new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc, yc, zSi ), siLV[L], "SSD"+s+"PV",    worldLV, false, cp, false);
            }
        }
    }
    fSSD0_LV = siLV[0];  fSSD1_LV = siLV[1];  fSSD2_LV = siLV[2];

    if (DetectorConstruction::kUseStripHodoscope) {
        // X-Y strip hodoscope: two crossed layers of 1 cm x 1 cm x 1 cm scintillator
        // strips, replacing the poly filter + single scintillator slab. stripX_z/
        // stripY_z were already computed above. One G4LogicalVolume per layer, reused
        // across every strip placement (same pattern as the SSD stave grid) -- each
        // placement gets a unique copy number 0..N-1, which TrackerSD reads back as the
        // strip index via GetCopyNumber() (see its planeID 5/6 handling).
        auto* stripXVis = new G4VisAttributes(G4Colour(1.0, 0.6, 0.2, 0.9));
        stripXVis->SetForceSolid(true);
        auto* stripYVis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.9));
        stripYVis->SetForceSolid(true);

        fStripX_LV = new G4LogicalVolume(
            new G4Box("StripX", StripHodo::kHalfZ, Det::kActiveHalfY, StripHodo::kHalfZ),
            scintMat, "StripXLV");
        fStripX_LV->SetVisAttributes(stripXVis);
        for (G4int i = 0; i < StripHodo::kNStripsX; ++i) {
            const G4double xc = (i - (StripHodo::kNStripsX - 1) / 2.) * StripHodo::kPitch;
            new G4PVPlacement(detRot, G4ThreeVector(detOffset_x + xc, 0, stripX_z),
                              fStripX_LV, "StripXPV", worldLV, false, i, false);
        }

        fStripY_LV = new G4LogicalVolume(
            new G4Box("StripY", Det::kActiveHalfX, StripHodo::kHalfZ, StripHodo::kHalfZ),
            scintMat, "StripYLV");
        fStripY_LV->SetVisAttributes(stripYVis);
        for (G4int j = 0; j < StripHodo::kNStripsY; ++j) {
            const G4double yc = (j - (StripHodo::kNStripsY - 1) / 2.) * StripHodo::kPitch;
            new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, yc, stripY_z),
                              fStripY_LV, "StripYPV", worldLV, false, j, false);
        }

        G4cout << "[Det] strip hodoscope: " << StripHodo::kNStripsX << " x-strips + "
               << StripHodo::kNStripsY << " y-strips, " << StripHodo::kPitch / cm
               << " cm pitch/thickness, centre z = " << stripX_z / mm << " / "
               << stripY_z / mm << " mm" << G4endl;
    } else {
        // Polyethylene range filter (25 cm, threshold ~195 MeV) + scintillator.
        // poly_z / scint_z were already computed above since in v3 the filter sits inside
        // the tracker stack. Both get tiltClear against their neighbours for the same
        // reason as the intra-layer joins: a full-area slab rotated about its own centre
        // bites into the adjacent legs at a flush face. A true zero gap would need the
        // whole stack placed as one rotated assembly instead of each volume rotated
        // independently.
        fPoly_LV = new G4LogicalVolume(
            new G4Box("Poly", hx, hy, Det::kPoly_HalfZ), poly, "PolyLV");
        new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, 0, poly_z),
                          fPoly_LV, "PolyPV", worldLV, false, 0, true);
        auto* polyVis = new G4VisAttributes(G4Colour(0.5, 0.8, 1.0, 0.4));
        polyVis->SetForceSolid(true);
        fPoly_LV->SetVisAttributes(polyVis);

        fScint_LV = new G4LogicalVolume(
            new G4Box("Scint", hx, hy, Det::kScint_HalfZ), scintMat, "ScintLV");
        new G4PVPlacement(detRot, G4ThreeVector(detOffset_x, 0, scint_z),
                          fScint_LV, "ScintPV", worldLV, false, 0, true);
        auto* scintVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.2, 0.9));
        scintVis->SetForceSolid(true);
        fScint_LV->SetVisAttributes(scintVis);
    }

    return worldPV;
}

void DetectorConstruction::ConstructSDandField()
{
    // Collimator concept geometry is visualization-only, no sensitive detectors or field.
    if (sMode == Mode::kCollimatorConcept) return;

    auto* sd = new TrackerSD("TrackerSD", "TrackerHC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sd);

    SetSensitiveDetector(fSSD0_LV,  sd);
    SetSensitiveDetector(fSSD1_LV,  sd);
    if (fSSD2_LV) SetSensitiveDetector(fSSD2_LV,  sd);  // null in Det::kTwoLayerViz
    if (DetectorConstruction::kUseStripHodoscope) {
        SetSensitiveDetector(fStripX_LV, sd);
        SetSensitiveDetector(fStripY_LV, sd);
    } else {
        SetSensitiveDetector(fPoly_LV,  sd);
        SetSensitiveDetector(fScint_LV, sd);
    }
    for (auto* zoneLV : fCollTubeWallZones_LV) {   // 2026-09-15: multi-zone tube wall (tube style only)
        SetSensitiveDetector(zoneLV, sd);
    }

    // Uniform geomagnetic field, exact helical propagation. G4ExactHelixStepper solves
    // the Lorentz equation analytically for a uniform field, so bending has no step-size
    // integration error. Called once per worker thread (G4TransportationManager is
    // thread-local).
    auto* field    = new G4UniformMagField(G4ThreeVector(0., ImagingDet::kBField, 0.));
    auto* equation = new G4Mag_UsualEqRhs(field);
    auto* stepper  = new G4ExactHelixStepper(equation);
    auto* chordFinder = new G4ChordFinder(field, 1.0 * mm, stepper);
    auto* fldMgr = G4TransportationManager::GetTransportationManager()
                       ->GetFieldManager();
    fldMgr->SetDetectorField(field);
    fldMgr->SetChordFinder(chordFinder);
}
