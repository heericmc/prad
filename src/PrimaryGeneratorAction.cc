// Beam source term for a passive standoff detection sensor (orbital treaty
// verification), NOT a weapon/implosion simulation — see README.md. Primaries here are
// ambient Van Allen belt (AP9) trapped protons, the natural background flux the target
// object casts a shadow in; this is not an accelerator or explosive-driven source.

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <cstdlib>
#include <limits>
#include <utility>

// Vis-only flag: set by MAIN.cc before Initialize() (see header).
bool PrimaryGeneratorAction::sVisMode = false;

// Line-of-sight pitch angle: angle between the target-detector line and the local B
// field (here +y). The narrow FOV means every detected proton has this pitch angle, so
// the sampled energy spectrum is the directional flux J(alpha_LOS, E). 90 deg = LOS
// perpendicular to B (PAD peak). Pitch angles in the loss cone (< 60 or > 120 deg)
// give zero flux.
namespace { constexpr G4double kPitchLOS_deg = 90.0; }

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fGun = new G4ParticleGun(1);
    auto* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    fGun->SetParticleDefinition(proton);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    // Prefer the pitch-angle-resolved IRENE directional flux J(alpha, E), sampling
    // energies at the line-of-sight pitch angle. Fall back to the flat i316 spectrum
    // if the IRENE file isn't present next to the executable.
    fPitchLOS_deg = kPitchLOS_deg;
    const bool useVR = [] {
        const char* e = std::getenv("PROD_BEAM");
        return e && std::string(e) == std::string("vr");
    }();
    const G4String ireneFile = "run1.AP9.output_mean_flux.txt";
    std::ifstream ireneTest(ireneFile.c_str());
    if (ireneTest.is_open()) {
        ireneTest.close();
        LoadDirectionalSpectrum(ireneFile, fPitchLOS_deg);
        LoadPAD(ireneFile);   // pitch-angle distribution for the forward-hemisphere beam
    } else {
        G4cout << "[PGA] IRENE file '" << ireneFile
               << "' not found, falling back to diff_flux_AP9_i316.csv" << G4endl;
        LoadSpectrum("diff_flux_AP9_i316.csv");
        if (!useVR)
            G4Exception("PrimaryGeneratorAction", "ForwardNeedsIRENE", FatalException,
                        "Forward-hemisphere AP9 beam needs the IRENE file for the pitch-angle "
                        "distribution (set PROD_BEAM=vr for the legacy aiming beam without it).");
    }

    fMeanDxBend_mm = DetectorConstruction::GetDetOffsetX_mm();

    // Derive beam geometry from standoff.
    const G4double standoff = DetectorConstruction::GetStandoff_mm() * mm;
    fZdet = 0.5 * standoff;
    fZsrc = -(0.5 * standoff + 30.0 * mm);

    // Source (illumination-plane) half-width, auto-scaled with standoff in
    // DetectorConstruction so the back-projected open field stays flat across the
    // analysis cell. Read from there rather than hard-coded so the source, the world
    // sizing, and the analysis can't drift apart.
    fSrcHX = DetectorConstruction::GetSrcHX_mm() * mm;

    // Uranium slab geometry, for the WouldHitUranium() pre-check. z-centre derives from
    // the standoff alone (mirrors ImagingDet::kUraniumZ = -0.5*kStandoff); half-extents
    // come from DetectorConstruction so they can't drift from the built geometry.
    fUraniumZ     = -0.5 * standoff;
    fUraniumHalfX = DetectorConstruction::GetUraniumHalfX_mm() * mm;
    fUraniumHalfY = DetectorConstruction::GetUraniumHalfY_mm() * mm;
    fUraniumHalfZ = DetectorConstruction::GetUraniumHalfZ_mm() * mm;

    const char* collEnv = std::getenv("COLLIMATOR_DEG");
    // 2026-09-16: default fallback 5.0 -> 8.0 to match CollTube::kDesignDeg's widening
    // 5 -> 8 deg (re-derived from the pn=1e-3 confusion-probability bound, see that flag's
    // comment); a run launched without explicitly setting COLLIMATOR_DEG would otherwise
    // silently under-sample this bore.
    G4cout << "[PGA] Standoff = " << standoff / m << " m"
           << "  detector x = " << fMeanDxBend_mm << " mm"
           << "  source ±"      << fSrcHX / mm     << " mm"
           << "  COLLIMATOR_DEG = " << (collEnv ? std::atof(collEnv) : 8.0) << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

void PrimaryGeneratorAction::LoadSpectrum(const G4String& filename)
{
    std::ifstream f(filename.c_str());
    if (!f.is_open())
        G4Exception("PrimaryGeneratorAction::LoadSpectrum", "FileNotFound",
                    FatalException,
                    ("Cannot open spectrum file: " + filename).c_str());

    std::vector<G4double> energies, fluxes;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        G4double e, j;
        if (!(ss >> e >> j)) continue;
        if (e >= 200. && j > 0.) { energies.push_back(e); fluxes.push_back(j); }
    }

    if (energies.size() < 2)
        G4Exception("PrimaryGeneratorAction::LoadSpectrum", "TooFewPoints",
                    FatalException, "Fewer than 2 spectrum points at E >= 200 MeV.");

    std::vector<G4double> cdf(energies.size(), 0.);
    for (std::size_t i = 1; i < energies.size(); ++i) {
        G4double dE  = energies[i] - energies[i-1];
        G4double avg = 0.5 * (fluxes[i] + fluxes[i-1]);
        cdf[i] = cdf[i-1] + dE * avg;
    }
    G4double total = cdf.back();
    for (auto& c : cdf) c /= total;

    fEnergies = std::move(energies);
    fCDF      = std::move(cdf);

    G4cout << "[PGA] Loaded " << fEnergies.size()
           << " spectral points (E_min=" << fEnergies.front()
           << " MeV, E_max=" << fEnergies.back() << " MeV)" << G4endl;
}

void PrimaryGeneratorAction::LoadDirectionalSpectrum(const G4String& filename,
                                                     G4double pitchDeg)
{
    // Parse the IRENE mean-flux table: directional differential flux J(alpha, E) in
    // #/cm^2/sr/s/MeV, laid out as (13 pitch angles) x (N time steps), each row with
    // 32 differential-energy columns. Time-average over all steps, extract the row at
    // the line-of-sight pitch angle, then build the same E >= 200 MeV sampling CDF as
    // LoadSpectrum. Only the CDF shape enters the sim; absolute normalisation is
    // handled downstream in the rate budget.
    std::ifstream f(filename.c_str());
    if (!f.is_open())
        G4Exception("PrimaryGeneratorAction::LoadDirectionalSpectrum", "FileNotFound",
                    FatalException,
                    ("Cannot open IRENE flux file: " + filename).c_str());

    std::vector<G4double> egrid;                          // 32-point energy grid [MeV]
    std::map<G4double, std::vector<G4double>> accum;      // pitch[deg] -> summed flux[E]
    std::string line;

    while (std::getline(f, line)) {
        if (line.empty()) continue;

        // Energy grid comes from the header comment "# Energy levels (MeV): ...".
        if (line[0] == '#') {
            const auto pos = line.find("Energy levels");
            if (pos != std::string::npos && egrid.empty()) {
                const std::string nums = line.substr(line.find(':', pos) + 1);
                std::istringstream ss(nums);
                G4double e;
                while (ss >> e) egrid.push_back(e);
            }
            continue;
        }

        // Data row: datetime, posx, posy, posz, pitch, then one flux per energy.
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        G4double mjd, px, py, pz, pitch;
        if (!(ss >> mjd >> px >> py >> pz >> pitch)) continue;

        std::vector<G4double>& row = accum[pitch];
        if (row.empty()) row.assign(egrid.size(), 0.);
        G4double v;
        for (std::size_t j = 0; j < row.size() && (ss >> v); ++j) row[j] += v;
    }

    if (egrid.size() < 2 || accum.empty())
        G4Exception("PrimaryGeneratorAction::LoadDirectionalSpectrum", "ParseFailed",
                    FatalException, "Could not parse energy grid or flux rows.");

    // Select the flux row at pitchDeg: exact grid point if present, else linear
    // interpolation between the two bracketing pitch angles.
    std::vector<G4double> jflux(egrid.size(), 0.);
    auto hi = accum.lower_bound(pitchDeg - 1e-6);
    if (hi != accum.end() && std::abs(hi->first - pitchDeg) < 1e-6) {
        jflux = hi->second;                               // exact grid match (e.g. 90 deg)
    } else if (hi == accum.begin() || hi == accum.end()) {
        G4Exception("PrimaryGeneratorAction::LoadDirectionalSpectrum", "PitchOutOfRange",
                    FatalException,
                    ("Requested pitch angle outside the IRENE grid: "
                     + std::to_string(pitchDeg) + " deg").c_str());
    } else {
        auto lo = std::prev(hi);
        const G4double w = (pitchDeg - lo->first) / (hi->first - lo->first);
        for (std::size_t j = 0; j < jflux.size(); ++j)
            jflux[j] = (1. - w) * lo->second[j] + w * hi->second[j];
    }

    // Keep E >= 200 MeV (analysis threshold; matches the poly-filter cutoff).
    std::vector<G4double> energies, fluxes;
    for (std::size_t j = 0; j < egrid.size(); ++j)
        if (egrid[j] >= 200. && jflux[j] > 0.) {
            energies.push_back(egrid[j]);
            fluxes.push_back(jflux[j]);
        }

    if (energies.size() < 2)
        G4Exception("PrimaryGeneratorAction::LoadDirectionalSpectrum", "EmptySpectrum",
                    FatalException,
                    ("No flux at E >= 200 MeV for pitch = " + std::to_string(pitchDeg)
                     + " deg -- is the line-of-sight in the loss cone?").c_str());

    std::vector<G4double> cdf(energies.size(), 0.);
    for (std::size_t i = 1; i < energies.size(); ++i) {
        G4double dE  = energies[i] - energies[i-1];
        G4double avg = 0.5 * (fluxes[i] + fluxes[i-1]);
        cdf[i] = cdf[i-1] + dE * avg;
    }
    G4double total = cdf.back();
    for (auto& c : cdf) c /= total;

    fEnergies = std::move(energies);
    fCDF      = std::move(cdf);

    G4cout << "[PGA] Loaded IRENE directional flux J(alpha=" << pitchDeg << " deg, E): "
           << fEnergies.size() << " points at E >= 200 MeV (E_max="
           << fEnergies.back() << " MeV), time-averaged over "
           << (accum.count(pitchDeg) ? "grid rows" : "interpolated rows") << G4endl;
}

G4double PrimaryGeneratorAction::SampleEnergy() const
{
    const G4double u = G4UniformRand();
    auto it = std::lower_bound(fCDF.begin(), fCDF.end(), u);
    if (it == fCDF.end()) --it;
    const std::size_t idx = std::distance(fCDF.begin(), it);
    if (idx == 0) return fEnergies[0];

    const G4double dC   = fCDF[idx] - fCDF[idx-1];
    const G4double frac = (dC > 0.) ? (u - fCDF[idx-1]) / dC : 0.;
    return fEnergies[idx-1] + frac * (fEnergies[idx] - fEnergies[idx-1]);
}

void PrimaryGeneratorAction::LoadPAD(const G4String& filename)
{
    // Build the pitch-angle distribution J(alpha) = integral over E>=200 MeV of the
    // IRENE directional flux, time-averaged, then a fine inverse-CDF for sampling the
    // incident pitch angle. Same file/format as LoadDirectionalSpectrum.
    std::ifstream f(filename.c_str());
    if (!f.is_open())
        G4Exception("PrimaryGeneratorAction::LoadPAD", "FileNotFound", FatalException,
                    ("Cannot open IRENE flux file: " + filename).c_str());

    std::vector<G4double> egrid;
    std::map<G4double, std::vector<G4double>> accum;      // pitch[deg] -> summed flux[E]
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            const auto pos = line.find("Energy levels");
            if (pos != std::string::npos && egrid.empty()) {
                const std::string nums = line.substr(line.find(':', pos) + 1);
                std::istringstream ss(nums);
                G4double e;
                while (ss >> e) egrid.push_back(e);
            }
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        G4double mjd, px, py, pz, pitch;
        if (!(ss >> mjd >> px >> py >> pz >> pitch)) continue;
        std::vector<G4double>& row = accum[pitch];
        if (row.empty()) row.assign(egrid.size(), 0.);
        G4double v;
        for (std::size_t j = 0; j < row.size() && (ss >> v); ++j) row[j] += v;
    }

    // J(alpha) = trapezoidal integral of flux over E >= 200 MeV, per pitch grid point.
    std::map<G4double, G4double> Jofalpha;
    for (auto& kv : accum) {
        G4double J = 0.;
        for (std::size_t j = 1; j < egrid.size(); ++j)
            if (egrid[j] >= 200. && egrid[j-1] >= 200.)
                J += 0.5 * (kv.second[j] + kv.second[j-1]) * (egrid[j] - egrid[j-1]);
        Jofalpha[kv.first] = J;
    }

    auto Jinterp = [&](G4double a) -> G4double {
        auto hi = Jofalpha.lower_bound(a);
        if (hi == Jofalpha.begin() || hi == Jofalpha.end()) return 0.;
        auto lo = std::prev(hi);
        const G4double w = (a - lo->first) / (hi->first - lo->first);
        return (1. - w) * lo->second + w * hi->second;
    };

    fPitchGrid.clear(); fPitchCDF.clear();
    G4double cum = 0.; G4double prevA = 40., prevJ = 0.;
    for (G4double a = 40.; a <= 140.0001; a += 0.5) {
        // Weight J(alpha) by sin(alpha): a pitch ring subtends 2*pi*sin(alpha)*d(alpha)
        // of solid angle, so the number of protons per unit pitch is J(alpha)*sin(alpha).
        // Sampling alpha from this CDF with uniform gyrophase then gives a per-solid-
        // angle direction density proportional to J(alpha), as required.
        const G4double J = Jinterp(a) * std::sin(a * deg);
        cum += 0.5 * (J + prevJ) * (a - prevA);
        fPitchGrid.push_back(a);
        fPitchCDF.push_back(cum);
        prevA = a; prevJ = J;
    }
    if (cum <= 0.)
        G4Exception("PrimaryGeneratorAction::LoadPAD", "EmptyPAD", FatalException,
                    "Pitch-angle distribution integrated to zero.");
    for (auto& c : fPitchCDF) c /= cum;

    G4cout << "[PGA] Built PAD over " << Jofalpha.size() << " IRENE angles; sampling pitch in ["
           << fPitchGrid.front() << ", " << fPitchGrid.back() << "] deg" << G4endl;
}

G4double PrimaryGeneratorAction::SamplePitch() const
{
    const G4double u = G4UniformRand();
    auto it = std::lower_bound(fPitchCDF.begin(), fPitchCDF.end(), u);
    if (it == fPitchCDF.end()) --it;
    const std::size_t idx = std::distance(fPitchCDF.begin(), it);
    if (idx == 0) return fPitchGrid[0];
    const G4double dC   = fPitchCDF[idx] - fPitchCDF[idx-1];
    const G4double frac = (dC > 0.) ? (u - fPitchCDF[idx-1]) / dC : 0.;
    return fPitchGrid[idx-1] + frac * (fPitchGrid[idx] - fPitchGrid[idx-1]);
}

G4ThreeVector PrimaryGeneratorAction::SampleNarrowConeDir(G4double halfAngleRad) const
{
    // Uniform in solid angle within halfAngleRad of +z. For ~1 deg this makes dir.y()
    // (and hence pitch = acos(dir.y())) land within ~1 deg of 90 -- exactly the "narrow
    // FOV pins pitch near LOS" regime the energy CDF (built at kPitchLOS_deg = 90) and
    // the sin(alpha)*J(alpha) PAD weighting are both effectively constant over, so a
    // uniform draw here is a faithful, much cheaper stand-in for rejection-sampling the
    // full PAD down to this tiny a cone. cos-incidence weighting is skipped: at 1 deg
    // dir.z() >= cos(1 deg) = 0.99985, a <0.02% effect, negligible against everything
    // else in this model.
    const G4double cosMax   = std::cos(halfAngleRad);
    const G4double cosTheta = 1.0 - G4UniformRand() * (1.0 - cosMax);
    const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));
    const G4double phi      = 2.0 * std::acos(-1.0) * G4UniformRand();
    return G4ThreeVector(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
}

bool PrimaryGeneratorAction::WouldHitUranium(const G4ThreeVector& pos, const G4ThreeVector& dir) const
{
    // Standard slab method: intersect the ray with each pair of box faces in turn,
    // shrinking the surviving parametric range [tEnter, tExit]; empty range or an exit
    // behind the source means no hit.
    const G4double lo[3] = { -fUraniumHalfX, -fUraniumHalfY, fUraniumZ - fUraniumHalfZ };
    const G4double hi[3] = {  fUraniumHalfX,  fUraniumHalfY, fUraniumZ + fUraniumHalfZ };
    const G4double p[3]  = { pos.x(), pos.y(), pos.z() };
    const G4double d[3]  = { dir.x(), dir.y(), dir.z() };

    G4double tEnter = -std::numeric_limits<G4double>::max();
    G4double tExit  =  std::numeric_limits<G4double>::max();
    for (int i = 0; i < 3; ++i) {
        if (std::abs(d[i]) < 1e-12) {
            if (p[i] < lo[i] || p[i] > hi[i]) return false;
            continue;
        }
        G4double t1 = (lo[i] - p[i]) / d[i];
        G4double t2 = (hi[i] - p[i]) / d[i];
        if (t1 > t2) std::swap(t1, t2);
        tEnter = std::max(tEnter, t1);
        tExit  = std::min(tExit,  t2);
        if (tEnter > tExit) return false;
    }
    return tExit >= 0.;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    // Legacy variance-reduction (inverted-aiming) beam: env PROD_BEAM=vr. Aims every
    // proton at the detector (~100% efficient) but structurally cannot sample uranium
    // scatter-in; kept for fast direct-signal / rate-budget studies only.
    static const bool kVR = [] {
        const char* e = std::getenv("PROD_BEAM");
        return e && std::string(e) == std::string("vr");
    }();
    if (kVR) { GenerateInvertedBeam(event); return; }

    // Wall-thickness investigation: env WALL_TEST=1 bypasses everything else below.
    static const bool kWallTest = [] {
        const char* e = std::getenv("WALL_TEST");
        return e && std::string(e) == std::string("1");
    }();
    if (kWallTest) { GenerateWallTestBeam(event); return; }

    // Physical AP9 flux, collimator-aware source term (2026-09-10 redesign; widened to 2
    // deg and paired with a REAL HDPE lattice collimator, Det::kHasCollimator /
    // CollV3 in DetectorConstruction.cc, as of the "v3" redefinition same day). When the
    // real collimator is built, COLLIMATOR_DEG is now a pure sampling-efficiency cone
    // matched to (slightly wider than, for safety) its physical design angle -- Geant4
    // still tracks genuine absorption/wall-scattering through the lattice, this just
    // avoids wasting compute on directions the real geometry has ~no chance of passing.
    // With no collimator built (Det::kHasCollimator = false), it's the same emulation as
    // before: no absorber exists, this alone stands in for the angular acceptance one
    // would give.
    //   - open run (no uranium): every proton must arrive within the collimator, so
    //     sample directly from the narrow cone -- no point sampling the full PAD only
    //     to reject nearly all of it.
    //   - shadow run: sample the FULL forward hemisphere from the PAD (dir.z() > 0 is
    //     guaranteed by construction; only the cos-incidence draw can reject), so
    //     uranium scatter-in (which needs off-axis starting directions) is still caught.
    //     Then classify: if this trajectory would actually hit the uranium, keep it as
    //     drawn -- Geant4 tracks the real (possibly scattering) physics from there. If
    //     it would NOT hit the uranium, it can only ever reach the detector directly,
    //     which the real collimator would reject unless it was already within COLLIMATOR_DEG --
    //     so redraw its direction from the narrow cone instead of wastefully tracking a
    //     wide-angle direction that can never trigger. Position is kept either way.
    static const G4double kHalfPi = std::acos(-1.0) / 2.0;
    static const G4double kCollimatorRad = [] {
        const char* e = std::getenv("COLLIMATOR_DEG");
        return (e ? std::atof(e) : 8.0) * (std::acos(-1.0) / 180.0);
    }();

    const bool isOpen = (DetectorConstruction::GetMode() == DetectorConstruction::Mode::kImagingOpen);
    const G4double srcHX = fSrcHX;
    const G4double E     = SampleEnergy();

    const G4double xs = (2.*G4UniformRand() - 1.) * srcHX;
    const G4double ys = (2.*G4UniformRand() - 1.) * srcHX;
    const G4ThreeVector pos(xs, ys, fZsrc);

    G4ThreeVector dir;
    if (isOpen) {
        dir = SampleNarrowConeDir(kCollimatorRad);
    } else {
        int guard = 0; bool ok = false;
        do {
            const G4double a   = SamplePitch() * deg;
            const G4double phi = (2.*G4UniformRand() - 1.) * kHalfPi;   // forward hemisphere
            dir.set(std::sin(a)*std::sin(phi), std::cos(a), std::sin(a)*std::cos(phi));
            ok = (G4UniformRand() <= dir.z());   // cos-incidence accept
        } while (!ok && ++guard < 100000);
        if (!ok) return;

        if (!WouldHitUranium(pos, dir)) dir = SampleNarrowConeDir(kCollimatorRad);
    }

    fGun->SetParticlePosition(pos);
    fGun->SetParticleMomentumDirection(dir.unit());
    fGun->SetParticleEnergy(E * MeV);
    fGun->GeneratePrimaryVertex(event);
    DumpBeamSample(xs, ys, fZsrc, E, dir.unit());
}

// Prescaled dump of every primary that is actually launched (after the cone cut and the
// cos-incidence rejection), so the dumped set is an unbiased 1-in-N sample of the beam
// as fired. The trigger CSV's birth_* columns cover the same quantities but only for
// the near-axis subset that made the quadruple coincidence.
long PrimaryGeneratorAction::BeamDumpPrescale()
{
    static const long kN = [] {
        const char* e = std::getenv("BEAM_DUMP_PRESCALE");
        return e ? std::atol(e) : 100000L;
    }();
    return kN;
}

void PrimaryGeneratorAction::DumpBeamSample(G4double xs, G4double ys, G4double zs,
                                            G4double E_MeV, const G4ThreeVector& dir) const
{
    const long presc = BeamDumpPrescale();
    if (presc <= 0) return;

    // Per-thread counter, no atomics or shared state; the union over threads is still
    // an unbiased 1-in-N sample of the fired beam.
    static G4ThreadLocal long tCount = 0;
    if (++tCount % presc != 0) return;

    // B is along +y, so the pitch angle is the angle from +y and the gyrophase is
    // measured in the x-z plane -- exactly the (alpha, phi) the sampler drew from:
    //   dir = (sin(alpha) sin(phi), cos(alpha), sin(alpha) cos(phi))
    const G4double pitchDeg = std::acos(std::max(-1.0, std::min(1.0, dir.y()))) / deg;
    const G4double gyroDeg  = std::atan2(dir.x(), dir.z()) / deg;

    G4AutoLock lock(&RunAction::GetBeamMutex());
    auto* out = RunAction::GetBeamFile();
    if (!out || !out->is_open()) return;
    *out << xs/mm << ',' << ys/mm << ',' << zs/mm << ','
         << E_MeV << ','
         << dir.x() << ',' << dir.y() << ',' << dir.z() << ','
         << pitchDeg << ',' << gyroDeg << '\n';
}

void PrimaryGeneratorAction::GenerateInvertedBeam(G4Event* event)
{
    // Uniform-flux sampling: pick the desired landing position on the SSD face first,
    // then back-compute the aim point for this proton's energy-dependent bending. Every
    // fired proton lands on the detector face, so this is ~100% efficient but cannot
    // sample uranium scatter-in (protons that would miss but scatter in) -- that's why
    // the default beam above is the physical forward-hemisphere flux instead.
    const G4double kActHX  = DetectorConstruction::GetActiveHalfX_mm() * mm;
    const G4double kActHY  = DetectorConstruction::GetActiveHalfY_mm() * mm;
    static constexpr G4double kMargin =  100.0 * mm;
    static constexpr G4double kB_T    = 0.1e-4;
    static constexpr G4double kMp     = 938.272;
    const G4double kd = fZdet - fZsrc;  // drift [mm], derived from standoff

    // Sample energy; compute this proton's bending.
    const G4double E       = SampleEnergy();
    const G4double p       = std::sqrt(E * (E + 2.*kMp));
    const G4double r_g     = p / (299.792 * kB_T);            // m
    const G4double dx_bend = -(kd * kd) / (2000.0 * r_g);    // mm  (kd in mm, r_g in m)

    // Sample desired landing position in local detector frame (centred at 0).
    const G4double xl_loc = (2.*G4UniformRand() - 1.) * (kActHX + kMargin);
    const G4double yl     = (2.*G4UniformRand() - 1.) * (kActHY + kMargin);

    // Aim point: shift by (mean_bend - this_bend) so the proton lands at xl_loc
    // relative to the detector centre (which sits at world x = fMeanDxBend_mm).
    const G4double xd = (fMeanDxBend_mm + xl_loc) - dx_bend;
    const G4double yd = yl;

    // Source position: uniform +-fSrcHX, gives angular diversity.
    const G4double xs = (2.*G4UniformRand() - 1.) * fSrcHX;
    const G4double ys = (2.*G4UniformRand() - 1.) * fSrcHX;

    G4ThreeVector dir(xd - xs, yd - ys, fZdet - fZsrc);
    dir = dir.unit();

    fGun->SetParticlePosition(G4ThreeVector(xs, ys, fZsrc));
    fGun->SetParticleMomentumDirection(dir);
    fGun->SetParticleEnergy(E * MeV);
    fGun->GeneratePrimaryVertex(event);
    DumpBeamSample(xs, ys, fZsrc, E, dir.unit());
}

void PrimaryGeneratorAction::GenerateWallTestBeam(G4Event* event)
{
    // Source plane: parallel to the y-z plane (normal along x), positioned just outside
    // the collimator's own outer wall (+x side), spanning its FULL z-length (front graded
    // zone + back side-shield zone) and the same +-outerHalfXY in y. Direction: the real
    // PAD (SamplePitch/LoadPAD, full 0-2pi gyrophase -- illumination can arrive from any
    // side), accepted with probability -dir.x() (cos-incidence onto this face, inward =
    // -x) -- the general form of the SAME cos-incidence fluence weighting the normal beam
    // already uses relative to +z. No WouldHitUranium, no redirect: the sampled direction
    // IS the final primary, so Geant4's real tracking (absorption in the graded wall, or
    // punch-through to the tracker) is what determines the outcome, not a source-term
    // shortcut. By the phi-uniform symmetry of the PAD sampling, the +x face is a fair
    // representative sample of any of the tube's 4 sides.
    static const G4double zFront = DetectorConstruction::GetCollZFront_mm() * mm;
    static const G4double zBack  = DetectorConstruction::GetCollZBack_mm()  * mm;
    static const G4double outerHalfXY = DetectorConstruction::GetCollOuterHalfXY_mm() * mm;
    static const G4double margin = 20.0 * mm;
    static const G4double X0 = fMeanDxBend_mm + outerHalfXY + margin;

    const G4double E = SampleEnergy();

    G4ThreeVector dir;
    bool ok = false;
    int guard = 0;
    do {
        const G4double a   = SamplePitch() * deg;
        const G4double phi = 2.0 * std::acos(-1.0) * G4UniformRand();   // full gyrophase
        dir.set(std::sin(a)*std::sin(phi), std::cos(a), std::sin(a)*std::cos(phi));
        ok = (G4UniformRand() <= -dir.x());   // cos-incidence onto the +x face, inward -x
    } while (!ok && ++guard < 100000);
    if (!ok) return;

    const G4double y = (2.*G4UniformRand() - 1.) * outerHalfXY;
    const G4double z = zFront + G4UniformRand() * (zBack - zFront);
    const G4ThreeVector pos(X0, y, z);

    fGun->SetParticlePosition(pos);
    fGun->SetParticleMomentumDirection(dir.unit());
    fGun->SetParticleEnergy(E * MeV);
    fGun->GeneratePrimaryVertex(event);
    DumpBeamSample(pos.x(), pos.y(), pos.z(), E, dir.unit());
}
