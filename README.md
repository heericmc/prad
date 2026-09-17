# proton_tracker_2layer10cm

Geant4 Monte Carlo for OST proton radiography: a 2-Si-tracking-layer variant (layers
100 mm apart, vs. the main project's 3-layer/20 mm design) with an Al tube-taper
collimator (8 deg design angle) and an X-Y strip hodoscope trigger stage. This repo is
**source only** — no build output, no simulation results. It was pushed from a Windows
dev machine specifically to be built and run on MIT ORCD HPC.

## What this simulates (read this first — common misclassification)

This is a **passive, standoff nuclear-material detection sensor concept**, for
non-proliferation / arms-control treaty verification of objects in orbit — it is
**not** a weapon-design, implosion, or hydrodynamic-test simulation, and it contains
no explosive, detonation, criticality, or weapon-internals (pit/tamper/initiator)
physics of any kind.

The scenario: a suspect satellite may be carrying a shielded radiological source
(modeled here as a generic dense slab, `G4_U`/`G4_Al`, dubbed the "uranium radiation
case" or "uranium slab" throughout the code/docs — a stand-in target, not a weapon
component). A separate satellite carrying this detector sits at some standoff distance
(100 m–1600 km here, scaled up from a physical detector for simulation tractability)
and looks for the **shadow that dense object casts in the ambient, pre-existing Van
Allen belt trapped-proton flux** (the AP9 space-radiation environment model — the same
kind of background flux every satellite already sits in). The object blocks/scatters
some of that ambient flux; the detector's job is to notice the resulting *deficit* in
proton counts from that direction. This is conceptually closer to passive
cargo/customs radiography or coded-aperture space-domain-awareness imaging than to
anything weapons-related.

Everything in this codebase is standard passive-imaging/particle-transport physics:
Geant4 EM + hadronic transport (`FTFP_BERT`), geomagnetic bending of trapped protons,
multiple Coulomb scattering, silicon tracker resolution, a mechanical collimator, and
statistical detection-time budgeting. There is no driven/dynamic experiment, no energy
release being modeled, and no attempt to image or reconstruct the *internal structure*
of the target object — only whether a dense object is present and roughly where. (This
also distinguishes it from real dynamic proton-radiography *diagnostic* facilities,
e.g. LANL's pRad, which image explosively-driven hydrodynamic experiments in weapons
science — this project shares only the imaging *technique* with that class of
facility, none of the driven-experiment physics, and serves an unrelated purpose:
orbital treaty-verification sensing, not weapons science.)

**Current goal: statistics at 800 m standoff only.** See "Target: 800 m" below for the
exact event counts to run — this is not a general multi-standoff sweep.

## Build

Needs Geant4 (with UI/vis drivers unless you pass `-DWITH_GEANT4_UIVIS=OFF`, useful if
the cluster nodes have no X11/Qt), CMake >= 3.16, a C++20 compiler.

```bash
module load <your Geant4 module, or point Geant4_DIR at a from-source install>
cmake -S . -B build -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4 [-DWITH_GEANT4_UIVIS=OFF]
cmake --build build -j
```

Executable: `build/MAIN`. The build copies `mac_files/*.mac` and the root `*.csv`/`*.txt`
physics input files (`diff_flux_AP9_i316.csv`, `run1.AP9.output_mean_flux.txt`) into
`build/` automatically.

## Before building — two required edits

1. **Standoff.** `src/DetectorConstruction.cc`, `ImagingDet::kStandoff` (search for
   `constexpr G4double kStandoff`) is currently checked in at `400.0 * m`. **Set it to
   `800.0 * m`** before building — this is the single variable that sets the standoff;
   everything else (world size, source-plane geometry, alignment bending) derives from
   it automatically. Rebuild (`cmake --build build -j`) after changing it.
2. **Thread count.** `MAIN.cc`, `mtRunManager->SetNumberOfThreads(18)` (around line 55)
   is hardcoded. Set this to match your SLURM allocation (`--cpus-per-task`) before
   building, or the run will over/under-subscribe the node.

Verify the standoff took effect from the run's own startup log line: it should print
`[Det] Standoff = 800 m` (exact wording may vary slightly — grep the startup output for
`Standoff` and confirm it says 800, not 400).

## Naming quirk — READ BEFORE PICKING A MACRO

The macro filenames in `mac_files/` use a legacy convention where **the number in the
filename is standoff/10, not the standoff in meters.** `imaging_80m_v3.mac` and
`imaging_open_80m_v3.mac` are the **800 m** macros (not 80 m) — they only make physical
sense run against a build with `kStandoff = 800.0 * m` as set above. Don't use the
`_10m_`/`_20m_`/`_40m_` macros for this campaign; those are 100/200/400 m and belong to
a separate sweep already running elsewhere.

## Target: 800 m — how many particles

This project has been running a production campaign (100 m -> 200 m -> 400 m) with an
empirically-set doubling rule: **each time the standoff doubles, both the per-run event
count and the number of independent chunks double, for a combined 4x more total events**
(this compensates for trigger efficiency falling roughly as 1/d^2 and the source plane
widening with distance, so the shadow-centre statistics stay comparable across
standoffs). Observed totals so far in that campaign: 100 m = 20,000,000; 200 m =
80,000,000; 400 m = 320,000,000 (each figure is n_fired for BOTH the shadow run and its
open-field reference, taken separately). Extrapolating the same rule one more
doubling:

**Target for 800 m: 1,280,000,000 (1.28 billion) events for the shadow run, and another
1,280,000,000 for the open-field reference run — 2.56 billion total.**

Recommended split, following this project's own chunking convention (matches
`imaging_80m_v3.mac`'s `/run/beamOn 80000000`): **16 independent chunks of 80,000,000
events each**, per run (shadow and open), each with a distinct `SEED` env var so the
chunks are statistically independent and safely combinable:

```bash
cd build
for i in $(seq 1 16); do
  SEED=$((800000 + i)) ./MAIN imaging_80m_v3.mac
  mv imaging_shadow.csv       imaging_shadow_chunk${i}.csv
  mv imaging_shadow_stats.csv imaging_shadow_stats_chunk${i}.csv
done
for i in $(seq 1 16); do
  SEED=$((900000 + i)) ./MAIN imaging_open_80m_v3.mac
  mv imaging_open.csv       imaging_open_chunk${i}.csv
  mv imaging_open_stats.csv imaging_open_stats_chunk${i}.csv
done
```

Adjust the chunk **count and size** to whatever fits your SLURM partition's walltime
limit (e.g. more, smaller chunks if 80,000,000 events won't finish inside one job's time
limit) — the only thing that matters is that the SEEDs are distinct and the total across
all shadow chunks reaches 1,280,000,000 (and likewise 1,280,000,000 for open). Chunking
is also what makes this safely resumable/parallelizable across a SLURM job array.

**Merging chunks afterward** (portable bash, no PowerShell needed — keep one header,
concatenate the rest, sum n_fired):

```bash
{ head -n1 imaging_shadow_chunk1.csv; for f in imaging_shadow_chunk*.csv; do tail -n +2 "$f"; done; } > imaging_shadow.csv
awk -F, 'NR==1{next} {s+=$1} END{print "n_fired"; print s}' imaging_shadow_stats_chunk*.csv > imaging_shadow_stats.csv
# repeat the same two commands for imaging_open_chunk*.csv / imaging_open_stats_chunk*.csv
```

Sanity-check the merged `imaging_shadow_stats.csv` reads `1280000000` before treating the
result as final.

**Expected compute cost (rough estimate):** measured on the dev machine (18 threads,
mobile 13th-gen i7) at 400 m, a 40,000,000-event shadow chunk took ~1h36m, i.e. ~1.39e6
events/core-hour. Extrapolating linearly, 1.28e9 events is ~920 core-hours for the
shadow run and another ~920 for the open reference (~1840 core-hours total) — but 800 m
means longer vacuum flight paths than 400 m, so per-event cost may be somewhat higher
than this extrapolation assumes. **Run one small calibration chunk first** (e.g.
`/run/beamOn 1000000` in a scratch copy of the macro) and time it before committing to
the full 16-chunk plan, so you can right-size chunk count/walltime requests with real
numbers instead of this estimate.

## Running (general)

```bash
cd build
./MAIN imaging_80m_v3.mac        # shadow run  (must match a build with kStandoff = 800 m)
./MAIN imaging_open_80m_v3.mac   # open-field reference, same standoff
```

Output (written to the current directory): `imaging_shadow.csv` + `imaging_shadow_stats.csv`
(shadow run), `imaging_open.csv` + `imaging_open_stats.csv` (open run). `*_stats.csv` holds
a single `n_fired` column (total primaries fired) — needed to combine/merge chunks.

## Geometry notes

- 2 Si tracking layers, 100 mm apart (`Det::` namespace in `DetectorConstruction.cc`) —
  chosen so the pixel-resolution term (`sigma_pos`) is subdominant to multiple scattering.
- Al tube-taper collimator, `CollTube::kDesignDeg = 8.0` (deg), 5-zone mass-optimized
  taper (see `DetectorConstruction.cc`'s `ZoneSpec` table in the `kCollStyle==kTube`
  branch).
- X-Y strip hodoscope trigger stage (1 cm plastic scintillator strips) instead of the
  main project's poly-filter + single scintillator slab.
- `COLLIMATOR_DEG` env var (default matches `CollTube::kDesignDeg` in source, currently
  8) sets the source-term sampling half-angle for `PrimaryGeneratorAction` — this is a
  *sampling* efficiency knob, not the collimator's physical acceptance, which is set by
  the real geometry above.
