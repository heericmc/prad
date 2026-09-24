# proton_tracker_2layer10cm

Geant4 Monte Carlo for OST proton radiography: a 2-Si-tracking-layer variant (layers
100 mm apart, vs. the main project's 3-layer/20 mm design) with an Al tube-taper
collimator (11.7 deg design angle) and an X-Y strip hodoscope trigger stage. This repo is
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

**Current goal (2026-09-24): re-run every standoff (100/200/400/800/1600 m) with the new
counted sampling scheme below.** Every shadow run taken before 2026-09-24 under-weights
uranium scatter-in ~10x; its T is biased high and its t_5sigma biased low.

## 2026-09-24: counted sampling scheme (`SAMPLING` env var) — READ FIRST

### What was wrong with the old (legacy) shadow beam

The legacy shadow run draws each proton from the full forward-hemisphere PAD. If the
proton would hit the uranium it keeps that direction, which is correct. If it would miss,
its direction is **redrawn into the 11.7 deg cone and it is fired anyway**. The analysis
then treats every shadow event as one open-run event, i.e. one proton from the cone.
The event count therefore follows the cone flux, but which events hit the uranium was
decided with the hemisphere flux, which is Phi_hemi / (J_3D(90) * Omega_cone) = **10.0x**
larger. So the uranium sees only 1/10 of its real hits, and so do scatter-in and
punch-through. The direct (non-uranium) flux is correct, and so is the whole open run.

Measured at 100 m against a brute-force run (`SAMPLING=all`, every proton fired at its
true angle, nothing redirected or skipped):

| per real second | legacy archive | brute force | new scheme |
|---|---|---|---|
| all shadow triggers | 205.9 +- 1.5 | 211.0 +- 4.8 | 215.8 +- 1.0 |
| scatter-in (launch > 1.5 deg) | **0.76 +- 0.09** | 7.6 +- 0.9 | 6.2 +- 0.8 |
| scatter-in in the +-229 mm cell | **0.27 +- 0.05** | 3.3 +- 0.6 | 1.9 +- 0.4 |
| open field, all triggers | 232.6 +- 0.5 | 235.0 +- 5.0 | 233.9 +- 0.6 |

Consequence at 100 m: T in the cell goes 0.476 -> 0.555 +- 0.015, t_5sigma 2.9 -> 4.0 s.
This also explains the old 0.886 shadow/open global-efficiency warning at 100 m: ~11% of
legacy shadow events are wide-angle uranium hits that almost never trigger.

### The new scheme

Every candidate is drawn from ONE physical parent distribution: uniform position on the
+-src_hx source plane, pitch from J(alpha)*sin(alpha), gyrophase over the forward
hemisphere, cos-incidence accepted, AP9 energy (E >= 200 MeV). E, position and direction
are all redrawn per candidate. A candidate that cannot matter is **counted and skipped,
never redirected**. Each run counts all candidates (`n_cand`) and those inside a 1 deg
reference cone about +z (`n_ref`; pitch 90 deg, where the PAD is flat), which gives its
real-time equivalent from the trusted J_3D(90 deg) = 160.09 with no flux integral:

    t_real = n_ref / (160.09 * pi*sin^2(ref_cone_deg) * (2*src_hx_mm/10)^2)   [s]

| `SAMPLING=` | fires | use |
|---|---|---|
| `legacy` (default) | the old beam, unchanged (also counts n_cand/n_ref) | reproducing old archives only |
| `direct` | candidates that do NOT hit the uranium and whose bent path reaches the collimator+tracker envelope (+100 mm margin, `DIRECT_MARGIN_MM`); drawn from a narrow pre-cone (`DIRECT_CONE_DEG`, default 2x the geometric reach) for speed | **production, open AND shadow** |
| `uranium` | candidates that hit the uranium, at their true (any) angle, full hemisphere | **production, shadow only** (fatal in an open run) |
| `split` | both of the above in one job, full hemisphere parent | small tests |
| `all` | every candidate (brute force) | validation only |

A shadow measurement is therefore TWO independent jobs, `direct` and `uranium`, each with
its own `t_real`; the analysis weights each by t_open / t_own. Cost on a 20-thread laptop
at 100 m: `direct` ~1000x faster per real second than brute force (e.g. 400 m: 5000 events
= 8 s of real time); `uranium` ~0.5 core-hours per real second at ANY standoff (the
uranium-hit rate per real second does not depend on standoff), i.e. it dominates the cost.
Its contribution is ~10% of the shadow signal, so it can be run at a shorter real time
than the direct part (e.g. t_uranium ~ t_direct / 3) at little cost in precision.

### Stats file (every run)

`imaging_{shadow,open}_stats.csv` is now
`n_fired,n_cand,n_ref,ref_cone_deg,src_hx_mm,sampling,seed` (n_fired stays first, so old
readers still work). Counts are exact 64-bit integers.

### Running a counted campaign at one standoff

```bash
cd build
SEED=<unique> SAMPLING=direct  ./MAIN imaging_open_<N/10>m_v3.mac   # open field
SEED=<unique> SAMPLING=direct  ./MAIN imaging_<N/10>m_v3.mac        # shadow, direct part
SEED=<unique> SAMPLING=uranium ./MAIN imaging_<N/10>m_v3.mac        # shadow, uranium part
```

Each writes the usual `imaging_open*.csv` / `imaging_shadow*.csv`; rename the uranium
job's outputs to `imaging_shadow_uranium.csv` / `imaging_shadow_uranium_stats.csv`. The
event counts in the old macros are sized for the legacy beam and are FAR too many for
`direct` (nearly every `direct` event is a trigger candidate), so set `/run/beamOn` from a
calibration chunk: pick the real time you want, then events = t_real_target / (t_real per
event measured in the calibration).

**Archive layout** (what the analysis expects), per standoff:
`run_<N>m_2layer10cm_v3_collTubeTaper_11p7deg_counted/` holding
`imaging_open.csv`, `imaging_open_stats.csv` (direct), `imaging_shadow.csv`,
`imaging_shadow_stats.csv` (direct), `imaging_shadow_uranium.csv`,
`imaging_shadow_uranium_stats.csv` (uranium).

**Merging chunks:** concatenate CSV rows (one header), and in the stats file **sum
`n_fired`, `n_cand` and `n_ref`** column by column; `ref_cone_deg`, `src_hx_mm` and
`sampling` must be identical across chunks (refuse to merge otherwise). Never merge
`direct` and `uranium` chunks together — they are different strata with different real
times. **Before merging, check that every chunk's `seed` is distinct and that no two chunk
CSVs are byte-identical (md5)**: in the 2026-09 HPC campaign all 64 shadow chunks at 100 m
and at 1600 m were byte-identical (same seed), so those runs held only one chunk's worth
of statistics.

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
   the standoff you are about to run** (one build per standoff) before building — this is the single variable that sets the standoff;
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

## Target: 800 m — how many particles (SUPERSEDED 2026-09-24: legacy-beam campaign, kept for history)

> The event counts below are for the legacy beam and its shadow runs under-weight
> scatter-in ~10x. For new runs use the counted scheme above; size `/run/beamOn` from a
> calibration chunk instead.

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
- Al tube-taper collimator, `CollTube::kDesignDeg = 11.7` (deg), 20-zone mass-optimized
  taper in 200 mm steps, thickness linearly interpolated between the original 8 deg
  design's validated anchor points (see `DetectorConstruction.cc`'s `thickFracAt`
  lambda in the `kCollStyle==kTube` branch — zone count/lengths are derived from the
  built length, not hardcoded). Bore ±562 mm, total length ≈3.99 m, mass ≈**3.07 t**.
  (An earlier same-day pass re-segmented to uniform 1000 mm zones instead, which came
  out to 4.18 t — *heavier* than the retired 8 deg/5-zone design at 4.06 t, because
  holding the near-detector 100%-thickness zone at a full 1000 mm instead of the
  original 500 mm cost more mass than the shorter bore saved. The 200 mm/interpolated
  taper fixes this: 27% lighter than that 1000 mm version and 24% lighter than the
  original 8 deg design, despite the wider FOV.)
  **2026-09-17: widened from 8 deg** — the 8 deg point was originally justified by a
  two-plane track-confusion probability calc that used the wrong plane spacing
  (2 cm instead of this repo's actual 10 cm); re-derived correctly, bare two-plane
  confusion needs ≈3.56 deg to meet the 10⁻³ ceiling on its own, but the *actual* built
  design also requires the X-Y strip hodoscope to register a **position-matched** hit
  (not just "some strip fired") as a required third confirmation, and that compound
  probability meets the same ceiling out to ≈11.7 deg — see
  `maps-pileup-collimation.md` §3a for the full derivation. **The strip hodoscope's
  position match is therefore load-bearing for this design, not optional**: the
  analysis pipeline must check that the fired `stripx_id`/`stripy_id` cell is
  consistent with the extrapolated plane0→plane1 track, not merely require a
  hodoscope coincidence. **The re-segmented taper has not yet been re-validated with a
  fresh `WALL_TEST=1` run** — do this before trusting it for a production campaign.
- X-Y strip hodoscope trigger stage (1 cm plastic scintillator strips) instead of the
  main project's poly-filter + single scintillator slab.
- `COLLIMATOR_DEG` env var (default matches `CollTube::kDesignDeg` in source, currently
  11.7) sets the source-term sampling half-angle for `PrimaryGeneratorAction` — this is
  a *sampling* efficiency knob, not the collimator's physical acceptance, which is set
  by the real geometry above.
