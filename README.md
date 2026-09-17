# proton_tracker_2layer10cm

Geant4 Monte Carlo for OST proton radiography: a 2-Si-tracking-layer variant (layers
100 mm apart, vs. the main project's 3-layer/20 mm design) with an Al tube-taper
collimator (8 deg design angle) and an X-Y strip hodoscope trigger stage. This repo is
**source only** — no build output, no simulation results. It was pushed from a Windows
dev machine specifically to be built and run on MIT ORCD HPC.

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

## Before running on this cluster — two things to check first

1. **Thread count is hardcoded.** `MAIN.cc` line ~55:
   `mtRunManager->SetNumberOfThreads(18);` — set this to match the SLURM allocation
   (`--cpus-per-task`) before building, or the run will over/under-subscribe the node.
2. **Standoff is hardcoded and requires a rebuild to change.** `src/DetectorConstruction.cc`,
   `ImagingDet::kStandoff` (currently checked in at `400.0 * m`) is the single variable
   that sets the standoff — edit it, then `cmake --build build -j` again, for each
   standoff you want (this project's convention: 100 m, 200 m, 400 m).

## Running

Mode and standoff are both picked up from the macro filename, but **the macro names in
this project are a legacy naming quirk — the number in the filename is NOT the standoff
in meters.** `imaging_10m_v3.mac` / `imaging_20m_v3.mac` / `imaging_40m_v3.mac` correspond
to standoffs **100 m / 200 m / 400 m** respectively (and must be run against a build with
`kStandoff` set to match). Each has an `imaging_open_*` counterpart for the open-field
(no-uranium) reference run.

```bash
cd build
./MAIN imaging_10m_v3.mac        # shadow run  (standoff = whatever kStandoff was built with)
./MAIN imaging_open_10m_v3.mac   # open-field reference, same standoff
```

Output (written to the current directory): `imaging_shadow.csv` + `imaging_shadow_stats.csv`
(shadow run), `imaging_open.csv` + `imaging_open_stats.csv` (open run). `*_stats.csv` holds
a single `n_fired` column (total primaries fired) — needed to combine/merge chunks.

**Chunking / seeding for independent runs:** set the `SEED` env var to a distinct value
per chunk (`SEED=12345 ./MAIN imaging_10m_v3.mac`) for reproducible, statistically
independent chunks that can be concatenated (CSV rows) and summed (`n_fired`) afterward.
Omit `SEED` to seed off the clock instead.

## Geometry notes

- 2 Si tracking layers, 100 mm apart (`Det::` namespace in `DetectorConstruction.cc`) —
  chosen so the pixel-resolution term (`sigma_pos`) is subdominant to multiple scattering;
  see the comment above `LAYER_GAP_MM` in the (not-included) MATLAB analysis script for the
  derivation if you need it.
- Al tube-taper collimator, `CollTube::kDesignDeg = 8.0` (deg), 5-zone mass-optimized
  taper (see `DetectorConstruction.cc`'s `ZoneSpec` table in the `kCollStyle==kTube`
  branch).
- X-Y strip hodoscope trigger stage (1 cm plastic scintillator strips) instead of the
  main project's poly-filter + single scintillator slab.
- `COLLIMATOR_DEG` env var (default matches `CollTube::kDesignDeg` in source, currently
  8) sets the source-term sampling half-angle for `PrimaryGeneratorAction` — this is a
  *sampling* efficiency knob, not the collimator's physical acceptance, which is set by
  the real geometry above.
