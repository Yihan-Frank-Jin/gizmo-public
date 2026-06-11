# Repository Guidelines

## Project Overview

This repository is the public GIZMO codebase: a Makefile-built, MPI-parallel C astrophysics simulation code descended from GADGET. It supports meshless finite-volume/mass hydrodynamics, SPH, fixed-grid Godunov modes, self-gravity with tree/PM solvers, cosmological integration, cooling/chemistry, galaxy formation, black holes, radiation transport, solids, dust, turbulence, and dark matter extensions.

The main executable is `GIZMO`. The entry point is `main.c`, which initializes MPI, calls `begrun()` for setup and parameter parsing, then enters the main timestep loop in `run.c`.

## Current Local Context

The active research task is **updating AGN feedback in GIZMO with MACER3D's model**. In practice, this means porting the AGN feedback model from `scripts/yuan18.cpp` into GIZMO's black hole infrastructure under `galaxy_sf/blackholes/`. The main references are `scripts/yuan18.cpp`, `scripts/Yuan_2018_ApJ_857_121.pdf`, and `scripts/Zhang_MACER3D.pdf`.

The port has three feedback components: wind, radiation, and jet. The accretion portion is currently completed. The active implementation target is the wind feedback component. Radiation and jet feedback may be developed as parallel tasks when useful, but they should not distract from or silently change the wind-feedback implementation unless the user explicitly redirects the task.

Important local compile flags in the current `Config.sh` include:

- `BOX_PERIODIC`
- `HYDRO_MESHLESS_FINITE_MASS`
- `EOS_GAMMA=(5.0/3.0)`
- `ADAPTIVE_GRAVSOFT_FORALL=0+5`
- `BLACK_HOLES`
- `BH_SWALLOWGAS`
- `BH_YUAN18_ACCRETION`
- `BH_YUAN18_WIND`
- `BH_OUTPUT_MOREINFO`
- `OUTPUT_ADDITIONAL_RUNINFO`
- `USE_FFTW3`

Treat `Config.sh`, `Makefile.systype`, built objects, executables, snapshots, and run outputs as machine/local state unless the user explicitly asks to commit or standardize them. `Makefile.systype` currently selects `SYSTYPE="MacBookCellar"`.

## Project Structure

Core integration, orchestration, I/O, and global state live at the repository root:

- `main.c` starts MPI and dispatches setup/run.
- `begrun.c` reads runtime parameter files into the global `All` struct.
- `run.c` owns the main timestep loop.
- `allvars.h` defines global state, particle data structures, compile-time macro logic, and the `All` runtime parameter struct.
- `proto.h` is the central function prototype header.
- `io.c`, `read_ic.c`, and `restart.c` handle snapshots, ICs, and restarts.
- `domain.c`, `ngb.c`, `timestep.c`, `kicks.c`, `predict.c`, `driftfac.c`, and `merge_split.c` implement decomposition, neighbor search, time integration, and particle maintenance.

Physics and infrastructure modules are grouped by directory:

- `gravity/`: tree/FMM gravity, PM gravity, cosmology, potentials, adaptive softening.
- `hydro/`: hydro driver, density, gradients, meshless/SPH headers, Riemann interface, conduction/viscosity/MHD headers.
- `cooling/`: cooling, Grackle interface, simple chemistry, CHIMES support.
- `eos/`: equations of state, hydrogen molecule code, cosmic ray fluid, Helmholtz EOS support.
- `galaxy_sf/`: star formation and stellar feedback.
- `galaxy_sf/blackholes/`: black hole accretion, feedback, swallowing/kicking/spawning, Yuan18 work.
- `radiation/`: radiation hydrodynamics utilities, chemistry, cooling, source injection, CG methods.
- `structure/`: FOF, SubFind, line-of-sight, two-point analysis.
- `turb/`: turbulent diffusion, turbulent driving, spectra.
- `solids/`: elastic, grain, and dust chemistry physics.
- `sidm/`: self-interacting, fuzzy, and scalar-field dark matter.
- `nuclear/`: nuclear reaction network code.
- `system/`: custom allocation, MPI utilities, sorting, Peano-Hilbert ordering, pinning, vector helpers.
- `scripts/`: docs, test problem parameter files, IC helpers, analysis helpers, and Yuan18/MACER references.

`blackhole_run/` contains a local Bondi/Yuan18 test setup, scripts, snapshots, logs, and comparison output. Do not treat it as source unless the user asks. `XPHO_PROJECT/` is a separate nested project with its own `.git`; avoid touching it for GIZMO work unless explicitly requested.

## Scripts Directory Notes

The `scripts/` directory is a mix of documentation, reference physics, example parameter files, analysis utilities, and destructive archival helpers. Do not treat every script as current source behavior; check whether it is a reference, a template, a validation helper, or an archival tool before using it.

Important files and how to use them:

- `scripts/gizmo_documentation.md` is the main user guide. It is the best local reference for feature/module meaning, public vs. proprietary restrictions, compile-time flags, parameter-file semantics, snapshot fields, units, common run failures, and test-problem descriptions.
- `scripts/params.txt` is a broad runtime parameter template, not a recommended production parameter file. Use it to discover parameter names and defaults, then use problem-specific parameter files or source parsing to confirm what a run really needs.
- `scripts/test_problems/*.params` are focused validation setups. Their leading comments often list required `Config.sh` flags and special source edits, such as enabling analytic gravity calls for Keplerian or Rayleigh-Taylor tests. Do not assume these examples can be mixed blindly.
- `scripts/make_IC.py` demonstrates HDF5 IC layout: `Header` attributes must match the actual particle counts; `PartType0` is gas and needs fields such as `InternalEnergy`; other particle types are collisionless unless compile-time flags give them special meaning; `MassTable` should be zero when masses are stored per particle; including a dummy `MagneticField` block is harmless when `MAGNETIC` is off.
- `scripts/load_from_snapshot.py` is the flexible HDF5 reader and can list or load arbitrary dataset keys. Prefer it for new physics outputs. `scripts/readsnap.py` maps common GIZMO/GADGET fields into short names and can handle binary outputs, but it must be updated manually for new output blocks.
- `scripts/compress_gizmosnap.py` losslessly rewrites HDF5 snapshots in place with gzip/shuffle/fletcher32 and marks `Header.attrs['CompactLevel']`. `scripts/compress_cleanup_sim.sh` deletes logs/restarts and compresses diagnostics only when a `clean_me` marker exists, but it still uses broad `rm` commands. Do not run either on active or user-owned outputs without explicit permission.
- `scripts/visit/` contains an experimental VisIt HDF5 plugin and index-file helper. Treat it as optional visualization infrastructure, not core simulation code.
- `scripts/indent.sh` wraps the repository's GNU-style `indent` settings.
- `scripts/check_yuan18_hot_rtr_split.py` is a lightweight regression check for the Yuan18 hot-mode radius split in `galaxy_sf/blackholes/blackhole.c`; run it after editing that logic.

## Build System

Build configuration is split across:

- `Template_Config.sh`: canonical list of compile-time physics flags.
- `Config.sh`: local active compile-time configuration.
- `Makefile`: object lists, compiler/linker logic, machine blocks, generated config processing.
- `Makefile.systype`: local machine selection.
- `config-makefile` and `prepare-config.perl`: generate compile-time headers such as `GIZMO_config.h`.

Common commands:

```sh
cp Template_Config.sh Config.sh
make
make CONFIG=MyConfig.sh EXEC=GIZMO_test
make clean
mpirun -np <N> ./GIZMO scripts/test_problems/shocktube.params
```

Run `make clean; make` after changing compile-time flags. New flags should be added disabled/defaulted in `Template_Config.sh`; avoid committing local-only `Config.sh` changes unless they are intentional and reproducible.

The Makefile always includes most object groups and relies heavily on compile-time guards inside the source. Special external-library cases include Helmholtz EOS, CHIMES, Grackle, HDF5, and FFTW.

## Main Runtime Flow

The main loop in `run.c` follows this shape:

1. Initial domain decomposition, physics setup, gravity, hydro, and source terms.
2. Per timestep: statistics and CPU logs.
3. Select timesteps and perform the first half-kick.
4. Drift to the next sync point and write scheduled outputs.
5. Rebuild/update domain and tree structures as needed.
6. Compute gravity.
7. Mark feedback sources.
8. Compute hydro densities, gradients, and forces.
9. Optionally merge/split/rearrange particles.
10. Apply the second half-kick.
11. Apply non-standard physics/source terms.
12. Stop, restart, or continue until `TimeMax`.

Source terms and feedback hooks are generally handled through `calculate_non_standard_physics()` and module-specific calls guarded by compile-time flags.

## Yuan18 / MACER AGN Work

The local branch contains substantial Yuan18/MACER black hole accretion and wind work in `galaxy_sf/blackholes/`. The goal is not to invent a new AGN model; it is to reproduce the MACER3D/Yuan18 model behavior inside GIZMO's particle-based architecture.

Reference implementation and notes:

- `scripts/yuan18.cpp`: Athena++/MACER3D reference source.
- `scripts/Yuan_2018_ApJ_857_121.pdf`: Yuan et al. 2018 model paper.
- `scripts/Zhang_MACER3D.pdf`: MACER3D paper.
- `scripts/yuan18_feedback_plan.md` and `.claude/*.md`: local design notes and prior analysis. `scripts/yuan18_feedback_plan.md` is partly historical: it describes a smooth/kernel wind-injection plan, while the current local implementation is spawn-based.

Reference-model behavior learned from `scripts/yuan18.cpp`:

- The file is Athena++ problem/reference code, not GIZMO source. It enrolls user boundary conditions, explicit source terms, history outputs, and mesh data for a grid calculation.
- `AGNWindFlag` is enabled in the reference, `AGNRadFlag` is enabled only with `COOLING_ENABLED == 2`, and `AGNJetFlag` is present but commented out. Jet quantities are still computed as placeholders in hot mode.
- The reference computes the inflow rate from the radial mass flux through an inner-grid face near the Bondi radius, not from a GIZMO neighbor kernel.
- Mode selection uses `mdot_crit = MdotWind_Cold(0.02 * mdot_edd) + 0.02 * mdot_edd`; above that is cold mode, with a separate super-Eddington branch above `1.66 * mdot_edd`; below it is hot mode when the inflow is positive.
- Sub-Eddington cold mode solves `MdotWind_Cold(mdot_bh) = mdot_in - mdot_bh`, where `MdotWind_Cold` scales as `0.28 * (L_bh/1e45)^0.85 Msun/yr`. Cold wind velocity scales as `2.5e4 * (L_bh/1e45)^0.4 km/s` and is capped at `0.3c`.
- Super-Eddington mode uses fitted power laws for `mdot_bh` and `v_wind`; hot mode uses a transition radius `r_tr`, a bridge polynomial when needed, `mdot_bh = mdot_bondi * sqrt(3 r_s / r_tr)` in the low-hot branch, `v_wind = 0.2 * sqrt(G M_bh / r_tr)`, and `gamma_wind = 4/3`.
- The reference queues wind/jet outflows with travel delays (`time + x1min / v`) before they reach the inner boundary.
- Boundary injection is angular: hot winds occupy 30-70 degrees from both poles, sub-Eddington cold winds default to an all-sky `cos^2(theta)` distribution, super-Eddington winds occupy polar caps within 30 degrees, and jets would occupy 0-10 degrees if enabled.
- Wind boundary gas is tagged with scalar tracers: hot wind uses `AGNWH`, cold/super-Eddington wind uses `AGNWC`, and jet uses `AGNJ`; wind metallicity comes from the `z_wind` problem parameter.
- Radiation uses a piecewise radiative efficiency from Yuan et al. 2018 eq. 25, sets `L_BH = eps_rad * mdot_bh * c^2`, propagates a radial flux, and computes AGN heating through the cooling module.

Physics context from the PDFs:

- Yuan et al. 2018 emphasizes resolving the Bondi radius, discriminating hot and cold accretion modes, and including both radiation and wind in each mode. It neglects jets for the isolated single-galaxy problem because well-collimated jets may escape with weak coupling, while noting this is a modeling assumption.
- MACER3D extends the MACER model to 3D Athena++ simulations from inside the Bondi radius to halo scales. It stresses that 3D turbulence, SN feedback, metal enrichment, and CGM treatment matter for galaxy-scale feedback. Its listed caveats include no hot-mode AGN jets yet, no cosmological context, limited CGM resolution, and missing non-thermal physics such as magnetic fields and cosmic rays.

Completed accretion pipeline:

- `blackhole_environment_loop()` gathers existing BH kernel/environment quantities.
- `blackhole_bondi_radius_loop()` computes a weighted Bondi radius for Yuan18.
- `blackhole_mass_flux_loop()` estimates inward mass flux through Fibonacci points on the Bondi sphere.
- `set_blackhole_mdot()` applies Yuan18 mode logic and stores `mdot_bh`, `mdot_wind`, `v_wind`, wind mode, injection radius, and radiative luminosity.
- `set_blackhole_new_mass()` updates the Yuan18 fall/disk reservoirs.

Key persistent Yuan18 fields are in `allvars.h` under `BH_YUAN18_ACCRETION` and `BH_YUAN18_WIND`, including `Yuan18_BH_Mass_fall`, `Yuan18_BH_Mass_disk`, `Yuan18_BH_Mdot_Bondi`, `Yuan18_BH_Bondi_Radius`, `Yuan18_BH_unspawned_wind_mass`, `Yuan18_BH_v_wind`, `Yuan18_BH_eps_wind`, `Yuan18_BH_r_inject`, `Yuan18_BH_mode_wind`, and `Yuan18_BH_J_dir`.

Current active task: wind feedback. Wind feedback is currently implemented as a spawn-based path parallel to `BH_WIND_SPAWN`, with its dispatcher and particle initialization in `blackhole_swallow_and_kick.c`. The intended relationship is that `BH_YUAN18_WIND` and `BH_WIND_SPAWN` are mutually exclusive. There is an explicit `#error` guard in `allvars.h`.

Important caution: comments in `Template_Config.sh` still describe `BH_YUAN18_WIND` as smooth/kernel-weighted wind injection, but the local implementation and `CLAUDE.md` describe the current design as spawn-based. When editing docs or flags, reconcile this wording with the actual code.

Current Yuan18 wind implementation decisions and deferred work:

- The hot-mode `r_tr` treatment in GIZMO intentionally does not impose the `x1min` upper clamp used in the grid-based `scripts/yuan18.cpp`; this was a prior modeling decision to avoid carrying over a grid-boundary artifact into the particle-based implementation. Do not "fix" this, since it is not a problem in particle-based codes.
- A physical Yuan18/MACER travel-time outflow buffer (`time + x1min / v_wind`) is still planned but intentionally deferred. The current spawn reservoir is a mass-accumulation mechanism, not the final travel-time queue.
- Splitting the wind reservoir by mode/arrival event is still planned but intentionally deferred. The current implementation can mix HOT/SUB/SUP material in one reservoir if the mode changes before a spawn event.
- Yuan18/MACER wind metallicity and tracer fields are still planned but intentionally deferred. The current implementation does not yet assign `z_wind`-based composition or hot/cold wind tracer equivalents.
- Continuous Yuan18 wind/jet feedback has an unresolved energy-partition question when multiple injection-surface samples couple to one gas cell: momentum is deposited using the net direction vector, while retaining the scalar kinetic energy of all contributing directions effectively thermalizes unresolved angular cancellation. This is likely negligible when the surface is well resolved, but should be revisited for low-resolution runs or if continuous feedback becomes the production path.

Yuan18 wind mode convention:

- `0`: none
- `1`: hot mode, biconical shell around 30-70 degrees from the spin axis
- `2`: sub-Eddington cold mode, `cos^2(theta)` angular weighting
- `3`: super-Eddington mode, polar caps within 30 degrees
- `4`: jet mode in the reference model, not yet implemented here

Radiation and jet feedback are not the current primary implementation path, but they are acceptable parallel development tracks if the user asks to split work or explore them alongside wind feedback.

## Coding Style

All project files should be written in English, including source comments, documentation, scripts, parameter notes, and diagnostics.

Follow the existing GNU-like C style and local formatting. The repository documentation recommends:

```sh
indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bap -pmt -l110 *.c
```

Use the existing macro-heavy style and keep new physics behind compile-time guards. Add prototypes to `proto.h` or module headers as appropriate. Prefer local helper patterns such as `DMAX`, `DMIN`, `BPP(i)`, `P[i]`, `SphP[i]`, `All.*`, and existing MPI/reduction utilities instead of introducing new conventions.

Global variables and struct fields usually use uppercase-leading CamelCase or existing GADGET/GIZMO names, for example `All.TimeMax`, `P[i].Mass`, and `NumForceCalculations`. Functions are generally lowercase with underscores or existing mixed local style. Match the surrounding file.

## Testing and Validation

There is no standalone unit-test harness. Validate with focused builds and representative runs:

- Build with the relevant `Config.sh` flags.
- Rebuild from clean after changing compile-time flags.
- Run small problems from `scripts/test_problems/`.
- For black hole/Yuan18 work, use `blackhole_run/bondi_test.params` and helper scripts such as `blackhole_run/crosscheck_accretion.py`, `blackhole_run/plot_mdot_history.py`, and `blackhole_run/visualize_BH.py` when appropriate.
- Compare logs such as `blackholes.txt`, `blackhole_details`, `energy.txt`, `balance.txt`, `timebin.txt`, and snapshots.

For physics changes, test both enabled and disabled compile-time paths when practical. Watch for debug `printf` output in the Yuan18 wind path before considering it production-ready.

## Repository Hygiene

The worktree may be dirty and contains generated build artifacts (`*.o`, `GIZMO`, `GIZMO_config.h`) and simulation outputs. Never revert or delete user changes without explicit instruction. Before changing source, inspect `git status --short --branch`.

Do not commit machine-local configuration or large generated outputs unless the user explicitly requests it. Avoid modifying `blackhole_run/output*`, snapshot files, `.DS_Store`, or the nested `XPHO_PROJECT/` content unless those are the actual task target.

## Commit and Pull Request Guidance

Use focused commits with short imperative summaries, such as `Fix Yuan18 wind reservoir accounting.` Pull requests should describe the physics or infrastructure change, list compile flags used, include build/run commands, mention changed outputs or diagnostics, and cite relevant issues or papers.

## Security and Permissions

Respect public vs. proprietary module boundaries described in `README.md`, `Template_Config.sh`, and `scripts/gizmo_documentation.md`. The public repository is GPL, but some modules and methods have collaboration/citation restrictions. If a module is labeled proprietary or permission-gated, do not assume access implies permission to use or share it.
