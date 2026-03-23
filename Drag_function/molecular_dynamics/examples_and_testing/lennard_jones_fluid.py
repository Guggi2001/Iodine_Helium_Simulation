"""
lj_lammps_ovito.py

Reference-driven LAMMPS + Python + MPI + OVITO workflow for the
"Lennard-Jones fluid" tutorial, translated into lmp.command(...) calls.

Target environment:
- Python 3.13
- Windows or Linux
- LAMMPS built in shared mode
- mpi4py installed for MPI runs
- ovito installed for post-processing

Run examples
------------
Serial:
    python lj_lammps_ovito.py

MPI:
    mpiexec -n 4 python lj_lammps_ovito.py

Notes
-----
1) This script uses the modern `lammps` Python module interface, not PyLammps.
2) On Windows, the cleanest installation is to put BOTH the Python package `lammps`
   and the LAMMPS shared library DLL into the same Python site-packages location
   using `make install-python`.
3) OVITO analysis is file-based here (`dump.lammpstrj` + thermo CSV), because that is
   the most robust workflow for PyCharm and batch validation.
4) A direct in-memory handoff path from a live LAMMPS object to OVITO is included as
   an optional helper (`ovito_live_snapshot_from_lammps`) for rank 0 only.
"""

from __future__ import annotations

import csv
import ctypes
import importlib
import math
import os
import site
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# Windows DLL handling
# -----------------------------------------------------------------------------

def _configure_windows_dll_search() -> None:
    """
    Make Windows able to locate liblammps.dll before importing the Python package.
    This is only needed for in-place or custom installations.
    If `make install-python` put the DLL into site-packages/lammps, this usually
    becomes unnecessary, but keeping it is harmless and robust.
    """
    if os.name != "nt":
        return

    candidate_dirs: list[Path] = []

    # Current script folder and common local build folders
    script_dir = Path(__file__).resolve().parent
    candidate_dirs.extend(
        [
            script_dir,
            script_dir / "lammps",
            script_dir / "bin",
            script_dir / "build",
            script_dir / "build" / "bin",
            script_dir / "build" / "Release",
            script_dir / "Release",
        ]
    )

    # User and system site-packages/lammps
    for p in site.getsitepackages():
        candidate_dirs.append(Path(p) / "lammps")
    user_site = site.getusersitepackages()
    if user_site:
        candidate_dirs.append(Path(user_site) / "lammps")

    # Optional explicit env var
    extra = os.environ.get("LAMMPS_DLL_DIR")
    if extra:
        candidate_dirs.append(Path(extra))

    seen = set()
    for folder in candidate_dirs:
        try:
            folder = folder.resolve()
        except Exception:
            continue
        if folder in seen or not folder.exists():
            continue
        seen.add(folder)
        try:
            os.add_dll_directory(str(folder))
        except (FileNotFoundError, OSError, AttributeError):
            pass


_configure_windows_dll_search()

# -----------------------------------------------------------------------------
# MPI
# -----------------------------------------------------------------------------

try:
    from mpi4py import MPI
    COMM = MPI.COMM_WORLD
    RANK = COMM.Get_rank()
    SIZE = COMM.Get_size()
except Exception:
    MPI = None
    COMM = None
    RANK = 0
    SIZE = 1

# -----------------------------------------------------------------------------
# LAMMPS import
# -----------------------------------------------------------------------------

try:
    from lammps import lammps
except Exception as exc:
    raise RuntimeError(
        "Could not import the LAMMPS Python module. "
        "Make sure LAMMPS was built in shared mode and installed via `make install-python`, "
        "or point Windows to the DLL folder via LAMMPS_DLL_DIR."
    ) from exc

# Optional OVITO imports
try:
    from ovito.io import import_file
    from ovito.modifiers import RadialDistributionFunctionModifier
    from ovito.io.lammps import lammps_to_ovito
    OVITO_AVAILABLE = True
except Exception:
    OVITO_AVAILABLE = False
    import_file = None
    RadialDistributionFunctionModifier = None
    lammps_to_ovito = None


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class LJReferenceConfig:
    # Geometry and composition from the reference tutorial
    box_lo: float = -20.0
    box_hi: float = 20.0
    n_type1: int = 1500
    n_type2: int = 100

    # Random seeds from the reference
    seed_type1: int = 34134
    seed_type2: int = 12756
    seed_langevin: int = 10917

    # LJ parameters from the reference
    mass1: float = 1.0
    mass2: float = 5.0
    epsilon11: float = 1.0
    sigma11: float = 1.0
    epsilon22: float = 0.5
    sigma22: float = 3.0
    cutoff: float = 4.0

    # Minimization
    etol: float = 1.0e-6
    ftol: float = 1.0e-6
    maxiter: int = 1000
    maxeval: int = 10000

    # MD
    timestep: float = 0.005
    nve_steps: int = 50000
    nvt_steps: int = 15000
    thermo_every_nve: int = 50
    thermo_every_nvt: int = 50
    dump_every: int = 100

    # Langevin thermostat
    temp_start: float = 1.0
    temp_stop: float = 1.0
    damp: float = 0.1

    # Analysis
    rdf_cutoff: float = 10.0
    rdf_bins: int = 200
    coordination_cutoff: float = 2.0


CFG = LJReferenceConfig()


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def mpi_print(*args, **kwargs) -> None:
    if RANK == 0:
        print(*args, **kwargs)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


# -----------------------------------------------------------------------------
# LAMMPS driver
# -----------------------------------------------------------------------------

class LJBinaryMixtureRunner:
    def __init__(self, workdir: Path, cfg: LJReferenceConfig):
        self.workdir = workdir.resolve()
        self.cfg = cfg
        ensure_dir(self.workdir)

        cmdargs = [
            "-echo", "both",
            "-log", str(self.workdir / "log.lammps"),
            "-screen", "none" if RANK != 0 else "screen",
        ]

        # Use mpi4py communicator explicitly when available.
        if COMM is not None:
            self.lmp = lammps(cmdargs=cmdargs, comm=COMM)
        else:
            self.lmp = lammps(cmdargs=cmdargs)

    def close(self) -> None:
        self.lmp.close()

    def cmd(self, line: str) -> None:
        self.lmp.command(line)

    def setup_common(self) -> None:
        c = self.cfg
        self.cmd("clear")
        self.cmd("units lj")
        self.cmd("dimension 3")
        self.cmd("atom_style atomic")
        self.cmd("boundary p p p")

        self.cmd(f"region simbox block {c.box_lo} {c.box_hi} {c.box_lo} {c.box_hi} {c.box_lo} {c.box_hi}")
        self.cmd("create_box 2 simbox")

        self.cmd(f"create_atoms 1 random {c.n_type1} {c.seed_type1} simbox overlap 0.3")
        self.cmd(f"create_atoms 2 random {c.n_type2} {c.seed_type2} simbox overlap 0.3")

        self.cmd(f"mass 1 {c.mass1}")
        self.cmd(f"mass 2 {c.mass2}")

        self.cmd(f"pair_style lj/cut {c.cutoff}")
        self.cmd(f"pair_coeff 1 1 {c.epsilon11} {c.sigma11}")
        self.cmd(f"pair_coeff 2 2 {c.epsilon22} {c.sigma22}")
        # Keep the reference mixing rule implicit:
        # epsilon_12 = sqrt(epsilon_11 * epsilon_22)
        # sigma_12   = sqrt(sigma_11   * sigma_22)

        self.cmd("neighbor 0.3 bin")
        self.cmd("neigh_modify every 1 delay 0 check yes")

        self.cmd("thermo_modify flush yes")
        self.cmd("compute mytemp all temp")
        self.cmd("compute mypress all pressure mytemp")

    def setup_thermo_csv(self, filename: str) -> None:
        """
        Record thermodynamic output to a plain text file that can be parsed reliably.
        """
        path = self.workdir / filename
        self.cmd(
            "fix fthermo all print 1 "
            "\"${step} ${temp} ${pe} ${ke} ${etotal} ${press}\" "
            f"file {path.as_posix()} screen no title \"step temp pe ke etotal press\""
        )

    def remove_thermo_csv(self) -> None:
        self.cmd("unfix fthermo")

    def part_a_minimize(self) -> None:
        c = self.cfg
        self.cmd("thermo 10")
        self.cmd("thermo_style custom step etotal press pe")
        self.cmd(f"minimize {c.etol} {c.ftol} {c.maxiter} {c.maxeval}")
        self.cmd(f"write_data {(self.workdir / 'after_min.data').as_posix()}")

    def part_b_nve(self) -> None:
        c = self.cfg
        self.cmd(f"thermo {c.thermo_every_nve}")
        self.cmd("thermo_style custom step temp etotal pe ke press")
        self.setup_thermo_csv("thermo_nve.dat")
        self.cmd("reset_timestep 0")
        self.cmd("fix mynve all nve")
        self.cmd(f"dump mydmp all atom {c.dump_every} {(self.workdir / 'dump_nve.lammpstrj').as_posix()}")
        self.cmd(f"timestep {c.timestep}")
        self.cmd(f"run {c.nve_steps}")
        self.cmd("undump mydmp")
        self.cmd("unfix mynve")
        self.remove_thermo_csv()

    def part_c_nvt_langevin(self) -> None:
        c = self.cfg
        self.cmd(f"thermo {c.thermo_every_nvt}")
        self.cmd("thermo_style custom step temp etotal pe ke press")
        self.setup_thermo_csv("thermo_nvt.dat")
        self.cmd("reset_timestep 0")
        self.cmd("fix mynve all nve")
        self.cmd(f"fix mylgv all langevin {c.temp_start} {c.temp_stop} {c.damp} {c.seed_langevin}")
        self.cmd(f"dump mydmp all atom {c.dump_every} {(self.workdir / 'dump_nvt.lammpstrj').as_posix()}")
        self.cmd(f"timestep {c.timestep}")
        self.cmd(f"run {c.nvt_steps}")
        self.cmd("undump mydmp")
        self.cmd("unfix mylgv")
        self.cmd("unfix mynve")
        self.remove_thermo_csv()

    def run_reference_workflow(self) -> None:
        self.setup_common()
        self.part_a_minimize()
        self.part_b_nve()
        self.part_c_nvt_langevin()

    def extract_counts(self) -> dict[str, int]:
        natoms = self.lmp.get_natoms()
        return {"natoms": int(natoms)}

    def version_info(self) -> dict[str, int]:
        return {
            "python_module_version": int(importlib.import_module("lammps").__version__),
            "shared_library_version": int(self.lmp.version()),
        }


# -----------------------------------------------------------------------------
# Thermo parsing
# -----------------------------------------------------------------------------

def load_thermo_table(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().strip().split()
        if header != ["step", "temp", "pe", "ke", "etotal", "press"]:
            raise ValueError(f"Unexpected thermo header in {path}: {header}")
        for line in fh:
            s = line.strip()
            if not s:
                continue
            vals = [float(x) for x in s.split()]
            if len(vals) != 6:
                continue
            rows.append(vals)
    if not rows:
        raise ValueError(f"No thermo rows parsed from {path}")
    return np.asarray(rows, dtype=float)


def plateau_mean(arr: np.ndarray, fraction_last: float = 0.2) -> float:
    start = max(0, int((1.0 - fraction_last) * len(arr)))
    return float(np.mean(arr[start:]))


def plateau_std(arr: np.ndarray, fraction_last: float = 0.2) -> float:
    start = max(0, int((1.0 - fraction_last) * len(arr)))
    return float(np.std(arr[start:], ddof=1)) if len(arr[start:]) > 1 else 0.0


# -----------------------------------------------------------------------------
# OVITO analysis
# -----------------------------------------------------------------------------

def analyze_with_ovito(workdir: Path, cfg: LJReferenceConfig) -> None:
    if not OVITO_AVAILABLE:
        mpi_print("OVITO not available; skipping OVITO analysis.")
        return

    # Analyze the NVT dump because the reference gives explicit NVT plateau values.
    dump_file = workdir / "dump_nvt.lammpstrj"
    if not dump_file.exists():
        raise FileNotFoundError(dump_file)

    pipeline = import_file(str(dump_file))
    pipeline.modifiers.append(
        RadialDistributionFunctionModifier(
            cutoff=cfg.rdf_cutoff,
            number_of_bins=cfg.rdf_bins,
            partial=True
        )
    )

    # Evaluate the last frame.
    data = pipeline.compute(pipeline.num_frames - 1)

    # RDF table
    rdf_table = data.tables["coordination-rdf"]
    rdf_xy = np.asarray(rdf_table.xy())

    rdf_csv = workdir / "rdf_lastframe.csv"
    np.savetxt(
        rdf_csv,
        rdf_xy,
        delimiter=",",
        header="r,g_r",
        comments=""
    )

    # Coordination numbers per particle if present
    coord_csv = workdir / "coordination_lastframe.csv"
    if "Coordination" in data.particles:
        coord = np.asarray(data.particles["Coordination"])
        np.savetxt(
            coord_csv,
            coord,
            delimiter=",",
            header="coordination",
            comments=""
        )

    # RDF plot
    plt.figure(figsize=(7, 4.5))
    plt.plot(rdf_xy[:, 0], rdf_xy[:, 1])
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.title("RDF from OVITO (last NVT frame)")
    plt.tight_layout()
    plt.savefig(workdir / "rdf_lastframe.png", dpi=180)
    plt.close()


def ovito_live_snapshot_from_lammps(lmp_obj) -> Optional[object]:
    """
    Optional direct handoff from an in-memory LAMMPS object to an OVITO DataCollection.
    This returns data only on rank 0 in parallel runs, per OVITO documentation.
    """
    if not OVITO_AVAILABLE:
        return None
    return lammps_to_ovito(lmp_obj)


# -----------------------------------------------------------------------------
# Validation against reference text
# -----------------------------------------------------------------------------

def validate_against_reference(workdir: Path) -> dict[str, float]:
    nve = load_thermo_table(workdir / "thermo_nve.dat")
    nvt = load_thermo_table(workdir / "thermo_nvt.dat")

    # columns: step temp pe ke etotal press
    nve_step, nve_temp, nve_pe, nve_ke, nve_etot, nve_press = nve.T
    nvt_step, nvt_temp, nvt_pe, nvt_ke, nvt_etot, nvt_press = nvt.T

    results = {
        "nve_pe_plateau_mean": plateau_mean(nve_pe),
        "nve_ke_plateau_mean": plateau_mean(nve_ke),
        "nve_etotal_plateau_mean": plateau_mean(nve_etot),
        "nve_etotal_plateau_std": plateau_std(nve_etot),
        "nvt_temp_plateau_mean": plateau_mean(nvt_temp),
        "nvt_temp_plateau_std": plateau_std(nvt_temp),
        "nvt_pe_plateau_mean": plateau_mean(nvt_pe),
        "nvt_pe_plateau_std": plateau_std(nvt_pe),
        "nvt_ke_plateau_mean": plateau_mean(nvt_ke),
        "nvt_ke_plateau_std": plateau_std(nvt_ke),
    }

    # Reference tutorial values are approximate:
    # T ≈ 1.0, PE ≈ -0.25, KE ≈ 1.5 in NVT.
    ref = {
        "nvt_temp_reference": 1.0,
        "nvt_pe_reference": -0.25,
        "nvt_ke_reference": 1.5,
    }
    results.update(ref)

    # Absolute deviations
    results["abs_err_temp"] = abs(results["nvt_temp_plateau_mean"] - ref["nvt_temp_reference"])
    results["abs_err_pe"] = abs(results["nvt_pe_plateau_mean"] - ref["nvt_pe_reference"])
    results["abs_err_ke"] = abs(results["nvt_ke_plateau_mean"] - ref["nvt_ke_reference"])

    return results


def write_validation_report(workdir: Path, results: dict[str, float]) -> None:
    out = workdir / "validation_report.txt"
    lines = []
    lines.append("Validation against tutorial reference")
    lines.append("=" * 72)
    lines.append("")
    lines.append("NVT plateau comparison")
    lines.append(f"Temperature: user={results['nvt_temp_plateau_mean']:.6f} ± {results['nvt_temp_plateau_std']:.6f} | ref≈{results['nvt_temp_reference']:.6f} | abs.err={results['abs_err_temp']:.6f}")
    lines.append(f"Potential E: user={results['nvt_pe_plateau_mean']:.6f} ± {results['nvt_pe_plateau_std']:.6f} | ref≈{results['nvt_pe_reference']:.6f} | abs.err={results['abs_err_pe']:.6f}")
    lines.append(f"Kinetic E:   user={results['nvt_ke_plateau_mean']:.6f} ± {results['nvt_ke_plateau_std']:.6f} | ref≈{results['nvt_ke_reference']:.6f} | abs.err={results['abs_err_ke']:.6f}")
    lines.append("")
    lines.append("NVE plateau diagnostics")
    lines.append(f"PE plateau mean      = {results['nve_pe_plateau_mean']:.6f}")
    lines.append(f"KE plateau mean      = {results['nve_ke_plateau_mean']:.6f}")
    lines.append(f"Etot plateau mean    = {results['nve_etotal_plateau_mean']:.6f}")
    lines.append(f"Etot plateau std     = {results['nve_etotal_plateau_std']:.6e}")
    write_text(out, "\n".join(lines))


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def make_energy_plots(workdir: Path) -> None:
    for name in ("nve", "nvt"):
        data = load_thermo_table(workdir / f"thermo_{name}.dat")
        step, temp, pe, ke, etotal, press = data.T

        plt.figure(figsize=(7, 4.5))
        plt.plot(step, pe, label="PE")
        plt.plot(step, ke, label="KE")
        plt.plot(step, etotal, label="Etotal")
        plt.xlabel("Step")
        plt.ylabel("Energy (LJ units)")
        plt.title(f"{name.upper()} energies")
        plt.legend()
        plt.tight_layout()
        plt.savefig(workdir / f"energies_{name}.png", dpi=180)
        plt.close()

        plt.figure(figsize=(7, 4.5))
        plt.plot(step, temp)
        plt.xlabel("Step")
        plt.ylabel("Temperature (LJ units)")
        plt.title(f"{name.upper()} temperature")
        plt.tight_layout()
        plt.savefig(workdir / f"temperature_{name}.png", dpi=180)
        plt.close()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    workdir = Path.cwd() / "lj_python_api_run"
    ensure_dir(workdir)

    mpi_print(f"Working directory: {workdir}")
    mpi_print(f"MPI ranks: {SIZE}")

    runner = LJBinaryMixtureRunner(workdir=workdir, cfg=CFG)
    try:
        versions = runner.version_info()
        if RANK == 0:
            print("LAMMPS module version:", versions["python_module_version"])
            print("LAMMPS library version:", versions["shared_library_version"])
            if versions["python_module_version"] != versions["shared_library_version"]:
                raise RuntimeError(
                    "LAMMPS Python module version does not match shared library version."
                )

        runner.run_reference_workflow()

        if COMM is not None:
            COMM.Barrier()

        if RANK == 0:
            make_energy_plots(workdir)
            analyze_with_ovito(workdir, CFG)
            results = validate_against_reference(workdir)
            write_validation_report(workdir, results)

            print("\nValidation summary")
            print("-" * 72)
            print(f"NVT T  plateau: {results['nvt_temp_plateau_mean']:.6f} (ref ≈ {results['nvt_temp_reference']:.6f})")
            print(f"NVT PE plateau: {results['nvt_pe_plateau_mean']:.6f} (ref ≈ {results['nvt_pe_reference']:.6f})")
            print(f"NVT KE plateau: {results['nvt_ke_plateau_mean']:.6f} (ref ≈ {results['nvt_ke_reference']:.6f})")
            print(f"NVE Etot std:   {results['nve_etotal_plateau_std']:.6e}")

    finally:
        runner.close()
        if MPI is not None:
            MPI.Finalize()


if __name__ == "__main__":
    main()