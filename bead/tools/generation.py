"""Single-environment REINVENT and LibINVENT generation jobs for BEAD."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any


PACKAGE_DIR = Path(__file__).resolve().parents[1]
RESOURCE_DIR = PACKAGE_DIR / "resources" / "generation"
PRIORS_DIR = RESOURCE_DIR / "priors"
TEMPLATE_PATH = RESOURCE_DIR / "input-intercept.toml"
DEFAULT_RESULTS_ROOT = Path.home() / "bead_generation_results"
DEFAULT_GLIDE_GRID = Path(
    "/home/kemove/maize_generation_check/maize_demo/glide-grid_2.zip"
)
DEFAULT_SCHRODINGER = Path("/home/kemove/schrodinger2024-1")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON atomically so status readers never see partial content."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _resolve_input_file(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} was not found: {path}")
    return path


def _resolve_executable(path: Path, label: str) -> Path:
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"{label} executable was not found: {path}")
    return path


def _local_device(device: str) -> str:
    """Expose one physical GPU and return its process-local CUDA name."""
    match = re.fullmatch(r"cuda:(\d+)", str(device).strip().lower())
    if match:
        os.environ["CUDA_VISIBLE_DEVICES"] = match.group(1)
        return "cuda:0"
    if str(device).strip().lower() == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        return "cpu"
    raise ValueError("device must be 'cpu' or a physical GPU such as 'cuda:2'")


def _prepare_reinvent_config(spec: dict[str, Any], local_device: str) -> Path:
    method = spec["method"]
    job_dir = Path(spec["job_dir"])
    output_dir = job_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    prior = PRIORS_DIR / f"{method}.prior"
    _resolve_input_file(prior, f"{method} prior")

    text = TEMPLATE_PATH.read_text()
    text = re.sub(
        r"^use_cuda\s*=.*$",
        f"use_cuda = {'true' if local_device != 'cpu' else 'false'}",
        text,
        flags=re.M,
    )
    text = re.sub(r'^device\s*=.*$', f'device = "{local_device}"', text, flags=re.M)
    text = re.sub(
        r'^tb_logdir\s*=.*$',
        f'tb_logdir = "{output_dir / "tensorboard"}"',
        text,
        flags=re.M,
    )
    text = re.sub(
        r'^json_out_config\s*=.*$',
        f'json_out_config = "{output_dir / "resolved_config.json"}"',
        text,
        flags=re.M,
    )
    text = re.sub(
        r'^summary_csv_prefix\s*=.*$',
        f'summary_csv_prefix = "{output_dir / method}"',
        text,
        flags=re.M,
    )
    text = re.sub(r'^prior_file\s*=.*$', f'prior_file = "{prior}"', text, flags=re.M)
    text = re.sub(r'^agent_file\s*=.*$', f'agent_file = "{prior}"', text, flags=re.M)

    if method == "libinvent":
        scaffold = _resolve_input_file(spec["scaffold_path"], "LibINVENT scaffold file")
        text = re.sub(
            r'^\s*#?\s*smiles_file\s*=.*$',
            f'smiles_file = "{scaffold}"',
            text,
            flags=re.M,
        )
    else:
        text = re.sub(
            r'^\s*#?\s*smiles_file\s*=.*$',
            "# smiles_file is not used by REINVENT",
            text,
            flags=re.M,
        )

    max_epoch = int(spec["max_epoch"])
    batch_size = int(spec["batch_size"])
    text = re.sub(r'^batch_size\s*=.*$', f'batch_size = {batch_size}', text, flags=re.M)
    text = re.sub(r'^min_steps\s*=.*$', f'min_steps = {min(100, max_epoch)}', text, flags=re.M)
    text = re.sub(r'^max_steps\s*=.*$', f'max_steps = {max_epoch}', text, flags=re.M)
    text = re.sub(
        r"^chkpt_file\s*=.*$",
        f'chkpt_file = "{output_dir / "stage1.chkpt"}"',
        text,
        flags=re.M,
    )
    runtime_config = job_dir / f"runtime_{method}.toml"
    runtime_config.write_text(text)
    return runtime_config


def _prepare_maize_config(spec: dict[str, Any]) -> Path:
    job_dir = Path(spec["job_dir"])
    reinvent_path = Path(sys.executable).resolve().with_name("reinvent")
    reinvent = str(reinvent_path) if reinvent_path.is_file() else shutil.which("reinvent")
    if not reinvent:
        raise FileNotFoundError(
            "The reinvent executable is not installed in the active BEAD environment"
        )
    config = job_dir / "maize.toml"
    sections = []
    if spec["backend"] == "schrodinger":
        schrodinger = Path(spec["schrodinger_root"])
        glide = _resolve_executable(schrodinger / "glide", "Glide")
        ligprep = _resolve_executable(schrodinger / "ligprep", "LigPrep")
        sections.extend(
            [
                "[glide]\n" f"commands.glide = {json.dumps(str(glide))}\n",
                "[ligprep]\n" f"commands.ligprep = {json.dumps(str(ligprep))}\n",
            ]
        )
    sections.append(
        "[reinvent]\n"
        f"scripts.reinvent.interpreter = {json.dumps(sys.executable)}\n"
        f"scripts.reinvent.location = {json.dumps(reinvent)}\n"
    )
    config.write_text("\n".join(sections))
    return config


def _schrodinger_readiness(root: Path) -> tuple[bool, str]:
    """Check executables and license features before selecting the RL backend."""
    missing = [
        str(path)
        for path in (root / "glide", root / "ligprep")
        if not path.is_file() or not os.access(path, os.X_OK)
    ]
    if missing:
        return False, "missing executable(s): " + ", ".join(missing)
    lictest = root / "utilities" / "lictest"
    if not lictest.is_file() or not os.access(lictest, os.X_OK):
        return False, "Schrodinger lictest is unavailable, so licenses cannot be verified"
    environment = {
        **os.environ,
        "SCHRODINGER": str(root),
        "PATH": f"{root}:{root / 'utilities'}:{os.environ.get('PATH', '')}",
    }
    try:
        completed = subprocess.run(
            [str(lictest), "-x", "GLIDE_MAIN", "LIGPREP_MAIN"],
            capture_output=True,
            text=True,
            timeout=30,
            env=environment,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"Schrodinger license test could not run: {type(exc).__name__}: {exc}"
    if completed.returncode:
        details = (completed.stderr or completed.stdout or "license checkout failed").strip()
        return False, f"Schrodinger license test failed: {details[-500:]}"
    return True, "Glide and LigPrep executables and licenses are available"


def _check_autodock_requirements() -> None:
    for name in ("scrub.py", "mk_prepare_ligand.py"):
        executable = Path(sys.executable).resolve().with_name(name)
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(
                f"AutoDock RL requires {name} in the BEAD environment: {executable}"
            )
    try:
        import meeko  # noqa: F401
        import vina  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("AutoDock RL requires the vina and meeko Python packages") from exc


def _check_cuda(device: str) -> None:
    if device == "cpu":
        return
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError(
            f"GPU {device} was requested, but {sys.executable} uses "
            f"{torch.__version__} and CUDA is unavailable"
        )


def _patch_maize_type_matching() -> None:
    """Handle identical NumPy generic hints rejected by Beartype on Python 3.11."""
    import maize.core.graph as graph_module
    from beartype.roar import BeartypeDoorNonpepException

    original = graph_module.matching_types
    original_is_path_type = graph_module.is_path_type
    if getattr(original, "_bead_compatible", False):
        return

    def matching_types_compatible(first, second):
        try:
            return original(first, second)
        except BeartypeDoorNonpepException:
            if first == second or repr(first) == repr(second):
                return True
            raise

    matching_types_compatible._bead_compatible = True
    graph_module.matching_types = matching_types_compatible

    def is_path_type_compatible(hint):
        try:
            return original_is_path_type(hint)
        except BeartypeDoorNonpepException:
            # NumPy NDArray hints are data channels, never filesystem paths.
            if "numpy.ndarray" in repr(hint):
                return False
            raise

    graph_module.is_path_type = is_path_type_compatible


def _execute_job(spec: dict[str, Any]) -> None:
    """Build and execute the Maize reinforcement-learning feedback workflow."""
    local_device = _local_device(spec["device"])
    _check_cuda(local_device)
    maize_config = _prepare_maize_config(spec)
    runtime_config = _prepare_reinvent_config(spec, local_device)
    job_dir = Path(spec["job_dir"])
    scratch = job_dir / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)

    # Maize reads these settings while its modules are imported.
    os.environ["XDG_CONFIG_HOME"] = str(job_dir)
    os.environ["MAIZE_CONFIG"] = str(maize_config)

    from maize.core.workflow import Workflow
    from maize.steps.mai.misc import ReInvent
    import maize.steps.mai.misc.reinvent as maize_reinvent_module

    _patch_maize_type_matching()
    # maize-contrib 0.5.7 contains a developer-machine absolute log path.
    # Keep all REINVENT logs inside the isolated Maize node work directory.
    maize_reinvent_module.REINVENT_LOGFILE = Path("reinvent.log")

    flow = Workflow(
        name=f"bead-{spec['method']}-{spec['backend']}-rl",
        level="info",
        cleanup_temp=False,
    )
    flow.config.update(maize_config)
    flow.config.scratch = scratch

    generator = flow.add(ReInvent)
    if spec["backend"] == "schrodinger":
        from maize.steps.plumbing import MergeLists
        from maize.steps.mai.cheminformatics import ExtractScores, SortByTag, TagIndex
        from maize.steps.mai.docking.glide import Glide
        from maize.steps.mai.molecule.ligprep import Ligprep
        from maize.utilities.chem import IsomerCollection

        grid_path = _resolve_input_file(spec["grid_path"], "Glide grid")
        ligprep = flow.add(Ligprep, loop=True, cleanup_temp=False)
        tag_index = flow.add(TagIndex, loop=True)
        glide = flow.add(Glide, loop=True, cleanup_temp=False)
        merge = flow.add(MergeLists[IsomerCollection])
        sort_by_tag = flow.add(SortByTag, loop=True)
        scores = flow.add(ExtractScores, loop=True)

        flow.connect_all(
            (generator.out, ligprep.inp),
            (ligprep.out, tag_index.inp),
            (tag_index.out, glide.inp),
            (glide.out, merge.inp),
        )
        flow.connect_all(
            (merge.out, sort_by_tag.inp),
            (sort_by_tag.out, scores.inp),
            (scores.out, generator.inp),
        )
        glide.inp_grid.set(grid_path)
    else:
        from bead.tools.vina_rl import AutoDockVinaScore

        _check_autodock_requirements()
        receptor_path = _resolve_input_file(spec["receptor_path"], "AutoDock receptor")
        vina_score = flow.add(AutoDockVinaScore, loop=True, cleanup_temp=True)
        flow.connect_all(
            (generator.out, vina_score.inp),
            (vina_score.out, generator.inp),
        )
        vina_score.receptor.set(receptor_path)
        vina_score.search_center.set(tuple(spec["box_center"]))
        vina_score.search_range.set(tuple(spec["box_size"]))
        vina_score.scoring.set(spec["vina_scoring"])
        vina_score.exhaustiveness.set(int(spec["vina_exhaustiveness"]))
        vina_score.num_modes.set(int(spec["vina_num_modes"]))
        vina_score.seed.set(int(spec["vina_seed"]))
        vina_score.cpu.set(int(spec["vina_cpu"]))
        vina_score.ph.set(float(spec["ph"]))
        vina_score.ph_tolerance.set(float(spec["ph_tolerance"]))
        vina_score.max_states.set(int(spec["max_states"]))
        vina_score.failure_score.set(float(spec["failure_score"]))

    prior = PRIORS_DIR / f"{spec['method']}.prior"
    generator.configuration.set(runtime_config)
    generator.prior.set(prior)
    generator.agent.set(prior)
    generator.low.set(float(spec["score_low"]))
    generator.high.set(float(spec["score_high"]))
    generator.reverse.set(bool(spec["reverse_score"]))
    generator.max_epoch.set(int(spec["max_epoch"]))
    generator.batch_size.set(int(spec["batch_size"]))

    flow.check()
    if not spec.get("check_only", False):
        flow.execute()


def _collect_artifacts(job_dir: Path) -> list[str]:
    suffixes = {".chkpt", ".csv", ".json", ".log", ".toml"}
    return sorted(
        str(path)
        for path in job_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    )


def _worker(spec_path: Path) -> int:
    spec = _read_json(spec_path)
    job_dir = Path(spec["job_dir"])
    status_path = job_dir / "status.json"
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")
    _write_json(
        status_path,
        {
            "status": "running",
            "method": spec["method"],
            "requested_backend": spec["requested_backend"],
            "backend": spec["backend"],
            "backend_reason": spec["backend_reason"],
            "pid": os.getpid(),
            "python": sys.executable,
            "started_at": started_at,
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    try:
        _execute_job(spec)
    except Exception as exc:
        _write_json(
            status_path,
            {
                "status": "failed",
                "method": spec["method"],
                "requested_backend": spec["requested_backend"],
                "backend": spec["backend"],
                "backend_reason": spec["backend_reason"],
                "pid": os.getpid(),
                "python": sys.executable,
                "started_at": started_at,
                "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "job_dir": str(job_dir),
                "log": str(job_dir / "job.log"),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
                "artifacts": _collect_artifacts(job_dir),
            },
        )
        traceback.print_exc()
        return 1

    _write_json(
        status_path,
        {
            "status": "checked" if spec.get("check_only", False) else "completed",
            "method": spec["method"],
            "requested_backend": spec["requested_backend"],
            "backend": spec["backend"],
            "backend_reason": spec["backend_reason"],
            "pid": os.getpid(),
            "python": sys.executable,
            "started_at": started_at,
            "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
            "artifacts": _collect_artifacts(job_dir),
        },
    )
    return 0


def start_generation_job(
    method: str,
    grid_path: str | None = None,
    scaffold_path: str | None = None,
    max_epoch: int = 100,
    batch_size: int = 32,
    device: str | None = None,
    output_root: str | None = None,
    schrodinger_root: str | None = None,
    score_low: float = -10.0,
    score_high: float = -7.0,
    reverse_score: bool = True,
    backend: str = "auto",
    receptor_path: str | None = None,
    center_x: float | None = None,
    center_y: float | None = None,
    center_z: float | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    vina_scoring: str = "vina",
    vina_exhaustiveness: int = 8,
    vina_num_modes: int = 1,
    vina_seed: int = 20260928,
    vina_cpu: int = 1,
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    max_states: int = 1,
    failure_score: float = 0.0,
    check_only: bool = False,
) -> str:
    """Start a REINVENT/LibINVENT job with Glide or AutoDock Vina feedback."""
    method = str(method).strip().lower()
    if method not in {"reinvent", "libinvent"}:
        raise ValueError("method must be 'reinvent' or 'libinvent'")
    max_epoch = int(max_epoch)
    if max_epoch <= 0:
        raise ValueError("max_epoch must be greater than zero")
    batch_size = int(batch_size)
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than zero")
    if score_low >= score_high:
        raise ValueError("score_low must be smaller than score_high")
    requested_backend = str(backend).strip().lower()
    if requested_backend not in {"auto", "schrodinger", "autodock"}:
        raise ValueError("backend must be 'auto', 'schrodinger', or 'autodock'")
    schrodinger = Path(schrodinger_root or DEFAULT_SCHRODINGER).expanduser().resolve()
    schrodinger_ready, readiness_reason = _schrodinger_readiness(schrodinger)
    if requested_backend == "auto":
        selected_backend = "schrodinger" if schrodinger_ready else "autodock"
        backend_reason = (
            readiness_reason
            if schrodinger_ready
            else f"AutoDock fallback selected because {readiness_reason}"
        )
    elif requested_backend == "schrodinger":
        if not schrodinger_ready:
            raise RuntimeError(f"Schrodinger backend is unavailable: {readiness_reason}")
        selected_backend = "schrodinger"
        backend_reason = readiness_reason
    else:
        selected_backend = "autodock"
        backend_reason = "AutoDock backend was explicitly requested"

    grid: Path | None = None
    receptor: Path | None = None
    box_center: list[float] | None = None
    box_size = [float(size_x), float(size_y), float(size_z)]
    if selected_backend == "schrodinger":
        grid = _resolve_input_file(grid_path or DEFAULT_GLIDE_GRID, "Glide grid")
    else:
        _check_autodock_requirements()
        if not receptor_path:
            raise ValueError(
                "receptor_path is required when AutoDock is selected; use the PDBQT "
                "created by Target_preparation"
            )
        receptor = _resolve_input_file(receptor_path, "AutoDock receptor")
        if receptor.suffix.lower() != ".pdbqt":
            raise ValueError(
                "AutoDock reinforcement learning requires a prepared PDBQT receptor"
            )
        center_values = [center_x, center_y, center_z]
        if not all(value is not None for value in center_values):
            raise ValueError(
                "center_x, center_y, and center_z are required for AutoDock reinforcement learning"
            )
        box_center = [float(value) for value in center_values]
        if any(value <= 0 for value in box_size):
            raise ValueError("All AutoDock box dimensions must be greater than zero")
        vina_scoring = str(vina_scoring).strip().lower()
        if vina_scoring not in {"vina", "vinardo"}:
            raise ValueError("vina_scoring must be 'vina' or 'vinardo'")
        if min(int(vina_exhaustiveness), int(vina_num_modes), int(vina_cpu), int(max_states)) < 1:
            raise ValueError(
                "vina_exhaustiveness, vina_num_modes, vina_cpu, and max_states must be at least 1"
            )
        if not 0.0 <= float(ph) <= 14.0 or float(ph_tolerance) < 0:
            raise ValueError("ph must be in [0, 14] and ph_tolerance cannot be negative")

    scaffold: Path | None = None
    if method == "libinvent":
        if not scaffold_path:
            raise ValueError("scaffold_path is required for LibINVENT")
        scaffold = _resolve_input_file(scaffold_path, "LibINVENT scaffold file")

    selected_device = device or ("cuda:2" if method == "reinvent" else "cuda:7")
    # Validate syntax without changing the parent agent's CUDA visibility.
    if selected_device != "cpu" and not re.fullmatch(r"cuda:\d+", selected_device):
        raise ValueError("device must be 'cpu' or a physical GPU such as 'cuda:2'")

    root = Path(output_root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    date_dir = root / datetime.now().astimezone().date().isoformat()
    job_id = (
        f"{method}-{selected_backend}-{datetime.now().strftime('%H%M%S')}-"
        f"{uuid.uuid4().hex[:8]}"
    )
    job_dir = date_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=False)

    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir()
    local_grid = None
    if grid is not None:
        local_grid = inputs_dir / grid.name
        shutil.copy2(grid, local_grid)
    local_receptor = None
    if receptor is not None:
        local_receptor = inputs_dir / receptor.name
        shutil.copy2(receptor, local_receptor)
    local_scaffold = None
    if scaffold is not None:
        local_scaffold = inputs_dir / scaffold.name
        shutil.copy2(scaffold, local_scaffold)

    spec = {
        "method": method,
        "requested_backend": requested_backend,
        "backend": selected_backend,
        "backend_reason": backend_reason,
        "job_id": job_id,
        "job_dir": str(job_dir),
        "grid_path": str(local_grid) if local_grid else None,
        "receptor_path": str(local_receptor) if local_receptor else None,
        "box_center": box_center,
        "box_size": box_size,
        "scaffold_path": str(local_scaffold) if local_scaffold else None,
        "max_epoch": max_epoch,
        "batch_size": batch_size,
        "device": selected_device,
        "schrodinger_root": str(schrodinger),
        "score_low": float(score_low),
        "score_high": float(score_high),
        "reverse_score": bool(reverse_score),
        "vina_scoring": vina_scoring,
        "vina_exhaustiveness": int(vina_exhaustiveness),
        "vina_num_modes": int(vina_num_modes),
        "vina_seed": int(vina_seed),
        "vina_cpu": int(vina_cpu),
        "ph": float(ph),
        "ph_tolerance": float(ph_tolerance),
        "max_states": int(max_states),
        "failure_score": float(failure_score),
        "check_only": bool(check_only),
        "python": sys.executable,
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    spec_path = job_dir / "job.json"
    _write_json(spec_path, spec)
    status_path = job_dir / "status.json"
    _write_json(
        status_path,
        {
            "status": "queued",
            "method": method,
            "requested_backend": requested_backend,
            "backend": selected_backend,
            "backend_reason": backend_reason,
            "python": sys.executable,
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )

    log_handle = (job_dir / "job.log").open("a")
    try:
        process = subprocess.Popen(
            [sys.executable, "-m", "bead.tools.generation", "--worker", str(spec_path)],
            cwd=job_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env={
                **os.environ,
                "PYTHONUNBUFFERED": "1",
                "PATH": (
                    f"{Path(sys.executable).resolve().parent}:"
                    f"{os.environ.get('PATH', '')}"
                ),
            },
        )
    finally:
        log_handle.close()

    response = _read_json(status_path)
    response.update({"pid": process.pid, "status_file": str(status_path)})
    return json.dumps(response, indent=2, ensure_ascii=False)


def generation_job_status(job_path: str, log_lines: int = 30) -> str:
    """Return status, artifacts, and the tail of a BEAD generation job log."""
    path = Path(job_path).expanduser().resolve()
    status_path = path if path.name == "status.json" else path / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"Generation status file was not found: {status_path}")
    status = _read_json(status_path)
    log_path = Path(status.get("log", status_path.parent / "job.log"))
    if log_path.is_file() and log_lines > 0:
        status["log_tail"] = "\n".join(log_path.read_text(errors="replace").splitlines()[-log_lines:])
    status["artifacts"] = _collect_artifacts(status_path.parent)
    return json.dumps(status, indent=2, ensure_ascii=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, help="Run one serialized generation job")
    args = parser.parse_args()
    if args.worker:
        return _worker(args.worker.resolve())
    parser.error("--worker is required")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
