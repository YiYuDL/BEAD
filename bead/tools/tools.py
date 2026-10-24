

import math
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from PIL import Image
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import rdMolEnumerator
from rdkit.Chem import rdFMCS
from bead.utils.sim import Sim3D


ADMET_NATIVE_PROPERTY_COLUMNS = [
    "Clearance_Hepatocyte_AZ",
    "Half_Life_Obach",
    "Bioavailability_Ma",
    "hERG",
]
ADMET_CLEARANCE_NATIVE_COLUMN = "Clearance_Hepatocyte_AZ"
ADMET_CLEARANCE_OUTPUT_COLUMN = "Clearance_Hepatocyte_AZ_L_h_kg"
ADMET_PROPERTY_COLUMNS = [
    ADMET_CLEARANCE_OUTPUT_COLUMN,
    "Half_Life_Obach",
    "Bioavailability_Ma",
    "hERG",
]
ADMET_OUTPUT_COLUMNS = ["smiles", *ADMET_PROPERTY_COLUMNS]
ADMET_PROBABILITY_COLUMNS = ["Bioavailability_Ma", "hERG"]
DEFAULT_HEPATOCELLULARITY_MILLION_CELLS_PER_G = 120.0
DEFAULT_LIVER_WEIGHT_G_PER_KG = 21.4
DEFAULT_FRACTION_UNBOUND_IN_INCUBATION = 1.0
_BEAD_PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADMET_RESOURCE_DIR = os.path.join(_BEAD_PACKAGE_DIR, "resources", "admet_ai")
ADMET_MODELS_DIR = os.path.join(ADMET_RESOURCE_DIR, "models")
ADMET_DRUGBANK_PATH = os.path.join(
    ADMET_RESOURCE_DIR, "data", "drugbank_approved.csv"
)
ADMET_HALF_LIFE_MODEL_PATH = os.path.join(
    ADMET_RESOURCE_DIR,
    "half_life_model",
    "half_life_descriptor_ensemble.joblib",
)
ADMET_CLEARANCE_CALIBRATOR_PATH = os.path.join(
    ADMET_RESOURCE_DIR,
    "clearance_model",
    "clearance_priority_calibrator.joblib",
)


class ADMETAIUnavailableError(RuntimeError):
    """Raised when neither requested local ADMET-AI interface is available."""


@lru_cache(maxsize=1)
def _get_admet_model():
    """Load BEAD's repository-local ADMET-AI models once per process."""
    try:
        from admet_ai import ADMETModel
    except ImportError as exc:
        raise ADMETAIUnavailableError(
            "ADMET-AI is not installed in this Python environment. Install a "
            "compatible release (ADMET-AI v2 requires Python >= 3.11), or use "
            "backend='cli' and point BEAD_ADMET_CLI to an admet_predict executable."
        ) from exc

    models_dir = os.environ.get("BEAD_ADMET_MODELS_DIR", ADMET_MODELS_DIR)
    drugbank_path = os.environ.get("BEAD_ADMET_DRUGBANK_PATH", ADMET_DRUGBANK_PATH)
    if not os.path.isdir(models_dir):
        raise ADMETAIUnavailableError(
            f"Repository-local ADMET-AI models were not found at {models_dir!r}"
        )
    if not os.path.isfile(drugbank_path):
        raise ADMETAIUnavailableError(
            f"Repository-local DrugBank reference was not found at {drugbank_path!r}"
        )
    return ADMETModel(
        models_dir=Path(models_dir),
        drugbank_path=Path(drugbank_path),
    )


@lru_cache(maxsize=1)
def _get_specialized_half_life_model():
    """Load the locally trained, half-life-only descriptor ensemble."""
    try:
        import joblib
    except ImportError as exc:
        raise ADMETAIUnavailableError(
            "joblib is required by BEAD's specialized half-life model"
        ) from exc

    model_path = os.environ.get(
        "BEAD_ADMET_HALF_LIFE_MODEL", ADMET_HALF_LIFE_MODEL_PATH
    )
    if not os.path.isfile(model_path):
        raise ADMETAIUnavailableError(
            f"Specialized half-life model was not found at {model_path!r}"
        )
    artifact = joblib.load(model_path)
    required = {"imputer", "models", "weights", "descriptor_names"}
    missing = required - set(artifact)
    if missing:
        raise ADMETAIUnavailableError(
            "Specialized half-life model is missing: " + ", ".join(sorted(missing))
        )
    return artifact


def _predict_specialized_half_life(smiles_values):
    """Predict elimination half-life in hours with the dedicated ensemble."""
    import numpy as np
    from rdkit.Chem import Descriptors, rdFingerprintGenerator

    artifact = _get_specialized_half_life_model()
    descriptor_functions = dict(Descriptors._descList)
    try:
        selected_descriptors = [
            descriptor_functions[name] for name in artifact["descriptor_names"]
        ]
    except KeyError as exc:
        raise ADMETAIUnavailableError(
            f"Current RDKit does not provide descriptor {exc.args[0]!r} required "
            "by the specialized half-life model"
        ) from exc

    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=int(artifact.get("morgan_radius", 2)),
        fpSize=int(artifact.get("morgan_fp_size", 2048)),
    )
    feature_rows = []
    for smiles in smiles_values:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES passed to half-life model: {smiles!r}")
        fingerprint = np.asarray(
            fingerprint_generator.GetFingerprintAsNumPy(mol), dtype=np.float32
        )
        descriptor_values = np.asarray(
            [function(mol) for function in selected_descriptors], dtype=np.float64
        )
        descriptor_values[
            ~np.isfinite(descriptor_values)
            | (np.abs(descriptor_values) > np.finfo(np.float32).max)
        ] = np.nan
        feature_rows.append(
            np.concatenate((fingerprint, descriptor_values.astype(np.float32)))
        )

    features = artifact["imputer"].transform(np.asarray(feature_rows))
    model_predictions = np.column_stack(
        [model.predict(features) for model in artifact["models"].values()]
    )
    weights = np.asarray(artifact["weights"], dtype=float)
    if model_predictions.shape[1] != len(weights):
        raise ADMETAIUnavailableError(
            "Specialized half-life ensemble weights do not match its models"
        )
    predictions = model_predictions @ weights
    return np.maximum(predictions, 0.0)


@lru_cache(maxsize=1)
def _get_priority_clearance_calibrator():
    """Load the local residual calibrator for priority clearance records."""
    try:
        import joblib
    except ImportError as exc:
        raise ADMETAIUnavailableError(
            "joblib is required by BEAD's priority clearance calibrator"
        ) from exc

    model_path = os.environ.get(
        "BEAD_ADMET_CLEARANCE_CALIBRATOR", ADMET_CLEARANCE_CALIBRATOR_PATH
    )
    if not os.path.isfile(model_path):
        raise ADMETAIUnavailableError(
            f"Priority clearance calibrator was not found at {model_path!r}"
        )
    artifact = joblib.load(model_path)
    required = {
        "train_fingerprints",
        "residual_coefficients",
        "similarity_threshold",
        "kernel_power",
    }
    missing = required - set(artifact)
    if missing:
        raise ADMETAIUnavailableError(
            "Priority clearance calibrator is missing: "
            + ", ".join(sorted(missing))
        )
    return artifact


def _calibrate_priority_clearance(smiles_values, base_predictions):
    """Apply structure-local corrections to native hepatocyte clearance."""
    import numpy as np
    from rdkit.Chem import rdFingerprintGenerator

    artifact = _get_priority_clearance_calibrator()
    base = np.asarray(base_predictions, dtype=float).reshape(-1)
    if len(base) != len(smiles_values):
        raise ValueError("Clearance calibrator input lengths do not match")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=int(artifact.get("radius", 2)),
        fpSize=int(artifact.get("fp_size", 2048)),
    )
    query_fingerprints = []
    for smiles in smiles_values:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(
                f"Invalid SMILES passed to clearance calibrator: {smiles!r}"
            )
        query_fingerprints.append(
            np.asarray(generator.GetFingerprintAsNumPy(mol), dtype=np.uint8)
        )
    query = np.asarray(query_fingerprints, dtype=np.uint8)
    reference = np.asarray(artifact["train_fingerprints"], dtype=np.uint8)
    intersection = query.astype(np.float32) @ reference.astype(np.float32).T
    union = (
        query.sum(axis=1, dtype=np.float32)[:, None]
        + reference.sum(axis=1, dtype=np.float32)[None, :]
        - intersection
    )
    similarity = np.divide(
        intersection,
        union,
        out=np.zeros_like(intersection, dtype=np.float32),
        where=union > 0,
    )
    threshold = float(artifact["similarity_threshold"])
    power = float(artifact["kernel_power"])
    kernel = np.clip(
        (similarity - threshold) / (1.0 - threshold), 0.0, 1.0
    ) ** power
    corrections = kernel @ np.asarray(
        artifact["residual_coefficients"], dtype=float
    )
    return np.maximum(base + corrections, 0.0)


def _select_admet_ai_columns(predictions, expected_count):
    """Validate and select the model-native endpoints used by BEAD."""
    if isinstance(predictions, dict):
        predictions = pd.DataFrame([predictions])
    elif not isinstance(predictions, pd.DataFrame):
        try:
            predictions = pd.DataFrame(predictions)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "ADMET-AI returned an unsupported prediction result type: "
                f"{type(predictions).__name__}"
            ) from exc

    predictions = predictions.reset_index(drop=True)
    if len(predictions) != expected_count:
        raise ValueError(
            "ADMET-AI returned a different number of predictions "
            f"({len(predictions)}) than valid input molecules ({expected_count})"
        )

    missing = [
        column
        for column in ADMET_NATIVE_PROPERTY_COLUMNS
        if column not in predictions.columns
    ]
    if missing:
        raise ValueError(
            "ADMET-AI result is missing required endpoints: " + ", ".join(missing)
        )

    selected = predictions.loc[:, ADMET_NATIVE_PROPERTY_COLUMNS].copy()
    for column in ADMET_NATIVE_PROPERTY_COLUMNS:
        selected[column] = pd.to_numeric(selected[column], errors="coerce")

    for column in ADMET_PROBABILITY_COLUMNS:
        invalid_probability = selected[column].notna() & ~selected[column].between(0, 1)
        if invalid_probability.any():
            raise ValueError(f"ADMET-AI returned {column} values outside [0, 1]")
    return selected


def _predict_admet_with_python(smiles_values, batch_size):
    """Run ADMET-AI through its official Python API."""
    model = _get_admet_model()
    prediction_batches = []
    for start in range(0, len(smiles_values), batch_size):
        batch = smiles_values[start:start + batch_size]
        prediction_batches.append(
            _select_admet_ai_columns(model.predict(smiles=batch), len(batch))
        )
    if not prediction_batches:
        return pd.DataFrame(columns=ADMET_NATIVE_PROPERTY_COLUMNS)
    return pd.concat(prediction_batches, ignore_index=True)


def _resolve_admet_cli(admet_cli):
    executable = os.fspath(admet_cli)
    if os.path.dirname(executable):
        resolved = os.path.abspath(executable)
        if os.path.isfile(resolved) and os.access(resolved, os.X_OK):
            return resolved
    else:
        resolved = shutil.which(executable)
        if resolved:
            return resolved
    raise ADMETAIUnavailableError(
        f"Could not find executable {executable!r}. Install ADMET-AI or set "
        "BEAD_ADMET_CLI to the absolute path of its admet_predict command."
    )


def _predict_admet_with_cli(smiles_values, admet_cli, timeout):
    """Run ADMET-AI through its official CSV command-line interface."""
    executable = _resolve_admet_cli(admet_cli)
    models_dir = os.environ.get("BEAD_ADMET_MODELS_DIR", ADMET_MODELS_DIR)
    drugbank_path = os.environ.get("BEAD_ADMET_DRUGBANK_PATH", ADMET_DRUGBANK_PATH)
    if not os.path.isdir(models_dir):
        raise ADMETAIUnavailableError(
            f"Repository-local ADMET-AI models were not found at {models_dir!r}"
        )
    if not os.path.isfile(drugbank_path):
        raise ADMETAIUnavailableError(
            f"Repository-local DrugBank reference was not found at {drugbank_path!r}"
        )
    with tempfile.TemporaryDirectory(prefix="bead-admet-ai-") as directory:
        input_path = os.path.join(directory, "input.csv")
        output_path = os.path.join(directory, "predictions.csv")
        pd.DataFrame({"smiles": smiles_values}).to_csv(input_path, index=False)
        command = [
            executable,
            "--data_path",
            input_path,
            "--save_path",
            output_path,
            "--models_dir",
            models_dir,
            "--drugbank_path",
            drugbank_path,
            "--smiles_column",
            "smiles",
        ]
        try:
            completed = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"ADMET-AI command timed out after {timeout} seconds"
            ) from exc
        except subprocess.CalledProcessError as exc:
            details = (exc.stderr or exc.stdout or "").strip()[-1000:]
            raise RuntimeError(f"ADMET-AI command failed: {details}") from exc

        if not os.path.isfile(output_path):
            details = (completed.stderr or completed.stdout or "").strip()[-1000:]
            raise RuntimeError(
                "ADMET-AI command did not create its output CSV. " + details
            )
        return _select_admet_ai_columns(pd.read_csv(output_path), len(smiles_values))


def ADMET_prediction(
    csv_path,
    smiles_column=None,
    output_csv=None,
    batch_size=100,
    backend=None,
    admet_cli=None,
    cli_timeout=None,
    hepatocellularity_million_cells_per_g=(
        DEFAULT_HEPATOCELLULARITY_MILLION_CELLS_PER_G
    ),
    liver_weight_g_per_kg=DEFAULT_LIVER_WEIGHT_G_PER_KG,
    fraction_unbound_in_incubation=DEFAULT_FRACTION_UNBOUND_IN_INCUBATION,
    use_specialized_half_life=None,
    use_priority_clearance_calibration=None,
):
    """Predict selected endpoints locally with ADMET-AI.

    The input is a CSV containing a ``smiles`` or ``canonical_smiles`` column.
    The output contains ``smiles``, scaled hepatocyte intrinsic clearance in
    L/h/kg, half-life in hours, the binary bioavailability probability, and
    the hERG-blocker probability. By default, half-life is overridden by
    BEAD's dedicated descriptor ensemble trained in hours on the cleaned
    3,508-molecule Fan et al. dataset. Native ADMET-AI clearance is corrected
    by a structure-local residual calibrator for user-supplied priority
    measurements; below Morgan similarity 0.35, that correction is exactly
    zero. The other endpoints remain ADMET-AI predictions.

    ADMET-AI natively predicts hepatocyte clearance in uL/min/10^6 cells. It
    is converted to L/h/kg as::

        CL_scaled = CL_cell * hepatocellularity * liver_weight * 60e-6 / fu_inc

    The defaults use human scaling factors (120 million cells/g liver and
    21.4 g liver/kg) and assume ``fu_inc = 1``. The converted value is a
    scaled intrinsic clearance estimate, not observed systemic clearance.

    ``backend`` may be ``python`` (official Python API), ``cli`` (official
    ``admet_predict`` CSV command), or ``auto``. Environment variables
    ``BEAD_ADMET_BACKEND`` and ``BEAD_ADMET_CLI`` provide deployment defaults.
    ``BEAD_ADMET_MODELS_DIR`` and ``BEAD_ADMET_DRUGBANK_PATH`` may override
    the bundled resource paths when needed.
    """
    batch_size = int(batch_size)
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than zero")
    if cli_timeout is not None:
        cli_timeout = float(cli_timeout)
        if cli_timeout <= 0:
            raise ValueError("cli_timeout must be greater than zero")

    scaling_values = {
        "hepatocellularity_million_cells_per_g": hepatocellularity_million_cells_per_g,
        "liver_weight_g_per_kg": liver_weight_g_per_kg,
        "fraction_unbound_in_incubation": fraction_unbound_in_incubation,
    }
    for name, value in scaling_values.items():
        try:
            scaling_values[name] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be a positive number") from exc
        if not math.isfinite(scaling_values[name]) or scaling_values[name] <= 0:
            raise ValueError(f"{name} must be a positive number")

    backend = str(backend or os.environ.get("BEAD_ADMET_BACKEND", "auto")).lower()
    if backend not in {"auto", "python", "cli"}:
        raise ValueError("backend must be one of: auto, python, cli")
    admet_cli = admet_cli or os.environ.get("BEAD_ADMET_CLI", "admet_predict")
    if use_specialized_half_life is None:
        use_specialized_half_life = os.environ.get(
            "BEAD_ADMET_USE_SPECIALIZED_HALF_LIFE", "1"
        ).strip().lower() not in {"0", "false", "no", "off"}
    elif not isinstance(use_specialized_half_life, bool):
        raise ValueError("use_specialized_half_life must be a boolean or None")
    if use_priority_clearance_calibration is None:
        use_priority_clearance_calibration = os.environ.get(
            "BEAD_ADMET_USE_PRIORITY_CLEARANCE_CALIBRATION", "1"
        ).strip().lower() not in {"0", "false", "no", "off"}
    elif not isinstance(use_priority_clearance_calibration, bool):
        raise ValueError(
            "use_priority_clearance_calibration must be a boolean or None"
        )

    input_path = resolve_dataset_path(csv_path)
    input_df = pd.read_csv(input_path)
    if smiles_column is None:
        columns_by_lower_name = {str(column).strip().lower(): column for column in input_df.columns}
        for candidate in ("smiles", "canonical_smiles"):
            if candidate in columns_by_lower_name:
                smiles_column = columns_by_lower_name[candidate]
                break
    if smiles_column not in input_df.columns:
        raise ValueError(
            f"Could not find SMILES column {smiles_column!r}; available columns: "
            f"{list(input_df.columns)}"
        )

    smiles_values = [
        "" if pd.isna(value) else str(value).strip()
        for value in input_df[smiles_column].tolist()
    ]
    native_output_columns = ["smiles", *ADMET_NATIVE_PROPERTY_COLUMNS]
    result_rows = [
        {
            column: (smiles if column == "smiles" else pd.NA)
            for column in native_output_columns
        }
        for smiles in smiles_values
    ]
    valid_indices = [
        index
        for index, smiles in enumerate(smiles_values)
        if smiles and Chem.MolFromSmiles(smiles) is not None
    ]
    invalid_count = len(smiles_values) - len(valid_indices)
    if invalid_count:
        print(
            f"Warning: {invalid_count} empty or invalid SMILES entr"
            f"{'y was' if invalid_count == 1 else 'ies were'} left without predictions."
        )

    if valid_indices:
        valid_smiles = [smiles_values[index] for index in valid_indices]
        if backend == "python":
            predictions = _predict_admet_with_python(valid_smiles, batch_size)
        elif backend == "cli":
            predictions = _predict_admet_with_cli(valid_smiles, admet_cli, cli_timeout)
        else:
            try:
                predictions = _predict_admet_with_python(valid_smiles, batch_size)
            except ADMETAIUnavailableError as python_error:
                try:
                    predictions = _predict_admet_with_cli(valid_smiles, admet_cli, cli_timeout)
                except ADMETAIUnavailableError as cli_error:
                    raise ADMETAIUnavailableError(
                        f"Python API unavailable: {python_error} CLI unavailable: {cli_error}"
                    ) from cli_error

        if use_priority_clearance_calibration:
            predictions[ADMET_CLEARANCE_NATIVE_COLUMN] = (
                _calibrate_priority_clearance(
                    valid_smiles,
                    predictions[ADMET_CLEARANCE_NATIVE_COLUMN].to_numpy(),
                )
            )
        if use_specialized_half_life:
            predictions["Half_Life_Obach"] = _predict_specialized_half_life(
                valid_smiles
            )

        for input_index, (_, prediction) in zip(valid_indices, predictions.iterrows()):
            result_rows[input_index].update(prediction.to_dict())

    result_df = pd.DataFrame(result_rows, columns=native_output_columns)
    clearance_scale = (
        scaling_values["hepatocellularity_million_cells_per_g"]
        * scaling_values["liver_weight_g_per_kg"]
        * 60e-6
        / scaling_values["fraction_unbound_in_incubation"]
    )
    result_df[ADMET_CLEARANCE_NATIVE_COLUMN] = (
        pd.to_numeric(result_df[ADMET_CLEARANCE_NATIVE_COLUMN], errors="coerce")
        * clearance_scale
    )
    result_df = result_df.rename(
        columns={ADMET_CLEARANCE_NATIVE_COLUMN: ADMET_CLEARANCE_OUTPUT_COLUMN}
    ).loc[:, ADMET_OUTPUT_COLUMNS]
    if output_csv is None:
        save_folder_path = os.path.join(os.getcwd(), "save")
        os.makedirs(save_folder_path, exist_ok=True)
        output_path = create_unique_filename(
            os.path.join(save_folder_path, "admet_result_0.csv")
        )
    else:
        output_path = os.path.abspath(os.fspath(output_csv))
        output_parent = os.path.dirname(output_path)
        if output_parent:
            os.makedirs(output_parent, exist_ok=True)
        output_path = create_unique_filename(output_path)

    result_df.to_csv(output_path, index=False)
    missing_properties = [
        column for column in ADMET_OUTPUT_COLUMNS[1:] if result_df[column].isna().all()
    ]
    if missing_properties:
        print(
            "Warning: ADMET-AI returned no values for: "
            + ", ".join(missing_properties)
        )
    print(f"ADMET-AI results have been saved to {output_path}")
    return output_path


def _sim3d_worker(args):
    """Run one Sim3D calculation in a process-pool worker."""
    idx, query_smiles, target_smiles, num_confs, rdkit_threads = args
    try:
        if not isinstance(target_smiles, str) or not target_smiles.strip():
            return idx, 0.0, "Missing target SMILES"
        similarity = Sim3D(
            query_smiles,
            target_smiles,
            num_confs=num_confs,
            rdkit_threads=rdkit_threads,
        )
        return idx, similarity, None
    except Exception as exc:
        return idx, 0.0, str(exc)


def create_unique_filename(base_filename):
    """Return an unused path, incrementing an existing numeric suffix."""
    filename = os.fspath(base_filename)
    name, ext = os.path.splitext(filename)
    suffix_match = re.match(r"^(.*)_(\d+)$", name)
    if suffix_match:
        root = suffix_match.group(1)
        index = int(suffix_match.group(2))
    else:
        root = name
        index = 0

    while os.path.exists(filename):
        index += 1
        filename = f"{root}_{index}{ext}"
    return filename


def resolve_dataset_path(filename):

    aliases = {
        'mol_list.csv': 'mollist.csv',
        'mollist3.csv': 'mollist.csv',
        'target_molecule': 'mollist.csv',
        'target_molecules': 'mollist.csv',
        'target molecule': 'mollist.csv',
        'target molecules': 'mollist.csv',
        'dataset': 'mollist.csv',
        'database': 'mollist.csv',
        'molecule_dataset': 'mollist.csv',
        'molecule database': 'mollist.csv',
    }
    filename = os.fspath(filename)
    alias_key = filename.strip().lower()

    if os.path.exists(filename):
        return filename
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(current_dir, '..'))
    resolved_path = os.path.join(project_root, 'source', filename)
    if not os.path.exists(resolved_path) and alias_key in aliases:
        resolved_path = os.path.join(project_root, 'source', aliases[alias_key])
    
    if os.path.exists(resolved_path):
        return resolved_path
        
    raise FileNotFoundError(f"Could not find the dataset at '{filename}' or in the 'source' directory: {resolved_path}")


def align_bundle_coords(bndl):

    for m in bndl:
        Chem.SanitizeMol(m)
    mcs = rdFMCS.FindMCS(bndl,completeRingsOnly=True)
    q = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(q)
    for m in bndl:
        rdDepictor.GenerateDepictionMatching2DStructure(m,q)

def draw_mols_grid_with_highlights(
    mols,
    highlight_atom_lists=None,
    legends=None,
    sub_img_size=(300, 250),
    mols_per_row=3,
    max_mols=30,
):
    if mols_per_row <= 0:
        raise ValueError("mols_per_row must be greater than zero")
    if max_mols <= 0:
        return None


    mols = list(mols or [])[:max_mols]
    if not mols:
        return None

    if highlight_atom_lists is None:
        highlight_atom_lists = [[] for _ in mols]
    else:
        highlight_atom_lists = list(highlight_atom_lists)[:len(mols)]
        highlight_atom_lists += [[] for _ in range(len(mols) - len(highlight_atom_lists))]

    if legends is None:
        legends = ["" for _ in mols]
    else:
        legends = [str(item) for item in list(legends)[:len(mols)]]
        legends += ["" for _ in range(len(mols) - len(legends))]

    cell_w, cell_h = sub_img_size
    rows = (len(mols) + mols_per_row - 1) // mols_per_row
    grid = Image.new("RGB", (cell_w * mols_per_row, cell_h * rows), "white")

    for i, (mol, atoms, legend) in enumerate(zip(mols, highlight_atom_lists, legends)):
        if mol is None:
            continue
        valid_atoms = []
        for atom in atoms:
            try:
                atom_idx = int(atom)
            except (TypeError, ValueError):
                continue
            if 0 <= atom_idx < mol.GetNumAtoms():
                valid_atoms.append(atom_idx)
        img = Draw.MolToImage(
            mol,
            size=sub_img_size,
            highlightAtoms=valid_atoms,
            legend=legend,
        )
        if img.mode != "RGB":
            img = img.convert("RGB")
        x = (i % mols_per_row) * cell_w
        y = (i // mols_per_row) * cell_h
        grid.paste(img, (x, y))

    return grid


def display_and_save_image(img, basename):
    if img is None:
        return None
    try:
        display(img)
    except Exception as exc:
        print(f"Image display failed: {exc}")

    save_folder_path = os.path.join(os.getcwd(), "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
    image_path = create_unique_filename(os.path.join(save_folder_path, basename))
    img.save(image_path)
    print(f"Image has been saved to {image_path}")
    return image_path


def _safe_add_hs(mol, add_coords=True):
    """Add explicit hydrogens after refreshing RDKit's valence cache."""
    if mol is None:
        return None
    mol.UpdatePropertyCache(strict=False)
    return Chem.AddHs(mol, addCoords=add_coords)


def Substruc_Prepare(mol):
    '''   
    This function performs structural enumeration, adds explicit hydrogens,
    aligns coordinates, and displays images for all substructures.

    Parameters
    ------------
    mol : str
        Input molecule file path (e.g., .mol file) or SMARTS string.
    
    Returns
    -------
    qry_addHs_bundle : list
        Returns a list of hydrogenated RDKit substructure-query molecules.
    '''
    
    if not isinstance(mol, (str, os.PathLike)):
        raise TypeError("mol must be a SMILES/SMARTS string or a .mol file path")

    mol = os.fspath(mol)
    qry_addHs_bundle = []
    if mol.lower().endswith('.mol'):
        mol = resolve_dataset_path(mol)
        qry = Chem.MolFromMolFile(mol)
        if qry is None:
            print(f"Failed to load molecule file: {mol}")
            return []
    else:
        qry = Chem.MolFromSmiles(mol)
        if qry is None:
            qry = Chem.MolFromSmarts(mol)
        if qry is None:
            print(f"Failed to parse SMILES/SMARTS string: {mol}")
            return []
    
    try:
        qry_bundle = rdMolEnumerator.Enumerate(qry)
    except Exception as exc:
        raise ValueError(f"Failed to enumerate substructure {mol!r}: {exc}") from exc
    
    if len(qry_bundle) > 1:
        
        for molecule in qry_bundle:
            qry_addHs_bundle.append(_safe_add_hs(molecule, add_coords=True))
        align_bundle_coords(qry_bundle)
        print("Substructure list:")
        display(Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle)))
    
    elif len(qry_bundle) == 0:
        qry_addHs_bundle.append(_safe_add_hs(qry, add_coords=True))
        print("This substructure has only one conformation:")
        display(Draw.MolToImage(qry))
    else:
        molecule = qry_bundle[0]
        qry_addHs_bundle.append(_safe_add_hs(molecule, add_coords=True))
        print("This substructure has only one conformation:")
        display(Draw.MolToImage(molecule))
    
    return qry_addHs_bundle

def Sub(substructs, mollist):

    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    required_columns = {'canonical_smiles', 'component_synonym', 'standard_value'}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Dataset is missing required columns: {sorted(missing_columns)}")
    molecule_list = data['canonical_smiles'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    molecule_dict = {}

    for idx, molecule in zip(id_list, molecule_list):
        if not isinstance(molecule, str) or not molecule.strip():
            print(f"Skipping molecule {idx}: missing canonical_smiles")
            continue
        mol1 = Chem.MolFromSmiles(molecule)
        if mol1 is None:
            print(f"Skipping molecule {idx}: invalid SMILES {molecule!r}")
            continue
        molecule_dict[idx] = _safe_add_hs(mol1, add_coords=True)
    
    matches_no_H = []
    matched_ats = []
    legends = []
    match_info = []
    
    for one in substruct_list:
        for aa, (idx,x) in enumerate(molecule_dict.items()):
            match = x.GetSubstructMatch(one)
            if match:
                x_no_H = Chem.RemoveHs(x)
                new_match = [
                    atom.GetIdx()
                    for atom in x_no_H.GetAtoms()
                    if atom.GetIdx() in match
                ]
                matches_no_H.append(x_no_H)
                matched_ats.append(new_match)
                legends.append(f"{idx} | {data.loc[idx, 'component_synonym']} | activity_value: {data.loc[idx, 'standard_value']}")
            
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        img = draw_mols_grid_with_highlights(
            matches_no_H,
            highlight_atom_lists=matched_ats,
            sub_img_size=(300, 250),
            legends=legends,
            mols_per_row=5,
            max_mols=20,
        )
        if img is not None:
            display_and_save_image(img, "sub_matches_0.png")
    
    output_csv = 'search_result_0.csv'

    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    df = pd.DataFrame(
        match_info,
        columns=['component_synonym', 'standard_value', 'canonical_smiles', 'id'],
    )
    df.set_index('id', inplace=True)
    df.to_csv(file, index=True)
    print(f"Matching results have been saved to {file}")
    
    return file


def Subsearch(substructs, mollist):
    '''  
    Parameters
    ----------
    substructs : list/str
        Configuration SMILES/SMARTS string or file path for the substructure
    mollist : str
        Path to the CSV file containing the list of small molecules
    
    Returns
    -------
    file : str
        The path of the saved output CSV file
    '''
    
    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    required_columns = {'canonical_smiles', 'component_synonym', 'standard_value'}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Dataset is missing required columns: {sorted(missing_columns)}")
    molecule_list = data['canonical_smiles'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    molecule_dict = {}

    for idx, molecule in zip(id_list, molecule_list):
        if not isinstance(molecule, str) or not molecule.strip():
            print(f"Skipping molecule {idx}: missing canonical_smiles")
            continue
        mol1 = Chem.MolFromSmiles(molecule)
        if mol1 is None:
            print(f"Skipping molecule {idx}: invalid SMILES {molecule!r}")
            continue
        molecule_dict[idx] = _safe_add_hs(mol1, add_coords=True)
    
    matches_no_H = []
    matched_ats = []
    legends = []
    match_info = []
    
    for one in substruct_list:
        for aa, (idx,x) in enumerate(molecule_dict.items()):
            match = x.GetSubstructMatch(one)
            if match:
                x_no_H = Chem.RemoveHs(x)
                new_match = [
                    atom.GetIdx()
                    for atom in x_no_H.GetAtoms()
                    if atom.GetIdx() in match
                ]
                matches_no_H.append(x_no_H)
                matched_ats.append(new_match)
                legends.append(f"{idx} | {data.loc[idx, 'component_synonym']} | activity_value: {data.loc[idx, 'standard_value']}")
            
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        img = draw_mols_grid_with_highlights(
            matches_no_H,
            highlight_atom_lists=matched_ats,
            sub_img_size=(300, 250),
            legends=legends,
            mols_per_row=3,
            max_mols=30,
        )
        if img is not None:
            display(img)
        
    output_csv = 'search_result_0.csv'
    script_dir = os.getcwd()
   
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    df = pd.DataFrame(
        match_info,
        columns=['component_synonym', 'standard_value', 'canonical_smiles', 'id'],
    )
    df.set_index('id', inplace=True)
    df.to_csv(file, index=True)
    print(f"Matching results have been saved to {file}")
    
    return file


def Similarity_prediction(
    idx,
    csv_path=None,
    min_similarity=0.7,
    max_similarity=0.95,
    num_confs=50,
    parallel=False,
    workers=None,
    rdkit_threads=0
):
    idx = int(idx)
    min_similarity = float(min_similarity)
    max_similarity = float(max_similarity)
    if not 0.0 <= min_similarity < max_similarity <= 1.0:
        raise ValueError(
            "Similarity bounds must satisfy 0 <= min_similarity < max_similarity <= 1"
        )
    num_confs = int(num_confs)
    rdkit_threads = int(rdkit_threads)
    if num_confs < 1:
        raise ValueError("num_confs must be at least 1")
    if rdkit_threads < 0:
        raise ValueError("rdkit_threads cannot be negative")
    print(f"###The input ID is {idx}.###")
    
    if csv_path is None:
        csv_path = "mollist.csv"
    if not os.path.exists(csv_path):
        csv_path = resolve_dataset_path(csv_path)

    df = pd.read_csv(csv_path, index_col=0)
    required_columns = {'canonical_smiles', 'component_synonym', 'standard_value'}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Dataset is missing required columns: {sorted(missing_columns)}")
    df.index = df.index.astype(int)

    reference_df = df
    if idx not in reference_df.index:
        reference_path = resolve_dataset_path("mollist.csv")
        reference_df = pd.read_csv(reference_path, index_col=0)
        reference_df.index = reference_df.index.astype(int)
    if idx not in reference_df.index:
        raise ValueError(f"Reference molecule ID {idx} was not found in {csv_path} or mollist.csv.")

    print("**index1**: ", reference_df.loc[idx, "canonical_smiles"], " ", reference_df.loc[idx, "standard_value"])
    
    target_synonym = reference_df.loc[idx, 'component_synonym']
    target_df = df[df['component_synonym'] == target_synonym].copy()

    print("**target**: ", target_synonym, " candidate_count: ", len(target_df))

    if target_df.empty:
        raise ValueError(f"No molecules found associated with the target {target_synonym}.")

    query_smiles = reference_df.loc[idx, 'canonical_smiles']
    if not isinstance(query_smiles, str) or Chem.MolFromSmiles(query_smiles) is None:
        raise ValueError(f"Reference molecule ID {idx} has an invalid canonical_smiles value.")

    legends_target = []
    similarity_results = []
    legends_target_mols = []

    target_name = reference_df.loc[idx, 'component_synonym']
    activity = reference_df.loc[idx, 'standard_value']
    legends_target.append(f"{target_name} | activity_value: {activity}")

    print("Calculating 3D shape similarities. This might take a moment...")
    if parallel:
        worker_count = workers or min(os.cpu_count() or 1, len(target_df))
        worker_count = max(1, min(int(worker_count), len(target_df)))
        print(f"Using parallel Sim3D: workers={worker_count}, rdkit_threads={rdkit_threads}, num_confs={num_confs}")
        tasks = [
            (idx2, query_smiles, row['canonical_smiles'], num_confs, rdkit_threads)
            for idx2, row in target_df.iterrows()
        ]
        sim_by_idx = {}
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            for idx2, similarity, error in executor.map(_sim3d_worker, tasks):
                if error:
                    print(f"Failed to calculate 3D similarity for index {idx2}: {error}")
                sim_by_idx[idx2] = similarity
        similarity_results = [sim_by_idx[idx2] for idx2 in target_df.index]
    else:
        for idx2, row in target_df.iterrows():
            target_smiles = row['canonical_smiles']

            if isinstance(target_smiles, str) and target_smiles.strip():
                try:
                    similarity = Sim3D(
                        query_smiles,
                        target_smiles,
                        num_confs=num_confs,
                        rdkit_threads=rdkit_threads
                    )
                    similarity_results.append(similarity)
                except Exception as e:
                    print(f"Failed to calculate 3D similarity for {target_smiles}: {e}")
                    similarity_results.append(0.0)
            else:
                similarity_results.append(0.0)
    
    target_df = target_df.copy()
    target_df['similarity'] = similarity_results
    
    target_df2 = target_df[
        (target_df.index != idx)
        & (target_df['similarity'] > min_similarity)
        & (target_df['similarity'] < max_similarity)
    ].copy()
    target_df2 = target_df2.sort_values(by='similarity', ascending=False)
    
    legends_target_mols = target_df2.apply(lambda row: f"{row.name} | {target_name} | {row['standard_value']} | {row['similarity']:.3f}", axis=1).tolist()
    
    mols_target = [Chem.MolFromSmiles(reference_df.loc[idx, 'canonical_smiles'])]
    img_target = draw_mols_grid_with_highlights(
        mols_target,
        sub_img_size=(300, 250),
        legends=legends_target,
        mols_per_row=1,
        max_mols=1,
    )
    if img_target is not None:
        display(img_target)
    
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in target_df2.iterrows()]
    img_target_mols = draw_mols_grid_with_highlights(
        mols_target_mols,
        sub_img_size=(300, 250),
        legends=legends_target_mols,
        mols_per_row=3,
        max_mols=100,
    )
    if img_target_mols is not None:
        display_and_save_image(img_target_mols, "similarity_grid_0.png")

    output_csv = 'sim_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    target_df2.to_csv(file, index=True)

    print(f"Similarity results have been saved in descending order to the file: {file}")

    return file


def Filter_then_similarity(
    reference_id,
    scaffold_file='Scaffold_1.mol',
    min_similarity=0.7,
    max_similarity=0.95,
    num_confs=50,
    parallel=True,
    workers=None,
    rdkit_threads=1
):
    reference_id = int(reference_id)
    source_csv = resolve_dataset_path('mollist.csv')
    df = pd.read_csv(source_csv, index_col=0)
    df.index = df.index.astype(int)

    if reference_id not in df.index:
        raise ValueError(f"Reference molecule ID {reference_id} was not found in {source_csv}.")

    target_name = df.loc[reference_id, 'component_synonym']
    target_df = df[df['component_synonym'] == target_name].copy()

    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    target_csv = os.path.join(save_folder_path, f"{target_name}_target_0.csv")
    target_csv = create_unique_filename(target_csv)
    target_df.to_csv(target_csv, index=True)
    print(f"Target-specific molecules have been saved to {target_csv}")
    print(f"Target {target_name} candidate count before reverse subsearch: {len(target_df)}")

    reverse_csv = Rev_subsearch(scaffold_file, target_csv)
    reverse_df = pd.read_csv(reverse_csv, index_col=0)
    print(f"Candidate count after reverse subsearch: {len(reverse_df)}")

    sim_csv = Similarity_prediction(
        reference_id,
        csv_path=reverse_csv,
        min_similarity=min_similarity,
        max_similarity=max_similarity,
        num_confs=num_confs,
        parallel=parallel,
        workers=workers,
        rdkit_threads=rdkit_threads
    )
    print(f"Filter-then-similarity workflow completed. Results saved to {sim_csv}")
    return sim_csv


def Rev_subsearch(mol, csv_list):

    print(f'###{mol}###')
    print(f'####{csv_list}####')
    
    csv_list = resolve_dataset_path(csv_list)
    
    df1 =  pd.read_csv(csv_list,index_col=0)
    required_columns = {'canonical_smiles', 'component_synonym', 'standard_value'}
    missing_columns = required_columns.difference(df1.columns)
    if missing_columns:
        raise ValueError(f"Dataset is missing required columns: {sorted(missing_columns)}")

    file = Sub(mol, csv_list)
    df2 = pd.read_csv(file,index_col=0)
    print('subsearch performed')
    
    df3 = df1[~df1['canonical_smiles'].isin(df2['canonical_smiles'])]
    
    index_name = df3.index.name or 'id'
    df3_reset = df3.reset_index().rename(columns={df3.index.name or 'index': index_name})
    agg = {
        'standard_value': 'mean',
        'component_synonym': 'first',
        index_name: 'first'
    }
    if 'similarity' in df3_reset.columns:
        agg['similarity'] = 'first'
    df4 = df3_reset.groupby('canonical_smiles', as_index=False).agg(agg)
    df5 = df4.set_index(index_name)

    legends_target_mols = df5.apply(lambda row: f"{row.name} | {row['component_synonym']} | {row['standard_value']}", axis=1).tolist()
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in df5.iterrows()]
    img_target_mols = draw_mols_grid_with_highlights(
        mols_target_mols,
        sub_img_size=(300, 250),
        legends=legends_target_mols,
        mols_per_row=5,
        max_mols=100,
    )
    if img_target_mols is not None:
        display_and_save_image(img_target_mols, "reverse_result_grid_0.png")

    output_csv = 'reverse_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    df5.to_csv(file, index=True)
    
    return file
    

def merge_and_deduplicate(variables_list):
    csv_paths = []
    
    for i, var in enumerate(variables_list, 1):
        csv_path = Subsearch(var, 'mol_list.csv')
        csv_paths.append(csv_path)
        print(f"f{i}: {csv_path}")
    
    dataframes = []
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
            dataframes.append(df)
            print(f"Successfully read: {path}")
        except Exception as e:
            print(f"Failed to read file {path}: {e}")
    
    if not dataframes:
        raise ValueError("No CSV files were successfully read")
    
    merged_df = pd.concat(dataframes, ignore_index=True)
    print(f"Total rows before merging: {len(merged_df)}")
    
    deduplicated_df = merged_df.drop_duplicates(subset=['id'], keep='first')
    print(f"Rows after deduplication: {len(deduplicated_df)}")
    print(f"Number of duplicate rows deleted: {len(merged_df) - len(deduplicated_df)}")

    output_csv = 'merge_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save1")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    deduplicated_df.to_csv(file, index=False)
    
    return deduplicated_df
