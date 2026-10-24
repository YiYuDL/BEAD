from pathlib import Path
import subprocess
import sys
import tempfile
from types import ModuleType
import unittest
from unittest.mock import patch

import pandas as pd

from bead.tools import tools


class FakeADMETModel:
    def __init__(self):
        self.calls = []

    def predict(self, smiles):
        self.calls.append(list(smiles))
        return pd.DataFrame(
            {
                "Clearance_Hepatocyte_AZ": [10.0 + index for index in range(len(smiles))],
                "Half_Life_Obach": [20.0 + index for index in range(len(smiles))],
                "Bioavailability_Ma": [0.6 + 0.1 * index for index in range(len(smiles))],
                "hERG": [0.1 + 0.1 * index for index in range(len(smiles))],
                "AUC0-t": [999.0] * len(smiles),
            },
            index=smiles,
        )


class ADMETPredictionTests(unittest.TestCase):
    def test_python_model_receives_path_objects(self):
        captured = {}

        class CapturingADMETModel:
            def __init__(self, **kwargs):
                captured.update(kwargs)

        fake_module = ModuleType("admet_ai")
        fake_module.ADMETModel = CapturingADMETModel
        tools._get_admet_model.cache_clear()
        try:
            with patch.dict(sys.modules, {"admet_ai": fake_module}):
                model = tools._get_admet_model()
        finally:
            tools._get_admet_model.cache_clear()

        self.assertIsInstance(model, CapturingADMETModel)
        self.assertEqual(captured["models_dir"], Path(tools.ADMET_MODELS_DIR))
        self.assertEqual(captured["drugbank_path"], Path(tools.ADMET_DRUGBANK_PATH))

    def test_python_api_selects_native_endpoints_and_batches(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame({"SMILES": ["CCO", "CCC", "CCN"]}).to_csv(input_path, index=False)
            model = FakeADMETModel()

            with patch.object(tools, "_get_admet_model", return_value=model):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    batch_size=2,
                    backend="python",
                    use_specialized_half_life=False,
                    use_priority_clearance_calibration=False,
                )

            result = pd.read_csv(output_path)
            self.assertEqual(result.columns.tolist(), tools.ADMET_OUTPUT_COLUMNS)
            self.assertEqual(result["smiles"].tolist(), ["CCO", "CCC", "CCN"])
            self.assertEqual(model.calls, [["CCO", "CCC"], ["CCN"]])
            self.assertAlmostEqual(
                result.loc[0, "Clearance_Hepatocyte_AZ_L_h_kg"], 1.5408
            )
            self.assertAlmostEqual(
                result.loc[1, "Clearance_Hepatocyte_AZ_L_h_kg"], 1.69488
            )
            self.assertAlmostEqual(
                result.loc[2, "Clearance_Hepatocyte_AZ_L_h_kg"], 1.5408
            )
            self.assertEqual(result["Bioavailability_Ma"].tolist(), [0.6, 0.7, 0.6])
            self.assertEqual(result["hERG"].tolist(), [0.1, 0.2, 0.1])
            self.assertNotIn("AUC0-t", result.columns)
            self.assertNotIn("CLz", result.columns)
            self.assertNotIn("F(%)", result.columns)

    def test_invalid_and_empty_smiles_keep_their_input_rows(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame(
                {"canonical_smiles": ["CCO", "not-a-smiles", None, "CCC"]}
            ).to_csv(input_path, index=False)
            model = FakeADMETModel()

            with patch.object(tools, "_get_admet_model", return_value=model):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    backend="python",
                    use_specialized_half_life=False,
                    use_priority_clearance_calibration=False,
                )

            result = pd.read_csv(output_path)
            self.assertEqual(len(result), 4)
            self.assertEqual(model.calls, [["CCO", "CCC"]])
            self.assertTrue(result.loc[1, tools.ADMET_PROPERTY_COLUMNS].isna().all())
            self.assertTrue(result.loc[2, tools.ADMET_PROPERTY_COLUMNS].isna().all())
            self.assertAlmostEqual(result.loc[3, "Half_Life_Obach"], 21.0)

    def test_official_csv_cli_backend(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame({"smiles": ["CCO", "CCC"]}).to_csv(input_path, index=False)

            def fake_run(command, **kwargs):
                cli_input = pd.read_csv(command[command.index("--data_path") + 1])
                cli_output = command[command.index("--save_path") + 1]
                pd.DataFrame(
                    {
                        "smiles": cli_input["smiles"],
                        "Clearance_Hepatocyte_AZ": [12.0, 13.0],
                        "Half_Life_Obach": [4.0, 5.0],
                        "Bioavailability_Ma": [0.8, 0.4],
                        "hERG": [0.2, 0.9],
                    }
                ).to_csv(cli_output, index=False)
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

            with (
                patch.object(tools.shutil, "which", return_value="/opt/admet/bin/admet_predict"),
                patch.object(tools.subprocess, "run", side_effect=fake_run) as run_mock,
            ):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    backend="cli",
                    use_specialized_half_life=False,
                    use_priority_clearance_calibration=False,
                )

            result = pd.read_csv(output_path)
            self.assertEqual(result.columns.tolist(), tools.ADMET_OUTPUT_COLUMNS)
            self.assertAlmostEqual(
                result.loc[0, "Clearance_Hepatocyte_AZ_L_h_kg"], 1.84896
            )
            self.assertAlmostEqual(
                result.loc[1, "Clearance_Hepatocyte_AZ_L_h_kg"], 2.00304
            )
            self.assertEqual(result["hERG"].tolist(), [0.2, 0.9])
            command = run_mock.call_args.args[0]
            self.assertEqual(command[0], "/opt/admet/bin/admet_predict")
            self.assertIn("--smiles_column", command)
            self.assertEqual(
                command[command.index("--models_dir") + 1],
                tools.ADMET_MODELS_DIR,
            )
            self.assertEqual(
                command[command.index("--drugbank_path") + 1],
                tools.ADMET_DRUGBANK_PATH,
            )

    def test_rejects_missing_endpoint_and_invalid_probability(self):
        missing = pd.DataFrame(
            {
                "Clearance_Hepatocyte_AZ": [1.0],
                "Half_Life_Obach": [2.0],
                "Bioavailability_Ma": [0.5],
            }
        )
        with self.assertRaisesRegex(ValueError, "hERG"):
            tools._select_admet_ai_columns(missing, 1)

        invalid_probability = missing.assign(hERG=[1.2])
        with self.assertRaisesRegex(ValueError, "outside"):
            tools._select_admet_ai_columns(invalid_probability, 1)

    def test_custom_clearance_scaling_and_invalid_scaling_values(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame({"smiles": ["CCO"]}).to_csv(input_path, index=False)
            model = FakeADMETModel()

            with patch.object(tools, "_get_admet_model", return_value=model):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    backend="python",
                    use_specialized_half_life=False,
                    use_priority_clearance_calibration=False,
                    hepatocellularity_million_cells_per_g=100,
                    liver_weight_g_per_kg=20,
                    fraction_unbound_in_incubation=0.5,
                )

            result = pd.read_csv(output_path)
            self.assertAlmostEqual(
                result.loc[0, "Clearance_Hepatocyte_AZ_L_h_kg"], 2.4
            )

            with self.assertRaisesRegex(
                ValueError, "fraction_unbound_in_incubation"
            ):
                tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "invalid.csv",
                    backend="python",
                    use_specialized_half_life=False,
                    use_priority_clearance_calibration=False,
                    fraction_unbound_in_incubation=0,
                )

    def test_specialized_half_life_overrides_only_half_life(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame({"smiles": ["CCO", "CCC"]}).to_csv(input_path, index=False)
            model = FakeADMETModel()

            with (
                patch.object(tools, "_get_admet_model", return_value=model),
                patch.object(
                    tools,
                    "_predict_specialized_half_life",
                    return_value=[1.25, 2.5],
                ) as half_life_mock,
            ):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    backend="python",
                    use_priority_clearance_calibration=False,
                )

            result = pd.read_csv(output_path)
            self.assertEqual(result["Half_Life_Obach"].tolist(), [1.25, 2.5])
            self.assertEqual(result["Bioavailability_Ma"].tolist(), [0.6, 0.7])
            half_life_mock.assert_called_once_with(["CCO", "CCC"])

    def test_priority_clearance_calibration_precedes_unit_scaling(self):
        with tempfile.TemporaryDirectory(prefix="bead-admet-test-") as directory:
            tmp_path = Path(directory)
            input_path = tmp_path / "input.csv"
            pd.DataFrame({"smiles": ["CCO", "CCC"]}).to_csv(input_path, index=False)
            model = FakeADMETModel()

            with (
                patch.object(tools, "_get_admet_model", return_value=model),
                patch.object(
                    tools,
                    "_calibrate_priority_clearance",
                    return_value=[20.0, 10.0],
                ) as clearance_mock,
            ):
                output_path = tools.ADMET_prediction(
                    input_path,
                    output_csv=tmp_path / "output.csv",
                    backend="python",
                    use_specialized_half_life=False,
                )

            result = pd.read_csv(output_path)
            self.assertAlmostEqual(
                result.loc[0, "Clearance_Hepatocyte_AZ_L_h_kg"], 3.0816
            )
            self.assertAlmostEqual(
                result.loc[1, "Clearance_Hepatocyte_AZ_L_h_kg"], 1.5408
            )
            self.assertEqual(clearance_mock.call_args.args[0], ["CCO", "CCC"])


if __name__ == "__main__":
    unittest.main()
