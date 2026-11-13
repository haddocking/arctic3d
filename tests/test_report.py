import json
import shutil
from pathlib import Path

from arctic3d.modules.output import get_resnames_dict
from arctic3d.modules.report import generate_report

from . import golden_data


def test_generate_report(tmp_path, monkeypatch):
    """Test report.html and state.mvsj."""
    shutil.copy(Path(golden_data, "1rypB_r_b.pdb"), tmp_path)
    monkeypatch.chdir(tmp_path)
    cl_residues_probs = {2: {36: 1.0, 37: 0.5}, 1: {1: 1.0}}
    cl_dict = {2: ["int_1", "int_2"], 1: ["int_3"]}
    # int_4 was filtered out of the clustering
    interface_residues = {
        "int_1": [36],
        "int_2": [37],
        "int_3": [1],
        "int_4": [5],
    }

    generate_report(
        Path("1rypB_r_b.pdb"),
        interface_residues,
        cl_dict,
        cl_residues_probs,
        get_resnames_dict("1rypB_r_b.pdb")[1],
        {"threshold": 0.866},
    )

    with open("state.mvsj") as f:
        snapshots = json.load(f)["snapshots"]
    assert [s["metadata"]["key"] for s in snapshots] == ["all", "cl2", "cl1"]
    # Mol* loads the structure relative to the report
    for s in snapshots:
        assert s["root"]["children"][0]["params"]["url"] == "1rypB_r_b.pdb"
    html = Path("report.html").read_text()
    assert (
        '<td class="text">int_4</td><td>1</td><td class="text">not clustered'
        in html
    )
