import logging
import os
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import plotly.graph_objects as go

log = logging.getLogger("arctic3dlog")


def setup_output_folder(uniprot_id, input_files):
    """Sets up output folder

    Parameters
    ----------
    uniprot_id : string or None
        uniprot_id of the run
    input_files : dict of Paths
        dict of input files
    Returns
    -------
    copied_input_files : dict of Paths
        dict of copied input files
    """
    if uniprot_id is None:
        key = "custom"
    else:
        key = uniprot_id
    run_dir = Path(f"arctic3d-{key}")
    if os.path.exists(run_dir):
        raise Exception(f"{run_dir} already exists!")
    # setting up the directory
    os.mkdir(run_dir)
    datadir = Path(run_dir, "input_data")
    os.mkdir(datadir)
    for key in input_files:
        filename = input_files[key].name
        filepath = Path(datadir, filename)
        shutil.copy(input_files[key], filepath)
    os.chdir(run_dir)
    # retrieving info about copied input files
    copied_input_files = {}
    for key in input_files:
        filename = input_files[key].name
        copied_input_files[key] = Path("input_data", filename)
    return copied_input_files


def output_pdb(pdb_f, cl_residues_probs):
    """Outputs pdb containing probabilities

    Parameters
    ----------
    pdb_f : str or Path
        Path to PDB file.

    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        example { 1 : {1:0.7, 2:0.2, 3:0.4 ...}
                  ...
                }

    Returns
    -------
    output_files : list of Paths
        List of Paths to output PDB files
    """
    output_files = []  # list of paths
    log.info("Creating output pdb files...")
    st_beta = 50.0  # baseline value for b-factors of the identified residues
    # read original file
    original_content = open(pdb_f, "r").read().split(os.linesep)
    # iterate over clusters
    for cl_id in cl_residues_probs.keys():
        # creating new file
        new_filename = f"{pdb_f.stem}_cl{cl_id}{pdb_f.suffix}"
        if os.path.exists(new_filename):
            raise Exception(f"Existing pdb file {new_filename}")
        output_files.append(Path(new_filename))
        log.info(f"Creating output file {new_filename}")
        # writing new file
        with open(new_filename, "w") as wfile:
            for ln in original_content:
                if ln.startswith("ATOM"):
                    resid = int(ln[23:27].strip())
                    new_beta = 0.00
                    if resid in cl_residues_probs[cl_id].keys():
                        new_beta = (
                            st_beta + (100 - st_beta) * cl_residues_probs[cl_id][resid]
                        )
                    n_blank = 3 - len(str(new_beta).split(".")[0])
                    new_line = f"{ln[:60]}{' '*n_blank}{new_beta:.2f}{ln[66:]}"
                    wfile.write(f"{new_line}{os.linesep}")
                else:
                    wfile.write(f"{ln}{os.linesep}")
    return output_files


def make_plotly_plot(conv_resids, probs):
    """
    Makes a plotly interactive plot.

    Parameters
    ----------
    conv_resids : list
        list of residues (x axis)
    probs : dict
        dictionary of probability (y axis)
    """
    log.info("Creating interactive plot.")
    colors = [f"rgb{key}" for key in plt.cm.Set1.colors]
    fig = go.Figure(layout={"width": len(conv_resids) * 10, "height": 500})
    for key_idx, key in enumerate(probs.keys()):
        fig.add_trace(
            go.Bar(x=conv_resids, y=probs[key], name=key, marker_color=colors[key_idx])
        )
    # prettifying layout
    fig.update_layout(
        title="ARCTIC3D clustering",
        xaxis=dict(
            title="Residue ID",
            tickfont_size=14,
            titlefont_size=16,
            tick0=conv_resids[0],
            dtick=10,
        ),
        yaxis=dict(
            title="Probability",
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
        barmode="group",
        bargap=0.05,
        bargroupgap=0.05,
        hovermode="x unified",
        hoverlabel=dict(font_size=16, font_family="Helvetica"),
    )
    # save html
    output_filename = "sequence_probability.html"
    fig.write_html(output_filename)
    log.info(f"Interactive plot {output_filename} successfully created.")


def plot_interactive_probs(pdb_f, cl_residues_probs):
    """
    Interactive plot

    Parameters
    ----------
    pdb_f : str or Path
        Path to PDB file.

    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        example { 1 : {1:0.7, 2:0.2, 3:0.4 ...}
                  ...
                }

    Returns
    -------
    html_filename : Path
        Path to output html file
    """
    mdu = mda.Universe(pdb_f)
    calphas = mdu.select_atoms("name CA")
    resids = calphas.resids
    resnames = [mda.lib.util.convert_aa_code(x) for x in calphas.resnames]
    conv_resids = [f"{resids[n]}-{resnames[n]}" for n in range(len(resids))]
    # create probs
    probs = {}
    for cl_id in cl_residues_probs.keys():
        new_probs = np.zeros(len(resids))
        for n in range(len(resids)):
            if resids[n] in cl_residues_probs[cl_id].keys():
                new_probs[n] = round(cl_residues_probs[cl_id][resids[n]], 2)
        probs[f"Cluster {cl_id}"] = new_probs
    # plotly
    try:
        make_plotly_plot(conv_resids, probs)
    except Exception:
        log.warning("Could not create interactive plot")


def make_output(pdb_f, cl_residues_probs):
    """
    wrapper to call the output functions

    Parameters
    ----------
    pdb_f : str or Path
        Path to PDB file.

    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        example { 1 : {1:0.7, 2:0.2, 3:0.4 ...}
                  ...
                }
    """
    output_pdb(pdb_f, cl_residues_probs)

    plot_interactive_probs(pdb_f, cl_residues_probs)

    # TODO : place write calls here
