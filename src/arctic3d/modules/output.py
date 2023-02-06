"Output library."
import logging
import os
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import plotly.graph_objects as go

log = logging.getLogger("arctic3d.log")


def setup_output_folder(uniprot_id, input_files, output_dir):
    """Sets up output folder.

    Parameters
    ----------
    uniprot_id : string or None
        uniprot_id of the run
    input_files : dict of Paths
        dict of input files
    output_dir : str or None
        user-defined name of the run

    Returns
    -------
    copied_input_files : dict of Paths
        dict of copied input files
    """
    run_dir = output_dir
    if run_dir is None:
        if uniprot_id is None:
            key = "custom"
        else:
            key = uniprot_id
        run_dir = f"arctic3d-{key}"
    run_dir = Path(run_dir)

    if os.path.exists(run_dir):
        raise Exception(f"{run_dir} already exists!")

    # setting up the directory
    log.info(f"Setting up output_directory {run_dir}")
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


def write_dict(input_dict, out_filename, keyword, sep=" "):
    """
    Writes dictionary to file.

    Parameters
    ----------
    input_dict : dict
        input dictionary
    out_filename : str or Path
        name of the output filename
    keyword : str
        keyword to be used before each entry
    """
    log.info(f"Writing {keyword} information to file {out_filename}")
    with open(out_filename, "w") as wfile:
        for key in input_dict.keys():
            cl_string = sep.join([str(el) for el in input_dict[key]])
            wfile.write(f"{keyword} {key} -> " + cl_string + os.linesep)


def parse_clusters(cl_filename):
    """
    Reads clusters file.

    Parameters
    ----------
    cl_filename : str or Path
        name of the input filename

    Returns
    cl_dict : dict
        dictionary of clustered interfaces
    """
    cl_dict = {}
    with open(cl_filename, "r") as rfile:
        for ln in rfile:
            splt_ln = ln.split()
            cl_dict[splt_ln[1]] = splt_ln[3:]
    return cl_dict


def write_residues_probs(cl_residues_probs, res_probs_filename):
    """
    Writes clustered residues to file with their probability.

    Parameters
    ----------
    res_dict : dict
        dictionary of clustered residues
    res_probs_filename : str or Path
        output filename
    """
    # write to file
    with open(res_probs_filename, "w") as wfile:
        for key in cl_residues_probs.keys():
            cl_string = (
                f"Cluster {key} :"
                f" {len(cl_residues_probs[key].keys())} residues{os.linesep}"
            )
            cl_string += f"rank\tresid\tprobability{os.linesep}"
            sorted_probs = sorted(
                cl_residues_probs[key].items(),
                key=lambda x: x[1],
                reverse=True,
            )
            for pair_idx, pair in enumerate(sorted_probs, start=1):
                cl_string += (
                    f"{pair_idx}\t{pair[0]}\t{pair[1]:.3f}{os.linesep}"
                )
            cl_string += os.linesep
            wfile.write(cl_string)


def output_pdb(pdb_f, cl_residues_probs):
    """
    Outputs pdb containing probabilities.

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
                            st_beta
                            + (100 - st_beta) * cl_residues_probs[cl_id][resid]
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
    # create a qualitatively diverse set of colors
    ref_colors = plt.cm.tab20.colors
    odd_colors = [
        f"rgb{ref_colors[n]}" for n in range(len(ref_colors)) if n % 2 == 1
    ]
    even_colors = [
        f"rgb{ref_colors[n]}" for n in range(len(ref_colors)) if n % 2 == 0
    ]
    colors = odd_colors + even_colors

    # create figure
    fig = go.Figure(layout={"width": len(conv_resids) * 10, "height": 500})
    for key_idx, key in enumerate(probs.keys()):
        fig.add_trace(
            go.Bar(
                x=conv_resids,
                y=probs[key],
                name=key,
                marker_color=colors[key_idx],
            )
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
    html_output_filename = "sequence_probability.html"
    json_output_filename = "sequence_probability.json"
    fig.write_html(html_output_filename)
    fig.write_json(json_output_filename)
    log.info(f"Interactive plot {html_output_filename} successfully created.")


def plot_interactive_probs(pdb_f, cl_residues_probs):
    """
    Interactive plot.

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
    except Exception as e:
        log.warning(
            "Could not create interactive plot. The following error"
            f" occurred {e}"
        )


def make_output(
    interface_residues, pdb_f, cl_dict, cl_residues, cl_residues_probs
):
    """
    wrapper to call the different output functions.

    Parameters
    ----------
    interface_residues : dict
        dictionary of all interfaces
    pdb_f : str or Path
        Path to PDB file
    cl_dict : dict
        dictionary of clustered interfaces
    cl_residues : dict
        dictionary of clustered residues
    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        example { 1 : {1:0.7, 2:0.2, 3:0.4 ...}
                ...
        }
    """
    # writing full set of retrieved interfaces to file
    int_filename = "retrieved_interfaces.out"

    write_dict(interface_residues, int_filename, keyword="Interface")

    # writing cluster information to files
    cl_filename = "clustered_interfaces.out"
    write_dict(cl_dict, cl_filename, keyword="Cluster")

    res_filename = "clustered_residues.out"
    write_dict(cl_residues, res_filename, keyword="Cluster")

    res_probs_filename = "clustered_residues_probs.out"
    write_residues_probs(cl_residues_probs, res_probs_filename)

    # write output pdb with probabilities
    output_pdb(pdb_f, cl_residues_probs)

    # make interactive plot with probabilities
    plot_interactive_probs(pdb_f, cl_residues_probs)
