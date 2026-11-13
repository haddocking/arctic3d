"""HTML report of an ARCTIC-3D run, with a Mol* view (MolViewSpec).

The report is built in the run directory, from the clustering results and
the files written by make_output. It refers to them by relative names, so
the directory can be served as is, e.g. with
`python -m http.server -d <run_dir>`.
"""

import logging
from collections import defaultdict
from html import escape
from pathlib import Path
from string import Template

from matplotlib.colors import to_hex
from molviewspec import (
    MVSJ,
    ComponentExpression,
    GlobalMetadata,
    States,
    create_builder,
)

from arctic3d.modules.output import cluster_color

log = logging.getLogger("arctic3d.log")

MOLSTAR_VERSION = "5.13.1"
# Number of residues and partners listed per cluster.
TOP = 10
MAX_PARTNERS = 20
# Color of the probabilities in the tables (red in the structure views).
GREEN = (35, 139, 69)

TEMPLATE = Template(
    """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ARCTIC-3D Report</title>
<script src="https://cdn.jsdelivr.net/npm/molstar@${molstar}/build/mvs-stories/mvs-stories.js"></script>
<link rel="stylesheet" type="text/css"
  href="https://cdn.jsdelivr.net/npm/molstar@${molstar}/build/mvs-stories/mvs-stories.css">
<style>
  body { margin: 0; font-family: sans-serif; }
  header, section { padding: 0 1rem; }
  .container { display: flex; height: 600px; }
  #viewer { flex: 2 1 0; min-width: 0; position: relative; }
  #controls { flex: 1 1 0; min-width: 300px; padding: 1rem; overflow: auto;
              border-left: 1px solid #ccc; background: #f6f5f3; }
  table { border-collapse: collapse; margin: 0.5rem 0 1rem; }
  th, td { border: 1px solid #a0a0a0; padding: 4px 8px; text-align: right; }
  td.text { text-align: left; }
  table.sortable th { cursor: pointer; }
  .swatch { display: inline-block; width: 0.8em; height: 0.8em;
            margin-right: 0.3em; }
  img { max-width: min(100%, 800px); }
  iframe { width: 100%; height: 540px; border: none; }
  #file-hint { background: #fff3cd; padding: 0.5rem; }
  @media (max-width: 900px) {
    .container { flex-direction: column; height: auto; }
    #viewer { height: 50vh; }
    #controls { border-left: none; border-top: 1px solid #ccc; }
  }
</style>
</head>
<body>
<header>
<h2>ARCTIC-3D Report</h2>
${summary}
<details><summary>Options</summary><ul>${options}</ul></details>
<p>Download: ${downloads}</p>
<p id="file-hint" hidden>The browser does not let a page opened from disk load
the files of the viewer. Serve this directory with
<code>python -m http.server -d &lt;directory&gt;</code> and open
<a href="http://localhost:8000/report.html">http://localhost:8000/report.html</a>.</p>
</header>
<div class="container">
  <div id="viewer"><mvs-stories-viewer /></div>
  <div id="controls"><mvs-stories-snapshot-markdown /></div>
</div>
<section>
<p><i>Click a column header to sort a table by that column; click again to
reverse the order.</i></p>
${body}</section>
<script>
// Sort a table by a column on a click on its header, again to reverse.
document.querySelectorAll('table.sortable th').forEach(th => {
  th.title = 'Click to sort';
  th.addEventListener('click', () => {
    const table = th.closest('table');
    const col = th.cellIndex;
    const rows = [...table.rows].slice(1);
    const key = r => {
      const text = r.cells[col].textContent;
      const x = parseFloat(text);
      return isNaN(x) ? text : x;
    };
    const dir = th.dataset.dir === 'desc' ? 1 : -1;
    th.dataset.dir = dir === 1 ? 'asc' : 'desc';
    // Show the sorted column and its direction.
    for (const h of table.rows[0].cells) {
      h.textContent = h.textContent.replace(/ [▲▼]$$/, '');
    }
    th.textContent += dir === 1 ? ' ▲' : ' ▼';
    rows.sort((a, b) => {
      const x = key(a), y = key(b);
      return dir * (x < y ? -1 : x > y ? 1 : 0);
    });
    rows.forEach(r => table.tBodies[0].appendChild(r));
  });
});
if (location.protocol === 'file:') {
  document.getElementById('file-hint').hidden = false;
} else {
  mvsStories.loadFromURL('state.mvsj', { format: 'mvsj' });
}
</script>
</body>
</html>
"""  # noqa: E501
)


def _shade(t, rgb=(178, 24, 43)):
    """Return the hex color at fraction t from white to rgb (red)."""
    return "#{:02x}{:02x}{:02x}".format(
        *[round(255 + t * (c - 255)) for c in rgb]
    )


def _table(title, header, rows, shade_from=None):
    """Return a sortable HTML table.

    Cells are escaped, except those given as a 1-tuple (html,). With
    shade_from, the cells from that column on are probabilities, colored
    from white (0) to green (1).
    """
    html = [
        f"<h3>{escape(title)}</h3>",
        '<table class="sortable"><tr>',
    ]
    html += [f"<th>{escape(h)}</th>" for h in header]
    html.append("</tr>")
    for row in rows:
        html.append("<tr>")
        for j, c in enumerate(row):
            if shade_from is not None and j >= shade_from:
                # white text on the darkest cells
                style = f"background: {_shade(c, GREEN)}"
                if c > 0.6:
                    style += "; color: white"
                html.append(f'<td style="{style}">{c:.3f}</td>')
            elif isinstance(c, tuple):
                html.append(f'<td class="text">{c[0]}</td>')
            elif isinstance(c, str):
                html.append(f'<td class="text">{escape(c)}</td>')
            else:
                html.append(f"<td>{c}</td>")
        html.append("</tr>")
    html.append("</table>")
    return "".join(html)


def _top_residues(probs, resnames):
    """Return the TOP residues of a cluster, highest probability first."""
    ranked = sorted(probs.items(), key=lambda p: -p[1])[:TOP]
    return ", ".join(f"{resnames[r]}{r} ({p:.2f})" for r, p in ranked)


def _view(pdb_name, colors):
    """Return a builder with a cartoon and a transparent surface of the PDB.

    Both are white, with the residues colored as in colors ({resid: color}),
    also shown as ball-and-stick. Also returns the component of these
    residues.
    """
    builder = create_builder()
    structure = (
        builder.download(url=pdb_name).parse(format="pdb").model_structure()
    )
    whole = structure.component()
    residues = structure.component(
        selector=[ComponentExpression(auth_seq_id=r) for r in colors]
    )
    # one color node per color, not per residue: a smaller state.mvsj
    groups = defaultdict(list)
    for resid, color in colors.items():
        groups[color].append(ComponentExpression(auth_seq_id=resid))
    for representation in (
        whole.representation(type="cartoon"),
        whole.representation(type="surface").opacity(opacity=0.5),
        residues.representation(type="ball_and_stick"),
    ):
        representation.color(color="white")
        for color, selector in groups.items():
            representation.color(color=color, selector=selector)
    return builder, residues


def _overview(pdb_name, cl_residues_probs):
    """Return the snapshot of all clusters, in the sequence plot colors.

    A residue in several clusters takes the color of the cluster where its
    probability is highest.
    """
    best = {}
    for idx, cl_id in enumerate(cl_residues_probs):
        for resid, prob in cl_residues_probs[cl_id].items():
            if resid not in best or prob > best[resid][0]:
                best[resid] = (prob, to_hex(cluster_color(idx)))
    builder, _ = _view(
        pdb_name, {resid: color for resid, (prob, color) in best.items()}
    )
    return builder.get_snapshot(
        key="all",
        title="All binding surfaces",
        description=(
            "Residues of each cluster (binding surface) colored as in the "
            "sequence probability plot. A residue in several clusters takes "
            "the color of the cluster where its probability is highest."
        ),
        description_format="markdown",
    )


def _cluster_snapshot(pdb_name, cl_id, interfaces, probs, resnames):
    """Return the snapshot of a cluster, residues shaded by probability."""
    builder, residues = _view(
        pdb_name, {resid: _shade(prob) for resid, prob in probs.items()}
    )
    residues.focus(radius_factor=2)

    partners = ", ".join(interfaces[:MAX_PARTNERS])
    if len(interfaces) > MAX_PARTNERS:
        partners += f" and {len(interfaces) - MAX_PARTNERS} more"
    description = (
        f"**Cluster {cl_id}**: {len(interfaces)} interfaces, "
        f"{len(probs)} residues, colored by probability (fraction of the "
        "interfaces of the cluster that contain the residue) from white "
        "(0) to red (1).\n\nTop residues: "
        f"{_top_residues(probs, resnames)}.\n\nInterfaces: {partners}."
    )
    return builder.get_snapshot(
        key=f"cl{cl_id}",
        title=f"Cluster {cl_id}",
        description=description,
        description_format="markdown",
    )


def generate_report(
    pdb_f,
    interface_residues,
    cl_dict,
    cl_residues_probs,
    resnames,
    options,
    uniprot_id=None,
    localised=None,
):
    """Write report.html and state.mvsj in the run directory (cwd).

    Parameters
    ----------
    pdb_f : Path
        PDB file used for the clustering, relative to the run directory
    interface_residues : dict
        dictionary of all retrieved interfaces
    cl_dict : dict
        dictionary of clustered interfaces
    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
    resnames : dict
        dictionary of residues names (3 letter) for each residue id
    options : dict
        command-line options shown in the report
    uniprot_id : str or None
        uniprot id of the target
    localised : dict or None
        run folder of arctic3d-localise for each biological clustering type
    """
    log.info("Writing report.html")
    pdb_f = Path(pdb_f)
    pdb_name = pdb_f.as_posix()
    # numpy integers are not JSON serializable
    cl_residues_probs = {
        cl_id: {int(r): float(p) for r, p in probs.items()}
        for cl_id, probs in cl_residues_probs.items()
    }

    snapshots = [_overview(pdb_name, cl_residues_probs)]
    for cl_id, probs in cl_residues_probs.items():
        snapshots.append(
            _cluster_snapshot(pdb_name, cl_id, cl_dict[cl_id], probs, resnames)
        )
    state = MVSJ(
        data=States(
            snapshots=snapshots,
            metadata=GlobalMetadata(description="ARCTIC-3D results"),
        )
    )
    Path("state.mvsj").write_text(state.dumps(indent=None), encoding="utf-8")

    n_clustered = sum(len(ints) for ints in cl_dict.values())
    summary = (
        f"<p>Target: {escape(str(uniprot_id or 'custom'))}; structure: "
        f"{escape(pdb_name)}; {len(interface_residues)} interfaces "
        f"retrieved, {n_clustered} clustered in {len(cl_dict)} binding "
        "surfaces.</p>"
    )

    cl_ids = list(cl_residues_probs)
    tables = [
        _table(
            "Clusters (binding surfaces)",
            [
                "Cluster",
                "N interfaces",
                "N residues",
                "Top residues",
                "Interfaces",
            ],
            [
                [
                    (
                        '<span class="swatch" style="background: '
                        f'{to_hex(cluster_color(idx))}"></span>{cl_id}',
                    ),
                    len(cl_dict[cl_id]),
                    len(probs),
                    _top_residues(probs, resnames),
                    ", ".join(cl_dict[cl_id]),
                ]
                for idx, (cl_id, probs) in enumerate(cl_residues_probs.items())
            ],
        )
    ]
    resids = sorted({r for probs in cl_residues_probs.values() for r in probs})
    tables.append(
        _table(
            "Residue probabilities",
            ["Residue", "Name"] + [f"Cluster {cl_id}" for cl_id in cl_ids],
            [
                [r, resnames[r]]
                + [cl_residues_probs[cl_id].get(r, 0.0) for cl_id in cl_ids]
                for r in resids
            ],
            shade_from=2,
        )
    )
    int_cluster = {
        interface: str(cl_id)
        for cl_id, interfaces in cl_dict.items()
        for interface in interfaces
    }
    tables.append(
        _table(
            "Retrieved interfaces",
            ["Interface", "Residues", "Cluster"],
            [
                [
                    str(interface),
                    len(residues),
                    int_cluster.get(interface, "not clustered"),
                ]
                for interface, residues in interface_residues.items()
            ],
        )
    )

    figures = []
    dendrograms = sorted(str(f) for f in Path().glob("dendrogram_*.png"))
    for fname in dendrograms:
        figures.append(
            f'<h3>Dendrogram</h3><img src="{escape(fname)}" '
            'alt="Dendrogram of the interfaces">'
        )
    if Path("sequence_probability.html").exists():
        figures.append(
            "<h3>Residue probabilities along the sequence</h3>"
            '<iframe src="sequence_probability.html"></iframe>'
        )
    for clustering_type, folder in (localised or {}).items():
        figures.append(
            f"<h3>Biological clustering by {escape(clustering_type)}</h3>"
            f'<p>All files: <a href="{folder}/">{folder}/</a></p>'
        )
        # no plot for a cluster whose partners have no annotation
        for cl_id in cl_ids:
            if Path(folder, f"cluster_{cl_id}.png").exists():
                figures.append(
                    f"<h4>Cluster {cl_id} (<a href="
                    f'"{folder}/cluster_{cl_id}.html">interactive</a>)</h4>'
                    f'<img src="{folder}/cluster_{cl_id}.png" '
                    f'alt="{escape(clustering_type)} of cluster {cl_id}">'
                )

    downloads = [
        "arctic3d.log",
        pdb_name,
        "retrieved_interfaces.out",
        "interface_matrix.txt",
        "clustered_interfaces.out",
        "clustered_residues.out",
        "clustered_residues_probs.out",
        *[f"{pdb_f.stem}_cl{cl_id}{pdb_f.suffix}" for cl_id in cl_ids],
        *dendrograms,
        "sequence_probability.html",
    ]
    html = TEMPLATE.substitute(
        molstar=MOLSTAR_VERSION,
        summary=summary,
        options="".join(
            f"<li>{escape(k)}: {escape(str(v))}</li>"
            for k, v in sorted(options.items())
        ),
        downloads=" ".join(
            f'<a href="{escape(f)}">{escape(f)}</a>'
            for f in downloads
            if Path(f).exists()
        ),
        body="".join(tables + figures),
    )
    Path("report.html").write_text(html, encoding="utf-8")
