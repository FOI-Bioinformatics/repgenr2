"""The combined report: one self-contained HTML file from step fragments.

Built with the standard library: every value read from data passes through
``html.escape``. Figures and table downloads are embedded as data URIs, so
the file opens without a network connection and can be archived or sent as
it is. Tables and figures are numbered across the report in step order.
"""

from __future__ import annotations

import base64
import html
import struct
from dataclasses import dataclass, field
from pathlib import Path

from .figures import PNG_DPI
from .methods import CITATIONS, parameter_rows, sentence, software_versions, unique_records
from .model import FigureRef, StepReport, Table
from .tables import read_rows

# Rows of a table shown in the report; the TSV holds all.
PREVIEW_ROWS = 25
# Embed a table for download up to this size, and an SVG figure (used when a
# figure has no PNG) up to this size.
MAX_EMBED_TABLE = 2 * 1024 * 1024
MAX_EMBED_SVG = 3 * 512 * 1024


@dataclass
class Provenance:
    """Where and when the report was made."""

    title: str
    generated: str
    repgenr_version: str
    workdir: str = ""
    command: str = ""
    # (stage, state, tool, completed)
    stages: list[tuple[str, str, str, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


_CSS = """
:root { --fg: #1a1a1a; --muted: #5f6368; --rule: #d0d4d9; --zebra: #f5f6f8; --accent: #0072B2; }
* { box-sizing: border-box; }
body { margin: 0; color: var(--fg); background: #fff;
  font: 15px/1.5 -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }
main { max-width: 62rem; margin: 0 auto; padding: 2rem 1.25rem 4rem; }
h1 { font-size: 1.7rem; margin: 0 0 .25rem; }
h2 { font-size: 1.3rem; margin: 2.5rem 0 .75rem; padding-bottom: .3rem;
  border-bottom: 1px solid var(--rule); }
h3 { font-size: 1.05rem; margin: 1.5rem 0 .5rem; }
p, li { max-width: 48rem; }
.meta { color: var(--muted); font-size: .9rem; }
.keys { display: grid; grid-template-columns: repeat(auto-fill, minmax(11rem, 1fr));
  gap: .5rem; margin: .75rem 0 1rem; padding: 0; }
.keys div { border: 1px solid var(--rule); border-radius: 4px; padding: .5rem .7rem; }
.keys dt { color: var(--muted); font-size: .8rem; }
.keys dd { margin: 0; font-size: 1.15rem; font-variant-numeric: tabular-nums; }
table { border-collapse: collapse; margin: .25rem 0 .5rem; font-size: .85rem; }
caption { caption-side: top; text-align: left; padding: 0 0 .35rem; font-size: .9rem; }
th, td { padding: .25rem .6rem; border-bottom: 1px solid var(--rule); text-align: left;
  vertical-align: top; }
th { border-bottom: 1.5px solid var(--fg); font-weight: 600; }
tbody tr:nth-child(even) { background: var(--zebra); }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
.tablewrap { overflow-x: auto; margin-bottom: 1.25rem; }
.tablenote { color: var(--muted); font-size: .8rem; margin: 0; }
figure { margin: 1rem 0 1.75rem; }
figure img { max-width: 100%; height: auto; display: block; }
figcaption { font-size: .9rem; margin-top: .4rem; max-width: 48rem; }
.label { font-weight: 600; }
.note { color: var(--muted); font-size: .9rem; }
nav ol { padding-left: 1.2rem; }
a { color: var(--accent); }
@page { size: A4; margin: 18mm; }
@media print {
  body { font-size: 10.5pt; }
  main { max-width: none; padding: 0; }
  h2 { break-before: page; }
  figure, table, .keys { break-inside: avoid; }
  nav, .download { display: none; }
}
"""


def _e(value: object) -> str:
    return html.escape(str(value), quote=True)


def _num(value: str) -> str:
    """A number in a table cell, with thousands separators; text unchanged."""
    try:
        if value.strip() and "." not in value and "e" not in value.lower():
            return f"{int(value):,}"
        f = float(value)
    except ValueError:
        return value
    return f"{f:,.2f}" if abs(f) < 1e6 else f"{f:,.0f}"


def _png_width(data: bytes) -> int | None:
    """Pixel width from a PNG header."""
    if data[:8] != b"\x89PNG\r\n\x1a\n" or len(data) < 24:
        return None
    return struct.unpack(">I", data[16:20])[0]


def _figure_img(folder: Path, fig: FigureRef) -> tuple[str, int | None] | None:
    """The data URI of a figure and its display width in CSS pixels.

    The PNG is embedded: an SVG keeps its text as text, and a browser
    without the figure font would set it in another face. The SVG (and PDF)
    stay in the step folder for manuscripts. The PNG is shown at its
    physical size (300 dpi), so a single-column figure is not stretched.
    """
    png = fig.files.get("png")
    if png and (folder / png).is_file():
        data = (folder / png).read_bytes()
        px = _png_width(data)
        width = round(px * 96 / PNG_DPI) if px else None
        return "data:image/png;base64," + base64.b64encode(data).decode(), width
    svg = fig.files.get("svg")
    if svg and (folder / svg).is_file() and (folder / svg).stat().st_size <= MAX_EMBED_SVG:
        data = (folder / svg).read_bytes()
        return "data:image/svg+xml;base64," + base64.b64encode(data).decode(), None
    return None


def _table_html(folder: Path, table: Table, number: int, step: str) -> str:
    path = folder / table.file
    parts = [f'<div class="tablewrap" id="{_e(step)}-{_e(table.id)}">']
    parts.append("<table>")
    parts.append(
        f'<caption><span class="label">Table {number}.</span> {_e(table.caption)}</caption>'
    )
    if not path.is_file():
        parts.append("</table>")
        parts.append(f'<p class="tablenote">The table file {_e(table.file)} is missing.</p></div>')
        return "".join(parts)
    header, rows = read_rows(path, PREVIEW_ROWS)
    numeric = set(table.numeric)
    parts.append("<thead><tr>")
    for col in header:
        cls = ' class="num"' if col in numeric else ""
        parts.append(f"<th{cls}>{_e(col)}</th>")
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        for col, value in zip(header, row, strict=False):
            if col in numeric:
                parts.append(f'<td class="num">{_e(_num(value))}</td>')
            else:
                parts.append(f"<td>{_e(value)}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table>")
    notes = []
    if table.n_rows > len(rows):
        notes.append(f"Showing {len(rows):,} of {table.n_rows:,} rows.")
    if table.note:
        notes.append(table.note)
    size = path.stat().st_size
    if size <= MAX_EMBED_TABLE:
        encoded = base64.b64encode(path.read_bytes()).decode()
        uri = "data:text/tab-separated-values;base64," + encoded
        notes.append(
            f'<a class="download" download="{_e(step)}_{_e(table.file)}" href="{uri}">'
            f"Download the table (TSV)</a>"
        )
    else:
        notes.append(f"The full table is in report/{_e(step)}/{_e(table.file)}.")
    parts.append(
        '<p class="tablenote">'
        + " ".join(n if n.startswith("<a ") else _e(n) for n in notes)
        + "</p></div>"
    )
    return "".join(parts)


def _figure_html(folder: Path, fig: FigureRef, number: int, step: str) -> str:
    img = _figure_img(folder, fig)
    cap = f'<figcaption><span class="label">Figure {number}.</span> {_e(fig.caption)}</figcaption>'
    if img is None:
        return f'<p class="note">Figure {number} could not be embedded: {_e(fig.caption)}</p>'
    uri, width = img
    size = f' width="{width}"' if width else ""
    return (
        f'<figure id="{_e(step)}-{_e(fig.id)}">'
        f'<img alt="{_e(fig.alt)}"{size} src="{uri}">{cap}</figure>'
    )


def build_html(fragments: list[tuple[Path, StepReport]], prov: Provenance) -> str:
    """The report page for ``fragments`` (folder, report), in the given order."""
    out: list[str] = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{_e(prov.title)}</title>",
        f"<style>{_CSS}</style></head><body><main>",
        f"<h1>{_e(prov.title)}</h1>",
        f'<p class="meta">Generated {_e(prov.generated)} by RepGenR2 {_e(prov.repgenr_version)}'
        + (f" from {_e(prov.workdir)}" if prov.workdir else "")
        + ".</p>",
    ]
    # Summary
    out.append('<h2 id="summary">Summary</h2>')
    if not fragments:
        out.append('<p class="note">No step of this run has results to report.</p>')
    for _, rep in fragments:
        if rep.summary:
            out.append(f"<p><strong>{_e(rep.title)}.</strong> {_e(' '.join(rep.summary))}</p>")
    keys = [(rep, k) for _, rep in fragments for k in rep.key_numbers]
    if keys:
        out.append('<dl class="keys">')
        for _, k in keys:
            value = k.value if not isinstance(k.value, int) else f"{k.value:,}"
            value = "not recorded" if value is None else value
            unit = f" {k.unit}" if k.unit else ""
            out.append(f"<div><dt>{_e(k.label)}</dt><dd>{_e(value)}{_e(unit)}</dd></div>")
        out.append("</dl>")
    for note in prov.notes:
        out.append(f'<p class="note">{_e(note)}</p>')
    # Contents
    out.append('<nav><h2 id="contents">Contents</h2><ol>')
    for _, rep in fragments:
        out.append(f'<li><a href="#step-{_e(rep.step)}">{_e(rep.title)}</a></li>')
    out.append('<li><a href="#methods">Methods</a></li>')
    out.append('<li><a href="#appendix">Software and parameters</a></li></ol></nav>')
    # Steps
    n_table = n_figure = 0
    for folder, rep in fragments:
        out.append(f'<h2 id="step-{_e(rep.step)}">{_e(rep.title)}</h2>')
        for para in rep.summary:
            out.append(f"<p>{_e(para)}</p>")
        for table in rep.tables:
            n_table += 1
            out.append(_table_html(folder, table, n_table, rep.step))
        for fig in rep.figures:
            n_figure += 1
            out.append(_figure_html(folder, fig, n_figure, rep.step))
        for note in rep.notes:
            out.append(f'<p class="note">{_e(note)}</p>')
    # Methods
    records = unique_records(m for _, rep in fragments for m in rep.methods)
    out.append('<h2 id="methods">Methods</h2>')
    cited: list[str] = []
    if not records:
        out.append('<p class="note">No stage record was available for the methods text.</p>')
    for m in records:
        text, keys_ = sentence(m)
        out.append(f"<p>{_e(text)}</p>")
        cited += [k for k in keys_ if k not in cited]
    if cited:
        out.append("<h3>References</h3><ol>")
        for key in cited:
            ref, doi = CITATIONS[key]
            out.append(f'<li>{_e(ref)} <a href="https://doi.org/{_e(doi)}">doi:{_e(doi)}</a></li>')
        out.append("</ol>")
    # Appendix
    out.append('<h2 id="appendix">Software and parameters</h2>')
    if prov.stages:
        out.append('<h3>Stages of the run</h3><div class="tablewrap"><table><thead><tr>')
        out.append(
            "<th>stage</th><th>state</th><th>tool</th><th>completed (UTC)</th></tr></thead><tbody>"
        )
        for stage, state, tool, completed in prov.stages:
            out.append(
                f"<tr><td>{_e(stage)}</td><td>{_e(state)}</td><td>{_e(tool)}</td>"
                f"<td>{_e(completed)}</td></tr>"
            )
        out.append("</tbody></table></div>")
    versions = software_versions(records)
    if versions:
        out.append('<h3>Software versions</h3><div class="tablewrap"><table><thead><tr>')
        out.append("<th>software</th><th>version</th></tr></thead><tbody>")
        for tool, ver in versions:
            out.append(f"<tr><td>{_e(tool)}</td><td>{_e(ver)}</td></tr>")
        out.append("</tbody></table></div>")
    params = parameter_rows(records)
    if params:
        out.append('<h3>Parameters</h3><div class="tablewrap"><table><thead><tr>')
        out.append("<th>stage</th><th>parameter</th><th>value</th></tr></thead><tbody>")
        for stage, key, value in params:
            out.append(f"<tr><td>{_e(stage)}</td><td>{_e(key)}</td><td>{_e(value)}</td></tr>")
        out.append("</tbody></table></div>")
    if prov.command:
        out.append(f'<p class="meta">Command: <code>{_e(prov.command)}</code></p>')
    out.append("</main></body></html>")
    return "\n".join(out) + "\n"
