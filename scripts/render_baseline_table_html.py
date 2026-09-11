#!/usr/bin/env python3
"""Render the aggregate verified Table 1 CSV as a self-contained HTML report.

No patient-level data are read or emitted. Usage:
python3 scripts/render_baseline_table_html.py <input_csv> <output_html>
"""

from __future__ import annotations

import csv
import html
import sys
from pathlib import Path


def esc(value: str | None) -> str:
    return html.escape(value or "")


def find_row(rows: list[dict[str, str]], variable: str) -> dict[str, str]:
    return next(row for row in rows if row["variable"] == variable)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: render_baseline_table_html.py <input_csv> <output_html>")

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    with input_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows or rows[0].get("variable") != "Patients, n":
        raise ValueError("Expected the verified Table 1 CSV with a Patients, n header row.")

    n_row = rows[0]
    death_row = find_row(rows, "60-day mortality, n (%)")
    stat_cards = [
        ("Eligible cohort", n_row["overall"]),
        ("Non-resistant", n_row["Non-resistant (A=0)"]),
        ("Recorded carbapenem-resistant", n_row["Carbapenem-resistant (A=1)"]),
        ("60-day mortality", death_row["overall"]),
    ]

    body_rows: list[str] = []
    for row in rows:
        row_type = row["row_type"]
        if row_type == "section":
            body_rows.append(f'<tr class="section"><td colspan="5">{esc(row["variable"])}</td></tr>')
            continue
        row_class = " outcome" if row["variable"] == "60-day mortality, n (%)" else ""
        body_rows.append(
            "<tr class=\"data%s\">"
            "<th scope=\"row\">%s</th>"
            "<td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
            % (
                row_class,
                esc(row["variable"]),
                esc(row["overall"]),
                esc(row["Non-resistant (A=0)"]),
                esc(row["Carbapenem-resistant (A=1)"]),
                esc(row["Max pairwise SMD"]),
            )
        )

    cards = "".join(
        f'<div class="stat"><span>{esc(label)}</span><strong>{esc(value)}</strong></div>'
        for label, value in stat_cards
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>REGARD-VAP | Verified baseline Table 1</title>
  <style>
    :root {{ --navy:#093b56; --teal:#087b85; --aqua:#eaf4f5; --ink:#172b36; --muted:#60717b; --line:#dbe4e7; --warm:#fbefea; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; color:var(--ink); background:#f4f7f8; font-family:Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    .hero {{ background:linear-gradient(125deg,var(--navy),#0b6172); color:#fff; padding:52px max(28px,calc((100vw - 1200px)/2)); }}
    .eyebrow {{ margin:0 0 10px; color:#bce8e8; font-weight:800; font-size:.78rem; letter-spacing:.12em; text-transform:uppercase; }}
    h1 {{ max-width:1000px; margin:0; font-size:clamp(1.8rem,4vw,2.65rem); line-height:1.12; letter-spacing:-.035em; }}
    .subtitle {{ max-width:950px; margin:15px 0 0; color:#d9eff1; font-size:1rem; line-height:1.55; }}
    main {{ max-width:1200px; margin:0 auto; padding:0 28px 48px; }}
    .stats {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:14px; margin-top:-23px; }}
    .stat {{ min-height:92px; padding:18px 20px; border:1px solid var(--line); border-radius:14px; background:#fff; box-shadow:0 10px 30px rgba(9,59,86,.08); }}
    .stat span {{ display:block; color:var(--muted); font-size:.78rem; font-weight:750; line-height:1.25; }}
    .stat strong {{ display:block; margin-top:7px; color:var(--navy); font-size:1.55rem; letter-spacing:-.03em; }}
    .panel {{ margin-top:26px; border:1px solid var(--line); border-radius:14px; overflow:hidden; background:#fff; box-shadow:0 8px 26px rgba(9,59,86,.05); }}
    .panel-head {{ padding:23px 26px 20px; border-bottom:1px solid var(--line); }}
    .panel-head h2 {{ margin:0; color:var(--navy); font-size:1.22rem; letter-spacing:-.02em; }}
    .panel-head p {{ margin:7px 0 0; color:var(--muted); font-size:.87rem; line-height:1.45; }}
    .table-scroll {{ overflow-x:auto; }}
    table {{ width:100%; min-width:900px; border-collapse:collapse; font-size:.88rem; }}
    thead th {{ padding:15px 14px; background:var(--navy); color:#fff; border-right:1px solid rgba(255,255,255,.16); text-align:right; font-size:.76rem; line-height:1.25; }}
    thead th:first-child {{ width:38%; text-align:left; }}
    tbody th, tbody td {{ padding:10px 14px; border-bottom:1px solid var(--line); vertical-align:top; }}
    tbody th {{ font-weight:650; text-align:left; }}
    tbody td {{ text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap; }}
    tbody tr.data:nth-child(odd) {{ background:#f8fbfb; }}
    tbody tr.section td {{ padding:9px 14px; background:var(--teal); color:#fff; border-bottom:0; font-size:.78rem; font-weight:800; letter-spacing:.025em; text-align:left; }}
    tbody tr.outcome {{ background:var(--warm) !important; }}
    .notes {{ padding:19px 26px 22px; color:var(--muted); background:#fbfcfc; font-size:.8rem; line-height:1.55; }}
    .notes strong {{ color:var(--ink); }}
    @media (max-width:760px) {{ .hero {{ padding:38px 22px; }} main {{ padding:0 14px 32px; }} .stats {{ grid-template-columns:repeat(2,minmax(0,1fr)); gap:10px; margin-top:-16px; }} .stat {{ min-height:78px; padding:14px; }} .stat strong {{ font-size:1.25rem; }} .panel-head {{ padding:19px 18px; }} }}
  </style>
</head>
<body>
  <header class="hero">
    <p class="eyebrow">REGARD-VAP · Study population and baseline characteristics</p>
    <h1>Baseline and onset characteristics by recorded carbapenem-resistance status</h1>
    <p class="subtitle">All eligible participants with recorded baseline exposure (A) and 60-day mortality outcome. This report uses only verified fields available in the current longitudinal analysis extract.</p>
  </header>
  <main>
    <section class="stats" aria-label="Cohort summary">{cards}</section>
    <section class="panel" aria-label="Complete baseline table">
      <div class="panel-head"><h2>Table 1</h2><p>Values are median [IQR] or n (%). The last column is the absolute standardized difference between A=0 and A=1.</p></div>
      <div class="table-scroll"><table>
        <thead><tr><th scope="col">Characteristic</th><th scope="col">Overall<br>N = {esc(n_row["overall"])}</th><th scope="col">Non-resistant<br>A = 0; N = {esc(n_row["Non-resistant (A=0)"])}</th><th scope="col">Recorded carbapenem-resistant<br>A = 1; N = {esc(n_row["Carbapenem-resistant (A=1)"])}</th><th scope="col">Absolute SMD</th></tr></thead>
        <tbody>{''.join(body_rows)}</tbody>
      </table></div>
      <div class="notes"><strong>Interpretation:</strong> the table describes the analysed cohort; the recorded resistance field is the exposure used in the current model. 3GCR subgroup labels, WGS variables and plasmid variables are intentionally excluded because they are not confirmed in the present analysis extract. SMD is descriptive and is not a hypothesis-test p value.</div>
    </section>
  </main>
</body>
</html>
"""
    output_path.write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
