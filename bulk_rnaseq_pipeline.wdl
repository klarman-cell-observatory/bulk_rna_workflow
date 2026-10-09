version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/kco:bcl_convert/versions/12/plain-WDL/descriptor" as bcl_convert_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/jgould:bulk_rna_seq/versions/17/plain-WDL/descriptor" as bulk_rna_seq_wdl

struct FastqTableEntry {
  String sample_name
  Array[String] read1
  Array[String] read2
}

workflow bulk_rnaseq_pipeline {
  input {
    String experiment_name

    String input_bcl_directory
    String output_directory
    File input_xlsx
    File bcl_samplesheet_template
    String bcl_convert_software_version = "4.2.7"
    Boolean rc_i5 = true

    Boolean run_bcl_convert
    Boolean force_rerun_bcl_convert = false
    Int bcl_convert_preemptible = 0

    Boolean upload_fastq_table_to_terra = false
    String? terra_workspace_namespace
    String? terra_workspace_name

    Boolean output_genome_bam = false

    Boolean run_alignment
    String? existing_results_directory

    Float deseq2_padj_thresh = 0.05
    Float deseq2_lfc_thresh = 1.0
    String build_count_matrix_script = "/bulk_rna_seq/repo/scripts/build_count_matrix.R"
    String qc_plots_script = "/bulk_rna_seq/repo/scripts/qc_plots.R"
    String run_deseq2_script = "/bulk_rna_seq/repo/scripts/run_deseq2.R"

    String downstream_docker = "gcr.io/genomics-xavier/bulkrnaseq:latest"
  }

  String input_bcl_directory_gcs = if sub(input_bcl_directory, "^gs://.*", "MATCH") == "MATCH" then input_bcl_directory else "gs://" + input_bcl_directory
  String output_directory_gcs = if sub(output_directory, "^gs://.*", "MATCH") == "MATCH" then output_directory else "gs://" + output_directory
  String existing_results_directory_gcs = if defined(existing_results_directory) then (if sub(select_first([existing_results_directory]), "^gs://.*", "MATCH") == "MATCH" then select_first([existing_results_directory]) else "gs://" + select_first([existing_results_directory])) else ""

  call generate_bcl_samplesheet {
    input:
      experiment_name = experiment_name,
      input_xlsx = input_xlsx,
      bcl_samplesheet_template = bcl_samplesheet_template,
      bcl_convert_software_version = bcl_convert_software_version,
      rc_i5 = rc_i5,
      output_directory = output_directory_gcs,
      docker = downstream_docker
  }

  call detect_existing_bcl_output {
    input:
      output_directory = output_directory_gcs,
      docker = downstream_docker
  }

  Boolean do_bcl_convert = run_bcl_convert && (force_rerun_bcl_convert || !detect_existing_bcl_output.exists)

  if (do_bcl_convert) {
    call bcl_convert_wdl.run_bcl_convert as bcl_convert {
      input:
        input_bcl_directory = input_bcl_directory_gcs,
        output_directory = output_directory_gcs,
        sample_sheet = generate_bcl_samplesheet.sample_sheet_path,
        preemptible = bcl_convert_preemptible
    }
  }

  Boolean bcl_convert_done = defined(bcl_convert.fastqs)

  call build_fastq_table {
    input:
      output_directory = output_directory_gcs,
      input_xlsx = input_xlsx,
      experiment_name = experiment_name,
      bcl_convert_completion_marker = bcl_convert_done,
      docker = downstream_docker
  }

  if (upload_fastq_table_to_terra) {
    call upload_table_to_terra {
      input:
        fastq_table_tsv = build_fastq_table.fastq_table_tsv,
        workspace_namespace = select_first([terra_workspace_namespace]),
        workspace_name = select_first([terra_workspace_name]),
        docker = downstream_docker
    }
  }

  Array[FastqTableEntry] fastq_table = read_json(build_fastq_table.fastq_table_json)

  if (run_alignment) {
    scatter (e in fastq_table) {
      call bulk_rna_seq_wdl.bulk_rna_seq as bulk_rna_seq {
        input:
          sample_name = e.sample_name,
          read1 = e.read1,
          read2 = e.read2,
          reference = generate_bcl_samplesheet.reference_build,
          aligner = "star",
          output_genome_bam = output_genome_bam
      }
    }
  }

  if (!run_alignment) {
    call discover_existing_rsem_results {
      input:
        results_directory = existing_results_directory_gcs,
        sample_names = build_fastq_table.sample_names_ordered,
        docker = downstream_docker
    }
  }

  Array[File] rsem_gene_results =
    if run_alignment
    then select_first([bulk_rna_seq.rsem_gene, []])
    else select_first([discover_existing_rsem_results.gene_results, []])

  Array[File] rsem_isoform_results =
    if run_alignment
    then select_first([bulk_rna_seq.rsem_isoform, []])
    else select_first([discover_existing_rsem_results.isoform_results, []])

  Array[File] aligner_logs =
    if run_alignment
    then select_first([bulk_rna_seq.aligner_log, []])
    else select_first([discover_existing_rsem_results.aligner_logs, []])
  Array[String] sample_names_final = build_fastq_table.sample_names_ordered

  call analysis {
    input:
      gene_results = rsem_gene_results,
      aligner_logs = aligner_logs,
      sample_names = sample_names_final,
      metadata_all_csv = generate_bcl_samplesheet.metadata_all_csv,
      set_names = generate_bcl_samplesheet.set_names,
      set_metadata_csvs = generate_bcl_samplesheet.set_metadata_csvs,
      set_comparisons_jsons = generate_bcl_samplesheet.set_comparisons_jsons,
      build_count_matrix_script = build_count_matrix_script,
      qc_plots_script = qc_plots_script,
      run_deseq2_script = run_deseq2_script,
      deseq2_padj_thresh = deseq2_padj_thresh,
      deseq2_lfc_thresh = deseq2_lfc_thresh,
      docker = downstream_docker
  }

  call generate_report {
    input:
      experiment_name = experiment_name,
      reference_build = generate_bcl_samplesheet.reference_build,
      metadata_all_csv = generate_bcl_samplesheet.metadata_all_csv,
      metrics_json = analysis.metrics_json,
      alignment_rate_histogram = analysis.alignment_rate_histogram,
      volcano_pngs = analysis.volcano_pngs,
      set_names = generate_bcl_samplesheet.set_names,
      set_metadata_csvs = generate_bcl_samplesheet.set_metadata_csvs,
      output_directory = output_directory_gcs,
      docker = downstream_docker
  }

  call delocalize_outputs {
    input:
      output_directory = output_directory_gcs,
      copy_rsem_outputs = run_alignment,
      fastq_table_tsv = build_fastq_table.fastq_table_tsv,
      rsem_gene_results = rsem_gene_results,
      rsem_isoform_results = rsem_isoform_results,
      aligner_logs = aligner_logs,
      count_ensembl_csv = analysis.count_ensembl_csv,
      count_geneID_csv = analysis.count_geneID_csv,
      pca_plot = analysis.pca_plot,
      correlation_heatmap = analysis.correlation_heatmap,
      qc_flags_csv = analysis.qc_flags_csv,
      alignment_rate_histogram = analysis.alignment_rate_histogram,
      deseq2_results = analysis.deseq2_results,
      volcano_pngs = analysis.volcano_pngs,
      metrics_json = analysis.metrics_json,
      warnings_log = analysis.warnings_log,
      docker = downstream_docker
  }

  output {
    File bcl_convert_samplesheet_file = generate_bcl_samplesheet.bcl_convert_samplesheet_file
    String bcl_convert_samplesheet_gcs_path = generate_bcl_samplesheet.sample_sheet_path
    File fastq_table_tsv = build_fastq_table.fastq_table_tsv

    Array[File] rsem_gene_results_out = rsem_gene_results
    Array[File] rsem_isoform_results_out = rsem_isoform_results

    File count_matrix_geneID_csv = analysis.count_geneID_csv
    File pca_plot = analysis.pca_plot
    File correlation_heatmap = analysis.correlation_heatmap
    File qc_flags_csv = analysis.qc_flags_csv
    File alignment_rate_histogram = analysis.alignment_rate_histogram
    Array[File] deseq2_results = analysis.deseq2_results
    Array[File] volcano_pngs = analysis.volcano_pngs
    File metrics_summary = analysis.metrics_json
    File analysis_warnings_log = analysis.warnings_log

    String results_location = delocalize_outputs.results_location
    File report_pdf = generate_report.report_pdf
  }
}

task detect_existing_bcl_output {
  input {
    String output_directory
    String docker
  }

  command <<<
    set -uo pipefail
    if gsutil -q stat "~{sub(output_directory, "/$", "")}/bcl_fastqs.txt"; then
      echo "true" > exists.txt
    else
      echo "false" > exists.txt
    fi
  >>>

  output {
    Boolean exists = read_boolean("exists.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task generate_bcl_samplesheet {
  input {
    String experiment_name
    File input_xlsx
    File bcl_samplesheet_template
    String bcl_convert_software_version
    Boolean rc_i5
    String output_directory
    String docker
  }

  command <<<
    set -euo pipefail
    python3 <<'PYEOF'
import json
import pandas as pd

XLSX = "~{input_xlsx}"
TEMPLATE = "~{bcl_samplesheet_template}"
SW_VERSION = "~{bcl_convert_software_version}"
RC_I5 = "~{rc_i5}" == "true"

_RC = str.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s): return s.translate(_RC)[::-1]

def read_samplesheet(xlsx_path):
    names = pd.ExcelFile(xlsx_path).sheet_names
    by_lower = {n.strip().lower(): n for n in names}
    if "samplesheet" in by_lower:
        return pd.read_excel(xlsx_path, sheet_name=by_lower["samplesheet"])
    if len(names) == 1:
        return pd.read_excel(xlsx_path, sheet_name=names[0])
    raise SystemExit(f"No sheet named 'Samplesheet' and the workbook has multiple sheets {names}")

ss = read_samplesheet(XLSX)
required = ["Sample_Name", "index", "index2", "Sample_Group", "Organism"]
missing = [c for c in required if c not in ss.columns]
if missing:
    raise SystemExit(f"Samplesheet sheet is missing required column(s): {missing}")
ss["Sample_Name"] = ss["Sample_Name"].astype(str).str.strip()

rows = []
for _, row in ss.iterrows():
    i5 = revcomp(str(row["index2"]).strip()) if RC_I5 else str(row["index2"]).strip()
    rows.append(f'{row["Sample_Name"]},{str(row["index"]).strip()},{i5}')

with open(TEMPLATE) as f:
    template = f.read()
filled = (template
          .replace("{{BCLCONVERT_SOFTWARE_VERSION}}", SW_VERSION)
          .replace("{{SAMPLE_ROWS}}", "\n".join(rows)))
with open("~{experiment_name}_bcl_convert_samplesheet.csv", "w") as f:
    f.write(filled)
PYEOF

    gsutil cp "~{experiment_name}_bcl_convert_samplesheet.csv" \
        "~{sub(output_directory, "/$", "")}/inputs/~{experiment_name}_bcl_convert_samplesheet.csv"
    echo -n "~{sub(output_directory, "/$", "")}/inputs/~{experiment_name}_bcl_convert_samplesheet.csv" > sample_sheet_path.txt

    python3 <<'PYEOF'
import json
import pandas as pd

def read_samplesheet(xlsx_path):
    names = pd.ExcelFile(xlsx_path).sheet_names
    by_lower = {n.strip().lower(): n for n in names}
    if "samplesheet" in by_lower:
        return pd.read_excel(xlsx_path, sheet_name=by_lower["samplesheet"])
    if len(names) == 1:
        return pd.read_excel(xlsx_path, sheet_name=names[0])
    raise SystemExit(f"No sheet named 'Samplesheet' and the workbook has multiple sheets {names}")

XLSX = "~{input_xlsx}"
ss = read_samplesheet(XLSX)
ss["Sample_Name"] = ss["Sample_Name"].astype(str).str.strip()

organisms = ss["Organism"].astype(str).str.strip().str.lower().unique()
if len(organisms) != 1:
    raise SystemExit(f"Organism column must be uniform across all samples; found: {list(organisms)}")

ORGANISM_TO_REFERENCE = {"mouse": "GRCm38_ens93filt", "human": "GRCh38_ens93filt"}
organism = organisms[0]
if organism not in ORGANISM_TO_REFERENCE:
    raise SystemExit(f"Unrecognized Organism '{organism}' — expected one of {list(ORGANISM_TO_REFERENCE)}")
with open("reference_build.txt", "w") as f:
    f.write(ORGANISM_TO_REFERENCE[organism])

meta_all = ss[["Sample_Name", "Sample_Group"]]
meta_all.to_csv("sample_metadata_full.csv", index=False)

try:
    comp = pd.read_excel(XLSX, sheet_name="Comparisons")
    have_comparisons = True
except ValueError:
    have_comparisons = False

set_names = []
if not have_comparisons:
    meta_all.to_csv("set_metadata_all.csv", index=False)
    with open("set_comparisons_all.json", "w") as f:
        json.dump([], f)
    set_names.append("all")
else:
    req = ["comparison_set", "numerator_group", "reference_group"]
    missing_c = [c for c in req if c not in comp.columns]
    if missing_c:
        raise SystemExit(f"Comparisons sheet is missing required column(s): {missing_c}")
    if "exclude_samples" not in comp.columns:
        comp["exclude_samples"] = ""
    blank_set = comp["comparison_set"].isna() | (comp["comparison_set"].astype(str).str.strip() == "")
    if blank_set.any():
        bad_rows = comp[blank_set][["numerator_group", "reference_group"]].values.tolist()
        raise SystemExit(f"{blank_set.sum()} row(s) have a blank comparison_set: {bad_rows}")
    for set_name, sub in comp.groupby("comparison_set"):
        set_name = str(set_name)
        exclude = set()
        for v in sub["exclude_samples"].fillna(""):
            exclude.update(x.strip() for x in str(v).split(";") if x.strip())
        meta_set = meta_all[~meta_all["Sample_Name"].isin(exclude)]
        meta_set.to_csv(f"set_metadata_{set_name}.csv", index=False)
        pairs = [[str(r["numerator_group"]).strip(), str(r["reference_group"]).strip()] for _, r in sub.iterrows()]
        with open(f"set_comparisons_{set_name}.json", "w") as f:
            json.dump(pairs, f)
        set_names.append(set_name)

with open("set_names.txt", "w") as f:
    f.write("\n".join(set_names))
PYEOF
  >>>

  output {
    File bcl_convert_samplesheet_file = "~{experiment_name}_bcl_convert_samplesheet.csv"
    String sample_sheet_path = read_string("sample_sheet_path.txt")
    String reference_build = read_string("reference_build.txt")
    File metadata_all_csv = "sample_metadata_full.csv"
    Array[String] set_names = read_lines("set_names.txt")
    Array[File] set_metadata_csvs = glob("set_metadata_*.csv")
    Array[File] set_comparisons_jsons = glob("set_comparisons_*.json")
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

task build_fastq_table {
  input {
    String output_directory
    File input_xlsx
    String experiment_name
    Boolean bcl_convert_completion_marker
    String docker
  }

  command <<<
    set -euo pipefail

    for attempt in 1 2 3; do
      gsutil ls -r "~{sub(output_directory, "/$", "")}/**" > all_files.txt || true
      n_fastqs=$(grep -c '\.fastq\.gz$' all_files.txt || true)
      if [ "${n_fastqs:-0}" -gt 0 ]; then
        break
      fi
      echo "Attempt ${attempt}: no FASTQs visible yet under ~{output_directory} — waiting 30s and retrying"
      sleep 30
    done

    python3 <<'PYEOF'
import json
import pandas as pd

def read_samplesheet(xlsx_path):
    names = pd.ExcelFile(xlsx_path).sheet_names
    by_lower = {n.strip().lower(): n for n in names}
    if "samplesheet" in by_lower:
        return pd.read_excel(xlsx_path, sheet_name=by_lower["samplesheet"])
    if len(names) == 1:
        return pd.read_excel(xlsx_path, sheet_name=names[0])
    raise SystemExit(f"No sheet named 'Samplesheet' and the workbook has multiple sheets {names}")

ss = read_samplesheet("~{input_xlsx}")
sample_names = ss["Sample_Name"].astype(str).str.strip().tolist()

with open("all_files.txt") as f:
    all_paths = [l.strip() for l in f if l.strip().endswith(".fastq.gz")]

table = []
missing = []
for name in sample_names:
    matches = [p for p in all_paths if f"/{name}_" in p or f"/{name}." in p]
    r1 = sorted(p for p in matches if "_R1_" in p)
    r2 = sorted(p for p in matches if "_R2_" in p)
    if not r1 or not r2:
        missing.append(name)
        continue
    table.append({"sample_name": name, "read1": r1, "read2": r2})

if missing:
    raise SystemExit(f"No paired FASTQs found for {len(missing)} sample(s): {missing}")

with open("~{experiment_name}_fastq_table.json", "w") as f:
    json.dump(table, f)

with open("~{experiment_name}_sample_names.txt", "w") as f:
    f.write("\n".join(row["sample_name"] for row in table))

with open("~{experiment_name}_fastq_table.tsv", "w") as f:
    f.write(f'entity:~{experiment_name}_id\tread1\tread2\n')
    for row in table:
        f.write(f'{row["sample_name"]}\t{json.dumps(row["read1"])}\t{json.dumps(row["read2"])}\n')

print(f"{len(table)} samples")
PYEOF
  >>>

  output {
    File fastq_table_json = "~{experiment_name}_fastq_table.json"
    File fastq_table_tsv = "~{experiment_name}_fastq_table.tsv"
    Array[String] sample_names_ordered = read_lines("~{experiment_name}_sample_names.txt")
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

task upload_table_to_terra {
  input {
    File fastq_table_tsv
    String workspace_namespace
    String workspace_name
    String docker
  }

  command <<<
    set -euo pipefail
    python3 <<'PYEOF'
import firecloud.api as fapi

with open("~{fastq_table_tsv}") as f:
    tsv_text = f.read()

resp = fapi.upload_entities("~{workspace_namespace}", "~{workspace_name}", tsv_text)
print(resp.status_code, resp.text)
if resp.status_code >= 300:
    print("WARNING: Terra table upload failed.")
PYEOF
  >>>

  output {
    Boolean attempted = true
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task discover_existing_rsem_results {
  input {
    String results_directory
    Array[String] sample_names
    String docker
  }

  command <<<
    set -euo pipefail
    RESULTS_DIR="~{sub(results_directory, "/$", "")}/rsem"
    gsutil ls "${RESULTS_DIR}/**" > rsem_files.txt || true

    python3 <<'PYEOF'
import os
import subprocess
import sys
import json

samples_json_path = "~{write_json(sample_names)}"
with open(samples_json_path) as f:
    samples = json.load(f)

with open("rsem_files.txt") as f:
    paths = [line.strip() for line in f if line.strip()]

# Files are delocalized into results/rsem/. Match by exact sample-name prefix
# plus a recognized suffix to avoid accidental partial-name matches.
def matches_for(sample, kind):
    if kind == "gene":
        suffixes = (".genes.results", "_genes.results", ".gene.results")
    elif kind == "isoform":
        suffixes = (".isoforms.results", "_isoforms.results", ".isoform.results")
    else:
        suffixes = (".log",)
    found = []
    for path in paths:  
        base = os.path.basename(path)
        if any(base == sample + suffix or base.startswith(sample + suffix) for suffix in suffixes):
            found.append(path)
        elif kind == "log" and base.startswith(sample) and "Log.final.out" in base:
            found.append(path)
    return sorted(set(found))

out = {"gene": [], "isoform": [], "log": []}
problems = []
for sample in samples:
    for kind in out:
        found = matches_for(sample, kind)
        if len(found) != 1:
            problems.append(f"{sample}: expected exactly one {kind} file under results/rsem, found {len(found)}: {found}")
        else:
            out[kind].append(found[0])

if problems:
    print("Could not resolve existing RSEM results from the supplied results directory:", file=sys.stderr)
    print("\n".join(problems), file=sys.stderr)
    print("Expected files under results/rsem/ with sample-name-prefixed names, e.g. SAMPLE.genes.results, SAMPLE.isoforms.results, and SAMPLE.Log.final.out.", file=sys.stderr)
    sys.exit(1)

# Localize selected GCS files so WDL can expose them as File outputs.
for key, values in out.items():
    with open(key + "_paths.txt", "w") as f:
        f.write("\n".join(values))
PYEOF
  >>>

  output {
    Array[File] gene_results = read_lines("gene_paths.txt")
    Array[File] isoform_results = read_lines("isoform_paths.txt")
    Array[File] aligner_logs = read_lines("log_paths.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task analysis {
  input {
    Array[File] gene_results
    Array[File] aligner_logs
    Array[String] sample_names

    File metadata_all_csv
    Array[String] set_names
    Array[File] set_metadata_csvs
    Array[File] set_comparisons_jsons

    String build_count_matrix_script
    String qc_plots_script
    String run_deseq2_script
    Float deseq2_padj_thresh
    Float deseq2_lfc_thresh

    String docker
  }

  command <<<
    set -euo pipefail
    : > warnings.log

    mkdir -p data_dir/counts
    for f in ~{sep=" " gene_results}; do
      cp "$f" data_dir/counts/
    done
    Rscript ~{build_count_matrix_script} data_dir

    Rscript ~{qc_plots_script} data_dir "~{metadata_all_csv}"

    mkdir -p deseq2_results volcano_plots
    set_names=(~{sep=" " set_names})
    metadata_csvs=(~{sep=" " set_metadata_csvs})
    comparisons_jsons=(~{sep=" " set_comparisons_jsons})

    for i in "${!set_names[@]}"; do
      name="${set_names[$i]}"
      mkdir -p "deseq2_results/${name}"
      if [ "$(cat "${comparisons_jsons[$i]}")" == "[]" ]; then
        echo "No comparisons defined for set '${name}' — skipping." >> warnings.log
        continue
      fi
      if ! Rscript ~{run_deseq2_script} data_dir "${metadata_csvs[$i]}" "${comparisons_jsons[$i]}" "deseq2_results/${name}" 2>>warnings.log; then
        echo "WARNING: run_deseq2.R failed for comparison set '${name}'" >> warnings.log
      fi
    done

    python3 <<'PYEOF'
import glob, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

PADJ_THRESH = ~{deseq2_padj_thresh}
LFC_THRESH  = ~{deseq2_lfc_thresh}

def small_pos(arr):
    arr = arr[arr > 0]
    return arr.min() if len(arr) > 0 else 1e-300

summary = []
for path in sorted(glob.glob("deseq2_results/*/*_shrink.csv")):
    set_name = os.path.basename(os.path.dirname(path))
    comp = os.path.basename(path).replace("_shrink.csv", "")
    if "_vs_" not in comp:
        continue
    numer, denom = comp.split("_vs_", 1)
    try:
        df = pd.read_csv(path)
        df = df.rename(columns={df.columns[0]: "gene"}).set_index("gene", drop=False)
        df["padj"] = df["padj"].replace(0, small_pos(df["padj"].dropna().values))
        df["-log10padj"] = -np.log10(df["padj"])
        df["color"] = "black"
        sig = (df["padj"] <= PADJ_THRESH) & (df["log2FoldChange"].abs() > LFC_THRESH)
        df.loc[df["padj"] <= PADJ_THRESH, "color"] = "blue"
        df.loc[sig, "color"] = "red"

        fig, ax = plt.subplots(figsize=(7, 8))
        for color, grp in df.groupby("color"):
            ax.scatter(grp["log2FoldChange"], grp["-log10padj"], c=color, s=5, alpha=0.7, rasterized=True)
        ax.axhline(-np.log10(PADJ_THRESH), color="black", linestyle="--", lw=0.8)
        ax.axvline(LFC_THRESH, color="black", linestyle="--", lw=0.8)
        ax.axvline(-LFC_THRESH, color="black", linestyle="--", lw=0.8)
        ax.set(xlabel="log2 Fold Change", ylabel="-log10(padj)", title=f"{numer} vs {denom}")

        top = df[df["color"] == "red"].copy()
        top["dist"] = np.sqrt((top["log2FoldChange"].abs() - LFC_THRESH) ** 2 +
                              (top["-log10padj"] - (-np.log10(PADJ_THRESH))) ** 2)
        top = top.nlargest(20, "dist")
        texts = [ax.text(r["log2FoldChange"], r["-log10padj"], r["gene"], fontsize=6) for _, r in top.iterrows()]
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="red", lw=0.5))

        n_sig = int(sig.sum())
        ax.text(0.02, 0.98, f"{n_sig} sig. genes", transform=ax.transAxes, va="top", fontsize=9, color="red")
        out_dir = f"volcano_plots/{set_name}"
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(f"{out_dir}/{numer}_vs_{denom}_volcano.png", dpi=200, bbox_inches="tight")
        plt.close()
        summary.append({"set_name": set_name, "comparison": comp, "numerator": numer, "reference": denom, "n_significant": n_sig})
    except Exception as e:
        with open("warnings.log", "a") as wf:
            wf.write(f"WARNING: volcano plot failed for {set_name}/{comp}: {e}\n")

with open("volcano_summary.json", "w") as f:
    json.dump(summary, f)
PYEOF

    python3 <<'PYEOF'
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

names = "~{sep=',' sample_names}".split(",")
log_files = "~{sep=',' aligner_logs}".split(",")

aln_rows = []
for name, path in zip(names, log_files):
    pct = None
    try:
        with open(path) as fh:
            for line in fh:
                if "Uniquely mapped reads %" in line:
                    pct = float(line.strip().split("|")[-1].replace("%", "").strip())
                    break
    except Exception as e:
        with open("warnings.log", "a") as wf:
            wf.write(f"WARNING: could not parse aligner log for {name}: {e}\n")
    aln_rows.append({"sample": name, "pct_uniquely_mapped": pct})

aln_df = pd.DataFrame(aln_rows).dropna(subset=["pct_uniquely_mapped"])
alignment_summary = {}
if len(aln_df) > 0:
    mean_pct = float(aln_df.pct_uniquely_mapped.mean())
    std_pct = float(aln_df.pct_uniquely_mapped.std()) if len(aln_df) > 1 else 0.0
    outlier_thresh = mean_pct - 2 * std_pct
    outliers = aln_df[aln_df.pct_uniquely_mapped < outlier_thresh]["sample"].tolist()
    alignment_summary = {
        "mean_pct": mean_pct,
        "std_pct": std_pct,
        "min_pct": float(aln_df.pct_uniquely_mapped.min()),
        "min_sample": aln_df.loc[aln_df.pct_uniquely_mapped.idxmin(), "sample"],
        "max_pct": float(aln_df.pct_uniquely_mapped.max()),
        "max_sample": aln_df.loc[aln_df.pct_uniquely_mapped.idxmax(), "sample"],
        "low_mean_warning": mean_pct < 70,
        "outlier_threshold_pct": outlier_thresh,
        "outlier_samples": outliers,
    }
    if alignment_summary["low_mean_warning"]:
        with open("warnings.log", "a") as wf:
            wf.write(f"WARNING: mean alignment rate {mean_pct:.1f}% is below 70%\n")
    for s in outliers:
        with open("warnings.log", "a") as wf:
            wf.write(f"WARNING: {s} alignment rate is an outlier (>2 SD below run mean)\n")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(aln_df.pct_uniquely_mapped, bins=20)
    ax.set(xlabel="Uniquely mapped reads %", title="STAR Alignment Rate")
    os.makedirs("data_dir/Figures", exist_ok=True)
    plt.savefig("data_dir/Figures/alignment_rates.png", dpi=150, bbox_inches="tight")
    plt.close()

qc_flags_path = "data_dir/Figures/qc_flags.csv"
qc_flags = pd.read_csv(qc_flags_path).to_dict(orient="records") if os.path.exists(qc_flags_path) else []

with open("volcano_summary.json") as f:
    volcano_summary = json.load(f)

metrics = {
    "alignment_rates": aln_rows,
    "alignment_summary": alignment_summary,
    "qc_flags": qc_flags,
    "deseq2_volcano_summary": volcano_summary,
}
with open("metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)
PYEOF
  >>>

  output {
    File count_ensembl_csv = "data_dir/star_expected_count_ensembl.csv"
    File count_geneID_csv = "data_dir/star_expected_count_geneID.csv"
    File pca_plot = "data_dir/Figures/pca_plot.png"
    File correlation_heatmap = "data_dir/Figures/correlation_heatmap.png"
    File qc_flags_csv = "data_dir/Figures/qc_flags.csv"
    File alignment_rate_histogram = "data_dir/Figures/alignment_rates.png"
    Array[File] deseq2_results = glob("deseq2_results/*/*.csv")
    Array[File] volcano_pngs = glob("volcano_plots/*/*.png")
    File metrics_json = "metrics.json"
    File warnings_log = "warnings.log"
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 50 HDD"
  }
}

task generate_report {
  input {
    String experiment_name
    String reference_build
    File metadata_all_csv
    File metrics_json
    File alignment_rate_histogram
    Array[File] volcano_pngs
    Array[String] set_names
    Array[File] set_metadata_csvs
    String output_directory
    String docker
  }

  command <<<
    set -euo pipefail
    python3 <<'PYEOF'
import json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

EXPERIMENT = "~{experiment_name}"
REFERENCE = "~{reference_build}"
ORGANISM = "Mouse" if "GRCm" in REFERENCE else "Human" if "GRCh" in REFERENCE else REFERENCE

meta = pd.read_csv("~{metadata_all_csv}")
with open("~{metrics_json}") as f:
    metrics = json.load(f)

set_names = "~{sep=',' set_names}".split(",")
set_metadata_paths = "~{sep=',' set_metadata_csvs}".split(",")
set_metadata = {n: pd.read_csv(p) for n, p in zip(set_names, set_metadata_paths)}

volcano_paths = "~{sep=',' volcano_pngs}".split(",")
volcano_by_comparison = {}
for p in volcano_paths:
    import os
    fname = os.path.basename(p)
    if fname.endswith("_volcano.png"):
        volcano_by_comparison[fname[:-len("_volcano.png")]] = p

pdf = PdfPages("~{experiment_name}_report.pdf")

fig = plt.figure(figsize=(11, 8.5))
fig.text(0.1, 0.65, "Bulk RNA-seq Analysis", fontsize=28, weight="bold")
fig.text(0.1, 0.57, EXPERIMENT, fontsize=16)
fig.text(0.1, 0.50, f"{ORGANISM} ({REFERENCE})  ·  {len(meta)} samples  ·  {meta['Sample_Group'].nunique()} groups", fontsize=12)
fig.text(0.1, 0.44, "Alignment: STAR/RSEM (jgould bulk_rna_seq)  ·  DE: DESeq2", fontsize=12)
pdf.savefig(fig)
plt.close(fig)

aln = metrics.get("alignment_summary", {})
fig = plt.figure(figsize=(11, 8.5))
fig.text(0.07, 0.93, "STAR Alignment Rate — Distribution", fontsize=18, weight="bold")
if aln:
    ax_img = fig.add_axes([0.15, 0.25, 0.7, 0.6])
    ax_img.imshow(plt.imread("~{alignment_rate_histogram}"))
    ax_img.axis("off")
    stats_line = f"n = {len(metrics['alignment_rates'])} samples   |   mean {aln['mean_pct']:.1f}%   |   SD {aln['std_pct']:.1f}%   |   range {aln['min_pct']:.1f}% - {aln['max_pct']:.1f}%"
    fig.text(0.07, 0.17, stats_line, fontsize=11)
    if aln.get("low_mean_warning"):
        fig.text(0.07, 0.12, "Overall rate is below the typical 70-85% range for this reference/pipeline.", fontsize=11, color="firebrick")
pdf.savefig(fig)
plt.close(fig)

if aln:
    rates = sorted(metrics["alignment_rates"], key=lambda r: r["pct_uniquely_mapped"] or 0)
    outlier_set = set(aln.get("outlier_samples", []))
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.93, "STAR Alignment Rate — Per Sample", fontsize=18, weight="bold")
    fig.text(0.07, 0.89, "Uniquely mapped reads %, sorted ascending. * = outlier (>2 SD below run mean)", fontsize=10)
    half = (len(rates) + 1) // 2
    for col, chunk in enumerate([rates[:half], rates[half:]]):
        rows = [[r["sample"] + (" *" if r["sample"] in outlier_set else ""), f"{r['pct_uniquely_mapped']:.2f}%"] for r in chunk]
        ax = fig.add_axes([0.07 + col * 0.47, 0.05, 0.42, 0.78])
        ax.axis("off")
        tbl = ax.table(cellText=rows, colLabels=["Sample", "% Uniquely Mapped"], loc="center", cellLoc="left")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
    pdf.savefig(fig)
    plt.close(fig)

for comp in metrics.get("deseq2_volcano_summary", []):
    key = comp["comparison"]
    if key not in volcano_by_comparison:
        continue
    numer, denom = comp["numerator"], comp["reference"]
    set_name = comp["set_name"]
    set_meta = set_metadata.get(set_name, meta)

    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.93, f"{numer}  vs  {denom}", fontsize=18, weight="bold")
    ax_img = fig.add_axes([0.05, 0.08, 0.55, 0.78])
    ax_img.imshow(plt.imread(volcano_by_comparison[key]))
    ax_img.axis("off")

    y = 0.80
    fig.text(0.65, y, "Comparison", fontsize=12, weight="bold"); y -= 0.04
    fig.text(0.65, y, f"{numer} vs {denom}", fontsize=10); y -= 0.03
    fig.text(0.65, y, f"(log2FC > 0 = higher in {numer})", fontsize=9); y -= 0.07

    fig.text(0.65, y, "Samples included", fontsize=12, weight="bold"); y -= 0.04
    for group in (numer, denom):
        members = set_meta[set_meta["Sample_Group"] == group]["Sample_Name"].tolist()
        fig.text(0.65, y, f"{group}: {', '.join(members)}  (n={len(members)})", fontsize=9); y -= 0.03
    y -= 0.04

    fig.text(0.65, y, "Significant genes", fontsize=12, weight="bold"); y -= 0.04
    fig.text(0.65, y, f"{comp['n_significant']:,}", fontsize=16); y -= 0.03
    fig.text(0.65, y, "padj <= 0.05 and |log2 fold change| > 1", fontsize=8); y -= 0.06

    excluded = sorted(set(meta["Sample_Name"]) - set(set_meta["Sample_Name"]))
    if excluded:
        fig.text(0.65, y, f"Excluded from this comparison: {', '.join(excluded)}", fontsize=8, color="firebrick")
    else:
        fig.text(0.65, y, "All samples included - none excluded from this comparison.", fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)

fig = plt.figure(figsize=(11, 8.5))
fig.text(0.07, 0.93, "QC Notes", fontsize=18, weight="bold")
y = 0.85
if aln:
    fig.text(0.07, y, "STAR alignment", fontsize=12, weight="bold"); y -= 0.04
    outliers = aln.get("outlier_samples", [])
    line = f"Overall mean uniquely-mapped rate was {aln['mean_pct']:.1f}%."
    fig.text(0.07, y, line, fontsize=10); y -= 0.03
    if outliers:
        fig.text(0.07, y, f"Outlier samples (>2 SD below run mean): {', '.join(outliers)}", fontsize=10); y -= 0.03
    y -= 0.04

qc_flags = metrics.get("qc_flags", [])
if qc_flags:
    fig.text(0.07, y, "Library size / PCA QC", fontsize=12, weight="bold"); y -= 0.04
    for flag in qc_flags:
        fig.text(0.07, y, f"{flag.get('sample','')}: {flag.get('check','')} — {flag.get('detail','')}", fontsize=9); y -= 0.03
    y -= 0.04
else:
    fig.text(0.07, y, "Library size / PCA QC: no flags raised.", fontsize=10); y -= 0.06

pdf.savefig(fig)
plt.close(fig)

pdf.close()
PYEOF

    gsutil cp "~{experiment_name}_report.pdf" "~{sub(output_directory, "/$", "")}/results/"
  >>>

  output {
    File report_pdf = "~{experiment_name}_report.pdf"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

task delocalize_outputs {
  input {
    String output_directory
    Boolean copy_rsem_outputs
    File fastq_table_tsv
    Array[File] rsem_gene_results
    Array[File] rsem_isoform_results
    Array[File] aligner_logs
    File count_ensembl_csv
    File count_geneID_csv
    File pca_plot
    File correlation_heatmap
    File qc_flags_csv
    File alignment_rate_histogram
    Array[File] deseq2_results
    Array[File] volcano_pngs
    File metrics_json
    File warnings_log
    String docker
  }

  command <<<
    set -euo pipefail
    DEST="~{sub(output_directory, "/$", "")}/results"

    gsutil cp "~{fastq_table_tsv}" "$DEST/"
    if [ "~{copy_rsem_outputs}" = "true" ]; then
      gsutil -m cp ~{sep=" " rsem_gene_results} "$DEST/rsem/"
      gsutil -m cp ~{sep=" " rsem_isoform_results} "$DEST/rsem/"
      gsutil -m cp ~{sep=" " aligner_logs} "$DEST/rsem/"
    fi
    gsutil cp "~{count_ensembl_csv}" "~{count_geneID_csv}" "$DEST/"
    gsutil cp "~{pca_plot}" "~{correlation_heatmap}" "~{alignment_rate_histogram}" "$DEST/figures/"
    gsutil cp "~{qc_flags_csv}" "~{metrics_json}" "~{warnings_log}" "$DEST/"
    gsutil -m cp ~{sep=" " deseq2_results} "$DEST/deseq2/"
    gsutil -m cp ~{sep=" " volcano_pngs} "$DEST/figures/"

    echo -n "$DEST" > results_location.txt
  >>>

  output {
    String results_location = read_string("results_location.txt")
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}
