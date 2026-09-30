Here's the full file inline:

```wdl
version 1.0

## bulk_rnaseq_pipeline.wdl
##
## ── The FISS/data-table question, resolved ──────────────────────────────────────────────────
## `bulk_rna_seq` is called as a NESTED call inside this same Cromwell run (import + `call`), not
## launched as a separate Terra submission against a data table. That means:
##   - No `this.<column>` syntax anywhere — that's Terra method-config templating, only meaningful
##     when Terra itself launches a workflow against table rows. It doesn't apply inside a `call`.
##   - No FISS upload in the critical path — `build_fastq_table` hands the scatter below a plain
##     JSON manifest directly. Nothing waits on Terra having ingested anything.
##   - No retry-depth globbing over `call-run_rsem/attempt-N/...` in `analysis` — that pattern
##     only exists when you're scraping a *separate* submission's raw execution bucket. Here,
##     `bulk_rna_seq.rsem_gene` is just a clean `Array[File]` Cromwell hands back directly,
##     already resolved past whatever retries happened internally.
## `upload_fastq_table_to_terra` below is kept as an OPTIONAL, best-effort side output (so the
## table is still browsable in Terra's DATA tab / rerunnable by hand later) — nothing downstream
## reads it back. That's what breaks the circularity: "table exists for humans to look at" and
## "table is required input for the next step" are no longer the same thing.
##
## ── Checkpoints ──────────────────────────────────────────────────────────────────────────────
## Two breakpoints, one per expensive external call, each skippable so a rerun doesn't redo work
## that already succeeded:
##   1. Post-bcl_convert: `output_directory` is a path YOU choose, so it's checkable.
##      `detect_existing_bcl_output` does a cheap `gsutil stat` on it before anything else runs;
##      bcl_convert is skipped automatically if `bcl_fastqs.txt` is already there. Override with
##      `force_rerun_bcl_convert = true`, or hard-disable with `run_bcl_convert = false`.
##      `build_fastq_table` always runs regardless — it just reads whatever's at
##      `output_directory`, whether freshly written this run or already there from a previous one.
##   2. Post-bulk_rna_seq: jgould's outputs land in Cromwell's own submission directories, not a
##      path you chose or can predict, so there's nothing to check the way bcl_convert's output
##      can be. This one is explicit-only: `run_alignment = false` + the four `existing_rsem_*`/
##      `existing_sample_names` inputs substitute RSEM output you already have. Leaving those unset
##      while skipping fails immediately (via `select_first`) rather than silently doing nothing.
##
## ── Other open items, carried over ──────────────────────────────────────────────────────────
## - reference_build mapping: using GRCm38_ens93filt / GRCh38_ens93filt (matches what I verified
##   of bulk_rna_seq's actual accepted `reference` values) — confirm this is right, not GRCm39.
## - Both import URLs still unverified (Terra auth wall) — same caveat as every earlier draft.
import "https://api.firecloud.org/ga4gh/v1/tools/kco:bcl_convert/versions/12/plain-WDL/descriptor" as bcl_convert_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/jgould:bulk_rna_seq/versions/17/plain-WDL/descriptor" as bulk_rna_seq_wdl

struct FastqTableEntry {
  String sample_name
  Array[String] read1
  Array[String] read2
}

workflow bulk_rnaseq_pipeline {
  input {
    # ============================================================================================
    # General
    # ============================================================================================
    # Used as a naming/prefix variable across generated artifacts (the BCLConvert samplesheet,
    # the exported fastq table, the Terra entity-id column) so everything from one run is easy to
    # trace back to it — not just a Terra-table column name anymore.
    String experiment_name

    # ============================================================================================
    # Stage 1 — BCL Convert
    # ============================================================================================
    String input_bcl_directory
    String output_directory
    File input_xlsx # "Samplesheet" (+ optional "Comparisons") sheets — same shape as before
    File bcl_samplesheet_template
    String bcl_convert_software_version = "4.2.7"
    Boolean rc_i5 = true

    # Checkpoint 1: skip bcl_convert if output_directory already has results.
    Boolean run_bcl_convert = true # hard-disable regardless of what's detected
    Boolean force_rerun_bcl_convert = false # rerun even if output already exists

    # ============================================================================================
    # Stage 2 — fastq table
    # ============================================================================================
    Boolean upload_fastq_table_to_terra = false # optional, best-effort, decoupled — see note above
    # These stay as plain String inputs — a WDL/Cromwell has no built-in notion of "which Terra
    # workspace is this running in", so there's nothing for the WDL itself to look up at runtime.
    # BUT Terra's method-configuration UI (separate from this file — the step where you attach/
    # configure this workflow as a method and set its inputs) lets you bind an input's value
    # source to the literal `workspace.namespace` / `workspace.name` instead of typing a string.
    # Set these two inputs that way once, in Terra, and it auto-fills them from wherever the
    # workflow actually runs — no per-run typing, and the WDL stays portable to non-Terra engines
    # (where you'd just pass the values directly instead).
    String? terra_workspace_namespace # required only if upload_fastq_table_to_terra = true
    String? terra_workspace_name

    # ============================================================================================
    # Stage 3 — alignment (jgould:bulk_rna_seq)
    # ============================================================================================
    Boolean output_genome_bam = false

    # Checkpoint 2: skip alignment entirely and substitute RSEM output you already have. Explicit
    # only — see the header note on why this one can't be auto-detected the way bcl_convert's can.
    Boolean run_alignment = true
    Array[File]? existing_rsem_gene_results
    Array[File]? existing_rsem_isoform_results
    Array[File]? existing_aligner_logs
    Array[String]? existing_sample_names

    # ============================================================================================
    # Stage 4 — downstream analysis
    # ============================================================================================
    Float deseq2_padj_thresh = 0.05
    Float deseq2_lfc_thresh = 1.0
    # Defaults point at the scripts baked into the image (see the Dockerfile — cloned into
    # /bulk_rna_seq/repo/scripts/). Override any of these if you want to test a modified script
    # without rebuilding the image.
    String build_count_matrix_script = "/bulk_rna_seq/repo/scripts/build_count_matrix.R"
    String qc_plots_script = "/bulk_rna_seq/repo/scripts/qc_plots.R"
    String run_deseq2_script = "/bulk_rna_seq/repo/scripts/run_deseq2.R"

    # Same image for every task now (generate_bcl_samplesheet / build_fastq_table / analysis) —
    # simpler to manage than two images, at the cost of a heavier container for the two small
    # linking tasks than they strictly need.
    String downstream_docker = "gcr.io/genomics-xavier/bulkrnaseq:latest"
  }

  # ================================================================================================
  # Stage 1 — BCL Convert
  # ================================================================================================
  call generate_bcl_samplesheet {
    input:
      experiment_name = experiment_name,
      input_xlsx = input_xlsx,
      bcl_samplesheet_template = bcl_samplesheet_template,
      bcl_convert_software_version = bcl_convert_software_version,
      rc_i5 = rc_i5,
      output_directory = output_directory,
      docker = downstream_docker
  }

  call detect_existing_bcl_output {
    input:
      output_directory = output_directory,
      docker = downstream_docker
  }

  Boolean do_bcl_convert = run_bcl_convert && (force_rerun_bcl_convert || !detect_existing_bcl_output.exists)

  if (do_bcl_convert) {
    call bcl_convert_wdl.bcl_convert as bcl_convert {
      input:
        input_bcl_directory = input_bcl_directory,
        output_directory = output_directory,
        sample_sheet = generate_bcl_samplesheet.sample_sheet_path # String GCS path
    }
  }

  # Dependency-forcing only — waits for bcl_convert when it actually ran; inert otherwise, since
  # output_directory's contents are already in place from a previous run.
  # TODO: swap `bcl_convert.output_directory` for its real declared output once you can see it in
  # Terra's Outputs panel.
  String bcl_dependency = select_first([bcl_convert.output_directory, output_directory])

  # ================================================================================================
  # Stage 2 — fastq table (plain WDL data; the Terra upload below is a side-output only).
  # Always runs, whether do_bcl_convert was true this time or output_directory already had
  # results from a previous run.
  # ================================================================================================
  call build_fastq_table {
    input:
      output_directory = output_directory,
      input_xlsx = input_xlsx,
      experiment_name = experiment_name,
      bcl_convert_completion_marker = bcl_dependency,
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

  # ================================================================================================
  # Stage 3 — alignment. Direct nested call — no `this.`, no data table read-back. Skippable per
  # checkpoint 2 above; select_first below substitutes existing_* output when run_alignment=false.
  # ================================================================================================
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

  Array[File] rsem_gene_results = select_first([bulk_rna_seq.rsem_gene, existing_rsem_gene_results])
  Array[File] rsem_isoform_results = select_first([bulk_rna_seq.rsem_isoform, existing_rsem_isoform_results])
  Array[File] aligner_logs = select_first([bulk_rna_seq.aligner_log, existing_aligner_logs])
  Array[String] sample_names_final = if run_alignment then fastq_table.sample_name else select_first([existing_sample_names])

  # ================================================================================================
  # Stage 4 — downstream analysis: one task, one container, run to completion.
  # ================================================================================================
  call analysis {
    input:
      gene_results = rsem_gene_results, # clean Array[File]; no globbing needed
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

  output {
    # Record-keeping exports — the pipeline itself never reads these back; they're here so a
    # human (or a separate tool) has a durable, inspectable copy of what was generated.
    File bcl_convert_samplesheet_file = generate_bcl_samplesheet.bcl_convert_samplesheet_file
    String bcl_convert_samplesheet_gcs_path = generate_bcl_samplesheet.sample_sheet_path
    File fastq_table_tsv = build_fastq_table.fastq_table_tsv # side-artifact; nothing above reads this back

    Array[File] rsem_gene_results_out = rsem_gene_results
    Array[File] rsem_isoform_results_out = rsem_isoform_results

    File count_matrix_geneID_csv = analysis.count_geneID_csv
    File pca_plot = analysis.pca_plot
    File correlation_heatmap = analysis.correlation_heatmap
    File qc_flags_csv = analysis.qc_flags_csv
    Array[File] deseq2_results = analysis.deseq2_results
    Array[File] volcano_pngs = analysis.volcano_pngs
    File metrics_summary = analysis.metrics_json
    File analysis_warnings_log = analysis.warnings_log
  }
}

# =================================================================================================
# detect_existing_bcl_output — cheap existence check that drives checkpoint 1's auto-skip.
# =================================================================================================
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

# =================================================================================================
# generate_bcl_samplesheet — template-filled BCLConvert samplesheet, reference_build from the
# Organism column, and the metadata/comparisons config the `analysis` task needs later.
# =================================================================================================
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

ss = pd.read_excel(XLSX, sheet_name="Samplesheet")
required = ["Sample_Name", "index", "index2", "Sample_Group", "Organism"]
missing = [c for c in required if c not in ss.columns]
if missing:
    raise SystemExit(f"Samplesheet sheet is missing required column(s): {missing}")
ss["Sample_Name"] = ss["Sample_Name"].astype(str).str.strip()

# ---- BCLConvert samplesheet, filled from the template ----
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

XLSX = "~{input_xlsx}"
ss = pd.read_excel(XLSX, sheet_name="Samplesheet")
ss["Sample_Name"] = ss["Sample_Name"].astype(str).str.strip()

# ---- reference_build, from Organism. Must be uniform — one `reference` value serves the whole run. ----
organisms = ss["Organism"].astype(str).str.strip().str.lower().unique()
if len(organisms) != 1:
    raise SystemExit(f"Organism column must be uniform across all samples; found: {list(organisms)}")

# TODO: confirm these are the exact accepted `reference` values for jgould:bulk_rna_seq — this
# matches what I verified from the cumulus docs, NOT the GRCm39/GRCh38 shorthand from an earlier
# whiteboard note. If GRCm39 support is real and different, this mapping needs to change.
ORGANISM_TO_REFERENCE = {"mouse": "GRCm38_ens93filt", "human": "GRCh38_ens93filt"}
organism = organisms[0]
if organism not in ORGANISM_TO_REFERENCE:
    raise SystemExit(f"Unrecognized Organism '{organism}' — expected one of {list(ORGANISM_TO_REFERENCE)}")
with open("reference_build.txt", "w") as f:
    f.write(ORGANISM_TO_REFERENCE[organism])

# ---- metadata_all.csv: Sample_Name + Sample_Group, for QC/DESeq2 ----
meta_all = ss[["Sample_Name", "Sample_Group"]]
meta_all.to_csv("sample_metadata_full.csv", index=False)

# ---- Comparisons sheet (optional) ----
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
    # pandas' groupby silently drops rows with a NaN/blank group key — that would mean a
    # comparison you actually wanted just vanishes with no error. Catch it here instead.
    blank_set = comp["comparison_set"].isna() | (comp["comparison_set"].astype(str).str.strip() == "")
    if blank_set.any():
        bad_rows = comp[blank_set][["numerator_group", "reference_group"]].values.tolist()
        raise SystemExit(f"{blank_set.sum()} row(s) have a blank comparison_set and would be silently dropped: {bad_rows}")
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
    String sample_sheet_path = read_string("sample_sheet_path.txt") # GCS path, matches bcl_convert's String-typed input
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

# =================================================================================================
# build_fastq_table — lists the predetermined FASTQ output layout, matches by sample_name, emits
# JSON (used directly by the scatter above) and a Terra-entity-shaped TSV (side-artifact only).
# =================================================================================================
task build_fastq_table {
  input {
    String output_directory
    File input_xlsx
    String experiment_name
    String bcl_convert_completion_marker # unused value; forces ordering after bcl_convert
    String docker
  }

  command <<<
    set -euo pipefail

    # Cromwell already guarantees this task doesn't START until bcl_convert's own task has
    # EXITED (that's what bcl_convert_completion_marker's dependency edge is for) — the normal
    # WDL call-dependency mechanism, not polling. This retry is only a defensive margin against
    # the one thing that guarantee *doesn't* cover: bcl_convert's own internal implementation
    # returning before it's actually done writing to output_directory (can't rule this out
    # without seeing its descriptor). If every listing attempt below comes up short, that's the
    # signal something upstream is genuinely wrong, not just slow.
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

ss = pd.read_excel("~{input_xlsx}", sheet_name="Samplesheet")
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
    # Genuinely missing samples after the retries above — surface this clearly rather than
    # silently building a partial table.
    raise SystemExit(f"No paired FASTQs found for {len(missing)} sample(s) under ~{output_directory}: {missing}")

with open("~{experiment_name}_fastq_table.json", "w") as f:
    json.dump(table, f)

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
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

# =================================================================================================
# upload_table_to_terra — OPTIONAL, best-effort. Not in the critical path: nothing above waits on
# this, and nothing downstream reads its output back. Purely so the table is browsable/rerunnable
# by hand in Terra's DATA tab later, if you want that.
# =================================================================================================
task upload_table_to_terra {
  input {
    File fastq_table_tsv
    String workspace_namespace
    String workspace_name
    String docker # needs: pip install firecloud
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
    # Non-fatal on purpose: this task's failure shouldn't be able to fail the pipeline, since
    # nothing downstream depends on it. Log and move on.
    print("WARNING: Terra table upload failed — table not available in the DATA tab this run.")
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

# =================================================================================================
# analysis — count matrix -> QC -> DESeq2 per comparison set -> volcano plots -> metrics.
# No retry-depth globbing: gene_results/aligner_logs are already clean, resolved Files.
# =================================================================================================
task analysis {
  input {
    Array[File] gene_results
    Array[File] aligner_logs
    Array[String] sample_names

    File metadata_all_csv
    Array[String] set_names
    Array[File] set_metadata_csvs
    Array[File] set_comparisons_jsons

    String build_count_matrix_script # path inside the image, or an overriding File-staged path
    String qc_plots_script
    String run_deseq2_script
    Float deseq2_padj_thresh
    Float deseq2_lfc_thresh

    String docker
  }

  command <<<
    set -euo pipefail
    : > warnings.log

    # ---- 1. count matrix ----
    mkdir -p data_dir/counts
    for f in ~{sep=" " gene_results}; do
      cp "$f" data_dir/counts/
    done
    Rscript ~{build_count_matrix_script} data_dir

    # ---- 2. QC plots ----
    Rscript ~{qc_plots_script} data_dir "~{metadata_all_csv}"

    # ---- 3. DESeq2, one call per comparison set; a failed set is logged and skipped ----
    mkdir -p deseq2_results volcano_plots
    set_names=(~{sep=" " set_names})
    metadata_csvs=(~{sep=" " set_metadata_csvs})
    comparisons_jsons=(~{sep=" " set_comparisons_jsons})

    for i in "${!set_names[@]}"; do
      name="${set_names[$i]}"
      mkdir -p "deseq2_results/${name}"
      # run_deseq2.R's own JSON parsing doesn't handle an empty `[]` cleanly (it attempts one
      # bogus comparison instead of running zero) — skip the call outright rather than let that
      # produce a confusing spurious error in the logs.
      if [ "$(cat "${comparisons_jsons[$i]}")" == "[]" ]; then
        echo "No comparisons defined for set '${name}' — skipping DESeq2 for this set." >> warnings.log
        continue
      fi
      if ! Rscript ~{run_deseq2_script} data_dir "${metadata_csvs[$i]}" "${comparisons_jsons[$i]}" "deseq2_results/${name}" 2>>warnings.log; then
        echo "WARNING: run_deseq2.R failed for comparison set '${name}'" >> warnings.log
      fi
    done

    # ---- 4. volcano plots, from whatever *_shrink.csv files exist ----
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

    # ---- 5. metrics ----
    python3 <<'PYEOF'
import json, os
import pandas as pd

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

qc_flags_path = "data_dir/Figures/qc_flags.csv"
qc_flags = pd.read_csv(qc_flags_path).to_dict(orient="records") if os.path.exists(qc_flags_path) else []

with open("volcano_summary.json") as f:
    volcano_summary = json.load(f)

metrics = {"alignment_rates": aln_rows, "qc_flags": qc_flags, "deseq2_volcano_summary": volcano_summary}
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
```