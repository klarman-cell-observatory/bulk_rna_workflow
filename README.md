# BULKRNASEQ Workflow

This repository is the end-to-end write of the standard BulkRNASeq workflow used within the Xavier Lab. Mainly containing the WDL and the files that constitute the pipeline.

Steps:
BCL_Convert -> RSEM Alignment -> Downstream Analysis & Figure Generation (QCs + Metrics)

Caveats for running the pipeline, there is the ability to rerun specific steps due to the removal of hard-checks and stops due to QC.

End-to-end run:

bulk_rnaseq_pipeline.run_bcl_convert = 'true'
bulk_rnaseq_pipeline.run_alignment = 'true'

BCL_convert only:

bulk_rnaseq_pipeline.run_bcl_convert = 'true'

Alignment, with figure generation:

bulk_rnaseq_pipeline.run_alignment = 'true'

No alignment, analysis steps only:

bulk_rnaseq_pipeline.run_alignment = 'true'
bulk_rnaseq_pipeline.existing_results_directory = /path/to/results/rsem/
