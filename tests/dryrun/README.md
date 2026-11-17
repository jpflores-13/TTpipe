# tests/dryrun: placeholder fixtures

Everything in this folder is a **placeholder**. The files exist only so that a Snakemake dry run
(`snakemake -n`) can build the job graph and find its inputs. They are **not real data** and no tool is ever
executed on them:

- `R1.fastq.gz`, `R2.fastq.gz`: one made-up read each (shared by all samples).
- `genome.fa`, `genes.gtf`: a 80 bp "chromosome" and two made-up genes (one tagged rRNA).
- `star_index/`: a fake STAR index directory (empty `SA`, a `chrNameLength.txt`); it is not a usable index.
- `samplesheet.txt`: 2 conditions (DMSO, 2H) x 2 replicates.
- `counts_placeholder.txt`: invented numbers in the count-matrix layout, for the analysis-only dry run.
- `placeholder.gmt`: one made-up gene set, to switch on the fgsea jobs.
- `config.yaml` (full run), `config_processing_only.yaml`, `config_analysis_only.yaml`: override paths and
  `outputDir` on top of `config/TTpipeConfig.yaml`.

Run from the repository root, e.g.

    snakemake -n -s workflows/TTpipe.snakefile --cores 1 --configfile config/TTpipeConfig.yaml tests/dryrun/config.yaml

`tests/test_dryrun.py` does this for several configurations. Results of the synthetic QC-report tests
(`tests/synth.py`) are likewise invented and carry a visible "synthetic" banner.
