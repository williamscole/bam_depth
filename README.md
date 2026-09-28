# bam_depth

Compute per-sample read depth from BAM files with [mosdepth](https://github.com/brentp/mosdepth), split into:

- **genome-wide** mean depth
- **on-target** mean depth (exome capture regions)
- **off-target** mean depth (everything outside the capture regions)

Off-target depth is the complement of the target BED across the reference contigs. It can be skipped with `--on-target-only` when only on-target coverage matters (faster: one mosdepth run per BAM instead of two).

Scripts live in `bam_depth/`. The earlier GeRI-specific scripts are kept for reference in `old_scripts/` and are not used by this pipeline.

## Layout

```
bam_depth/
  bam_depth.sh          main script (scheduler-agnostic)
  c4_wrapper.sh         SLURM array wrapper for the c4 cluster
  summarize_depth.sh    compile per-sample results into one table
  configs/
    c4_b37.config       reference + target paths for GRCh37 on c4
    c4_b38.config       reference + target paths for GRCh38 on c4
```

## Requirements

`mosdepth` and `bedtools` on `PATH`. The reference FASTA must have a `.fai` index next to it (`samtools faidx ref.fa`). Input is BAM only (no CRAM yet).

### Installing mosdepth

mosdepth is a compiled binary and is not on PyPI, so a plain Python venv (`pip install`) cannot install it. Two options:

**Conda (recommended).** On c4, load the miniforge module (check `module avail miniforge` for the version), then create an env:

```bash
module load miniforge/<version>
conda create -n bam_depth -c conda-forge -c bioconda mosdepth bedtools samtools
conda activate bam_depth
mosdepth --version
```

Activate the env in the wrapper before the call to `bam_depth.sh` (add `module load miniforge/<version>` and `conda activate bam_depth`), or activate it before running `bam_depth.sh` directly.

**Static binary (no conda).** Download `mosdepth` from the [releases page](https://github.com/brentp/mosdepth/releases), `chmod +x` it, and put it somewhere on `PATH`. It has no dependencies.

## Config files

Plain bash `KEY="value"` files that `bam_depth.sh` sources:

| Key | Meaning |
|---|---|
| `BUILD` | Genome build label (`b37`, `b38`, ...) |
| `REF_FASTA` | Reference FASTA (must have `.fai`) |
| `EXOME_TARGET_BED` | Exome target BED, `.bed` or `.bed.gz` (first 3 columns used) |
| `STRIP_CHR` | `true` if the BED has `chr`-prefixed contigs but the BAMs/reference do not (e.g. hg19 BED with b37 BAMs); default `false` |

Contigs in the BED that are not in the reference `.fai` are dropped. The script errors if no target regions remain (usually a `STRIP_CHR` or build mismatch).

## Inputs

- **BAM list**: text file, one BAM path per line.
- **Sample ID**: derived from the BAM file name (minus `.bam`), so BAM file names should be unique across the list.
- **Output dir**: per-sample results go to `<outdir>/tmp_depth/`; the summary goes to `<outdir>/depth_summary.tsv`.

## Running on c4

From `bam_depth/` (the wrapper finds `bam_depth.sh` via the submit directory, or set `BAM_DEPTH_DIR`):

```bash
cd bam_depth
mkdir -p logs
N=$(wc -l < /path/to/bams.txt)

# on- and off-target
sbatch --array=1-${N}%100 c4_wrapper.sh configs/c4_b38.config /path/to/bams.txt /path/to/outdir

# on-target only
sbatch --array=1-${N}%100 c4_wrapper.sh configs/c4_b38.config /path/to/bams.txt /path/to/outdir --on-target-only

# after the array finishes
bash summarize_depth.sh /path/to/outdir
```

One array task processes one line of the BAM list. The wrapper's default resources are 2 CPUs, 8G and 10 minutes per task; adjust the `#SBATCH` lines for your data. `mosdepth` must be on `PATH` inside the job (module or conda env); the wrapper loads `bedtools2` and `samtools`.

## Running without SLURM

```bash
bash bam_depth.sh --config configs/c4_b38.config --bam-list bams.txt --outdir out            # all BAMs, in turn
bash bam_depth.sh --config configs/c4_b38.config --bam-list bams.txt --outdir out --line 3   # only line 3
bash summarize_depth.sh out
```

`bam_depth.sh` options:

| Option | Meaning |
|---|---|
| `--config FILE` | config file (required) |
| `--bam-list FILE` | BAM list (required) |
| `--outdir DIR` | output directory (required) |
| `--line N` | process only line N (1-based) of the list; omit to process all |
| `--on-target-only` | skip the off-target run |
| `--threads T` | mosdepth threads (default `$SLURM_CPUS_PER_TASK` or 2) |

## Output

`<outdir>/tmp_depth/` holds, per sample, `<sample>_on.mosdepth.summary.txt` and `<sample>_on.mosdepth.global.dist.txt`, plus the `_off` equivalents unless `--on-target-only` was used. Per-base output is not written (`mosdepth -n`).

`<outdir>/depth_summary.tsv`:

| Column | Source |
|---|---|
| `sample_id` | BAM name minus `.bam` |
| `mean_genome` | `total` row of the on-target mosdepth summary (whole genome) |
| `mean_on` | `total_region` row of the on-target summary |
| `mean_off` | `total_region` row of the off-target summary; `NA` if not computed |

Each task rebuilds the on/off-target BEDs in its own scratch space (`$TMPDIR`, else `/tmp`), so parallel tasks never write to shared files.

## Planned

- Apptainer/Singularity image so other labs only need: a BAM list, the genome build, and the exome target BED.
- Generalized plotting script (the old one is tied to GeRI batches and GLIMPSE concordance).
