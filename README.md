# bam_depth

Compute per-sample read depth from BAM files with [mosdepth](https://github.com/brentp/mosdepth), split into:

- **genome-wide** mean depth
- **on-target** mean depth (exome capture regions)
- **off-target** mean depth (everything outside the capture regions)

Off-target depth is the complement of the target BED across the reference contigs. It can be skipped with `--on-target-only` when only on-target coverage matters (faster: one mosdepth run per BAM instead of two).

Scripts live in `bam_depth/`. The earlier GeRI-specific scripts are kept for reference in `old_scripts/` and are not used by this pipeline.

## Layout

```
configs/
  c4_b37.config         reference + target paths for GRCh37 on c4
  c4_b38.config         reference + target paths for GRCh38 on c4
bam_depth/
  bam_depth.sh          main script (scheduler-agnostic)
  c4_wrapper.sh         SLURM array wrapper for the c4 cluster
  summarize_depth.sh    compile per-sample results into one table
  rank_gene_variability.py  rank genes by coverage variability across samples; optional batch-effect test
old_scripts/            earlier GeRI-specific scripts (not used)
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
| `ADD_CHR` | `true` if the BED has no `chr` prefix (`1`, `2`, ...) but the BAMs/reference do (`chr1`, ...); default `false`. Cannot be combined with `STRIP_CHR` |

`REF_FASTA` is used only for its `.fai` (contig names and lengths), so it must use the same contig names as the BAMs. Check with `samtools view -H sample.bam | grep '^@SQ' | head`.

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
sbatch --array=1-${N}%100 c4_wrapper.sh ../configs/c4_b38.config /path/to/bams.txt /path/to/outdir

# on-target only
sbatch --array=1-${N}%100 c4_wrapper.sh ../configs/c4_b38.config /path/to/bams.txt /path/to/outdir --on-target-only

# after the array finishes
bash summarize_depth.sh /path/to/outdir
```

One array task processes one line of the BAM list. The wrapper's default resources are 2 CPUs, 8G and 10 minutes per task; adjust the `#SBATCH` lines for your data. `mosdepth` must be on `PATH` inside the job (module or conda env); the wrapper loads `bedtools2` and `samtools`.

## Running without SLURM

```bash
bash bam_depth.sh --config ../configs/c4_b38.config --bam-list bams.txt --outdir out            # all BAMs, in turn
bash bam_depth.sh --config ../configs/c4_b38.config --bam-list bams.txt --outdir out --line 3   # only line 3
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
| `--gene-bed FILE` | also compute mean depth per gene (see [Per-gene depth](#per-gene-depth)) |
| `--gene-thresholds LIST` | with `--gene-bed`: depths N (comma-separated, default `10,20,30`) for the fraction of bases at ≥Nx; `none` to skip |
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

## Per-gene depth

Pass `--gene-bed` to get the mean depth of each gene in each sample, e.g. from a canonical-transcript CDS BED:

```
1   69090   70008   OR4F5   0   +   ENSG00000186092   ENST00000335137   1
```

Columns used: 1-3 (chrom, start, end), 4 (gene symbol), and 7 (Ensembl gene ID, optional; `NA` if absent). One line per exon/CDS interval. The gene BED is separate from `EXOME_TARGET_BED`, uses the same `STRIP_CHR`/`ADD_CHR` handling, and intervals on contigs missing from the reference are dropped with a warning. The extra run costs about one more on-target mosdepth pass per BAM and works with or without `--on-target-only`.

```bash
sbatch --array=1-${N}%100 c4_wrapper.sh ../configs/c4_b37.config /path/to/bams.txt /path/to/outdir --on-target-only --gene-bed /path/to/genes.bed
bash summarize_depth.sh /path/to/outdir
```

Gene depth is the length-weighted mean over the gene's intervals, `sum(mean_i * len_i) / sum(len_i)`. The coverage fractions are `sum(bases_i at >= Nx) / sum(len_i)`, i.e. the share of the gene's bases covered by at least N reads. Intervals are assumed non-overlapping within a gene (use one transcript per gene); overlapping intervals would be counted twice. Per-interval means from mosdepth are rounded to 2 decimals.

Output: `<outdir>/tmp_depth/<sample>.gene_depth.tsv.gz` per sample, and after summarizing `<outdir>/gene_depth.tsv` (long format):

| Column | Meaning |
|---|---|
| `sample_id` | BAM name minus `.bam` |
| `gene`, `ensg` | gene symbol and Ensembl ID from the gene BED |
| `n_bases` | total bases in the gene's intervals |
| `mean_depth` | length-weighted mean depth |
| `frac_10x`, `frac_20x`, `frac_30x` | fraction of the gene's bases covered at ≥10x / 20x / 30x (columns follow `--gene-thresholds`; absent with `none`) |

## Ranking genes by variability and batch effects

`rank_gene_variability.py` is a standalone script (python3 with numpy, pandas, scipy) that reads `gene_depth.tsv` and ranks genes from most to least variable in coverage across samples:

```bash
python3 rank_gene_variability.py /path/to/outdir/gene_depth.tsv                        # variability only
python3 rank_gene_variability.py /path/to/outdir/gene_depth.tsv --batch-map batches.tsv  # plus batch effects
```

It loads the whole table, so on ~1,000+ samples run it on a compute node, not the login node.

**`gene_variability.tsv`**: one row per gene, ranked. `mean`, `sd`, `cv` are over samples of the raw metric. The `norm_*` columns are the same after rescaling each sample to the median sample's depth (the median gene depth in each sample sets its scale), so whole-sample depth differences don't make every gene look variable (use `--no-normalize` to turn this off). Columns also include the mean and sd of any `frac_*` columns. Genes are ranked by `--rank-by` (default `norm_cv`, i.e. sd/mean after normalization). Counting noise alone makes low-coverage genes have a high CV, so `cv_excess` (log CV minus the median log CV of genes with similar mean depth) is the better choice for finding genes that vary more than expected at their coverage level; use `--rank-by cv_excess` for that.

**Batch effects** (with `--batch-map`): a two-column file, `sample_id` and `batch` (tab, comma or space separated, header optional; a trailing `.bam` on the sample name is ignored). `sample_id` must match the `sample_id` in the depth table. For each gene the script tests whether normalized coverage differs between batches:

- `gene_batch_effects.tsv`: Kruskal-Wallis test (`kw_p`, `kw_fdr` = Benjamini-Hochberg), `eta2` / `eta2_adj` (fraction of the gene's variance explained by batch; adjusted for the number of batches), `fold_range` (highest / lowest batch mean), and the highest and lowest batches. Ranked by `eta2_adj`. `flag` is true when `kw_fdr < 0.05` and `eta2_adj >= 0.1` (`--fdr`, `--min-eta2`). With hundreds of samples nearly any gene gets a tiny p-value, so go by the effect size.
- `gene_batch_means.tsv`: mean coverage of each gene in each batch, same order, for plotting.

Batches with fewer than 3 samples are dropped (`--min-batch-size`), and samples missing from the map are excluded with a warning. The per-batch median sample depth is printed, since a batch-wide depth difference is removed by the normalization. A batch effect here means coverage differs by batch; if batch is confounded with something else (cohort, capture kit, ancestry), the script can't tell those apart.

## Planned

- Apptainer/Singularity image so other labs only need: a BAM list, the genome build, and the exome target BED.
- Generalized plotting script (the old one is tied to GeRI batches and GLIMPSE concordance).
