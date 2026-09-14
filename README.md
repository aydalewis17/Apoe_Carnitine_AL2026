# Apoe_Carnitine_AL2026
# ApoE Carnitine MetaG Pipeline

Modular, step-by-step metagenomics pipeline that takes raw paired-end reads from multiple
samples and produces a dereplicated, taxonomy-annotated MAG catalog with between-sample
relative abundance estimates.

Written by AL and Claude Sonnet 4.6, incorporating steps from the WCRC Multi-Omics
Crash Course workshop and KKAmundson's MetaG Processing LOOPS (April 2026).

---

## Overview

The pipeline runs in three phases:

**Phase 1 — Per-sample assembly and binning (steps 00-08)**
Each of the 13 ApoE samples is assembled, binned, and quality-assessed independently.
MQ/HQ MAGs (completeness >=50%, contamination <=10%) are copied to each sample's MAGs/
directory. SingleM pipe on raw reads also runs during this phase (no MAG dependency).

**Phase 2 — Co-assembly, iterative binning, and MAG pooling (steps 02b-02e)**
Five control samples and five carnitine samples are co-assembled to recover low-abundance
organisms. Two rounds of iterative assembly then recover organisms that escaped both
per-sample and co-assembly binning. All MQ/HQ MAGs from all strategies are pooled into
ALL_MQHQ_MAGs/ for DRAM2 annotation and dRep.

**Phase 3 — Cross-sample steps (steps 09-11 + SingleM appraise)**
Pooled MAGs are annotated with DRAM2, dereplicated at 99% ANI, and used as a shared
mapping reference for between-sample abundance estimation. SingleM appraise compares
the MAG catalog against read-based community profiles to quantify catalog completeness.

---

## Important note on pipeline execution order for this project

**Phase 3 was run twice.**

After Phase 1 was complete for all 13 samples, Phase 3 (09_cross_drep, 10_bowtie2_map,
11_coverm) was run once on per-sample MAGs only (693 MAGs → 162 primary clusters after
dRep) before Phase 2 was designed and executed. This was a deliberate interim step while
co-assembly and iterative assembly strategies were being planned.

Phase 2 (02b through 02e) was subsequently run, adding 245 additional MQ/HQ MAGs from
co-assembly and iterative strategies (938 total pre-dRep). Phase 3 must be re-run on the
complete ALL_MQHQ_MAGs/ pool to produce the final MAG catalog used for all downstream
analyses. The first dRep/Bowtie2/CoverM outputs from the per-sample-only run are
preserved on disk but should not be used for final analyses.

---

## Directory structure

```
/home/projects-phoenix/ApoE_Carnitine/MetaG/
│
├── scripts/                              ← all .sh scripts + sample lists live here
│   ├── 00_config.sh
│   ├── 00b_init_metadata.sh
│   ├── 00c_raw_qc.sh
│   ├── 01_trim.sh
│   ├── 02_assemble.sh
│   ├── 02b_coassemble.sh
│   ├── 02c_coassemble_bin.sh
│   ├── 02d_iterative_assembly.sh
│   ├── 02e_pool_mqhq_mags.sh
│   ├── 03_filter_contigs.sh
│   ├── 04_map_to_contigs.sh
│   ├── 05_metabat.sh
│   ├── 06_checkm2.sh
│   ├── 07_gtdbtk.sh
│   ├── 08_cleanup.sh
│   ├── 09_cross_drep.sh
│   ├── 10_bowtie2_map.sh
│   ├── 11_coverm.sh
│   ├── sample_list.txt                   ← 13 ApoE samples (no blanks/controls)
│   └── apoe_sample_list.txt              ← ApoE-specific list for cross-sample steps
│
├── coassemblies/scripts/
│   ├── control_list.txt                  ← 5 control sample names
│   └── carnitine_list.txt                ← 5 carnitine sample names
│
├── raw_reads/
├── fastqc_raw/
│   └── multiqc_raw/
├── pipeline_metadata.tsv
├── sample_MAG_database_mapping_summary.tsv
├── possort_bam_list.txt
│
├── {SAMPLE}/                             ← one per sample (x13)
│   ├── trimmed_reads/
│   │   └── fastqc/
│   ├── megahit_out/
│   │   ├── final.contigs.fa
│   │   ├── {SAMPLE}_2500.fa
│   │   ├── {SAMPLE}_2500.fa.metabat-bins/
│   │   └── checkm2_v1.1.0_fa/
│   │       └── quality_report.tsv
│   ├── MAGs/                             ← MQ/HQ MAGs from per-sample assembly
│   ├── gtdb_v2.7.0_r232/
│   ├── MAG_db_mapping/                   ← POSSORT BAMs from 10_bowtie2_map.sh
│   └── singlem_output/                   ← SingleM per-sample outputs
│       ├── {SAMPLE}_reads.otu_table.csv
│       ├── {SAMPLE}_reads.archive.otu_table.json.gz
│       ├── {SAMPLE}_reads.profile.tsv
│       ├── {SAMPLE}_prokaryotic_fraction.tsv
│       ├── {SAMPLE}_phylum_relabun.csv
│       ├── {SAMPLE}_krona.html
│       └── {SAMPLE}_appraise_unrecovered.csv  ← deferred until after dRep
│
├── coassemblies/
│   ├── MQ_HQ_genome_db_iter1.fa          ← pre-dRep pool used for iter1 BBMap filter
│   ├── MQ_HQ_genome_db_iter2.fa          ← pre-dRep pool used for iter2 BBMap filter
│   │
│   ├── control/
│   │   ├── concatenated_reads/           ← 02b: pooled R1/R2 from 5 control samples
│   │   ├── megahit_out/                  ← 02b: co-assembly
│   │   ├── checkm2_v1.1.0_fa/            ← 02c: CheckM2 on co-assembly bins
│   │   ├── MAGs/                         ← 02c: co-assembly MQ/HQ MAGs (83)
│   │   ├── iter1/                        ← 02d: iterative round 1
│   │   │   ├── unmapped_reads/
│   │   │   ├── megahit_out/
│   │   │   └── checkm2_v1.1.0_fa/
│   │   ├── MAGs_iter1/                   ← 02d: iter1 MQ/HQ MAGs (14)
│   │   ├── iter2/                        ← 02d: iterative round 2
│   │   │   ├── unmapped_reads/
│   │   │   ├── megahit_out/
│   │   │   └── checkm2_v1.1.0_fa/
│   │   └── MAGs_iter2/                   ← 02d: iter2 MQ/HQ MAGs (10)
│   │
│   └── carnitine/                        ← identical structure to control/
│       ├── MAGs/                         ← 102 MQ/HQ MAGs
│       ├── MAGs_iter1/                   ← 22 MQ/HQ MAGs
│       └── MAGs_iter2/                   ← 14 MQ/HQ MAGs
│
├── ALL_MQHQ_MAGs/                        ← 02e: all 938 pre-dRep MQ/HQ MAGs
│                                            input for DRAM2 and 09_cross_drep.sh
│
├── cross_sample_dRep/                    ← 09: dRep output (run twice — see note above)
│   └── dereplicated_genomes/
│       └── genome_renamed/
│
├── cross_sample_bowtie_DB/               ← 10: shared Bowtie2 index + FASTA
│
├── coverm_output/                        ← 11: between-sample abundance tables
│   ├── coverm_reads_per_base.txt
│   ├── coverm_min75.txt
│   └── coverm_trimmed_mean.txt
│
└── singlem_cross_sample/                 ← SingleM on dRep MAG catalog (post-09)
    ├── all_derep_MAGs.otu_table.csv
    └── all_derep_MAGs.profile.tsv
```

---

## Setup

### 1. Edit `00_config.sh`

This is the only file you should need to edit. Everything else derives from it.

At minimum update:
- `PROJECT_DIR` — root directory for the project
- `RAW_DIR` — where raw `.fastq.gz` files live
- `SCRIPTS_DIR` — where all scripts and `sample_list.txt` live
- `--mail-user` in every script header

### 2. Fix the hardcoded config path in every script

SLURM copies scripts to a temp directory at runtime, breaking `$(dirname "$0")`.
Every script has the config path hardcoded. Update all at once:

```bash
cd /your/scripts/directory
sed -i 's|/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts|/your/scripts/directory|g' \
    00b_init_metadata.sh 00c_raw_qc.sh 01_trim.sh 02_assemble.sh \
    02b_coassemble.sh 02c_coassemble_bin.sh 02d_iterative_assembly.sh \
    02e_pool_mqhq_mags.sh 03_filter_contigs.sh 04_map_to_contigs.sh \
    05_metabat.sh 06_checkm2.sh 07_gtdbtk.sh 08_cleanup.sh \
    09_cross_drep.sh 10_bowtie2_map.sh 11_coverm.sh
```

### 3. Generate sample lists

```bash
# Main sample list (all 13 ApoE samples, no blanks or undetermined)
ls /home/projects-phoenix/ApoE_Carnitine/MetaG/raw_reads/*_R1_001.fastq.gz | \
    xargs -n1 basename | sed 's/_R1_001\.fastq\.gz//' \
    > /home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/sample_list.txt

sed -i '/Undetermined_S0/d; /blank_1/d; /NTC/d' sample_list.txt

# Co-assembly group lists (in coassemblies/scripts/)
# Fill in actual sample names
echo -e "sample1\nsample2\nsample3\nsample4\nsample5" \
    > coassemblies/scripts/control_list.txt
echo -e "sample6\nsample7\nsample8\nsample9\nsample10" \
    > coassemblies/scripts/carnitine_list.txt
```

### 4. Initialize the metadata file

```bash
sbatch 00b_init_metadata.sh
```

---

## Running the pipeline

### Phase 1 — Per-sample steps

Run in order. Each step can be re-submitted safely if it fails.

```bash
sbatch 00c_raw_qc.sh
sbatch 01_trim.sh
sbatch 02_assemble.sh
sbatch 03_filter_contigs.sh
sbatch 04_map_to_contigs.sh
sbatch 05_metabat.sh
sbatch 06_checkm2.sh
sbatch 07_gtdbtk.sh
sbatch 08_cleanup.sh

# SingleM on raw reads — no MAG dependency, run any time after 01_trim.sh
# (uses raw reads, not trimmed, per singleM_pipe.py config)
python3 singleM_pipe.py   # reads apoe_sample_list.txt automatically
```

To re-run a single sample:
```bash
sbatch 01_trim.sh con1_c1r
python3 singleM_pipe.py -s con1_c1r
```

### Phase 2 — Co-assembly, iterative binning, and MAG pooling

Run after Phase 1 is complete for all samples.

```bash
sbatch 02b_coassemble.sh          # concatenate group reads + MEGAHIT co-assembly
sbatch 02c_coassemble_bin.sh      # BBMap + MetaBAT + CheckM2 on co-assemblies
sbatch 02d_iterative_assembly.sh  # two rounds iterative assembly + binning
sbatch 02e_pool_mqhq_mags.sh      # pool all 938 MQ/HQ MAGs into ALL_MQHQ_MAGs/
```

### Phase 3 — Cross-sample steps

Run after Phase 2 is complete. This phase should be run on the full ALL_MQHQ_MAGs/
pool — NOT on per-sample MAGs only (see note on pipeline execution order above).

```bash
# Step 1: Run DRAM2 on ALL_MQHQ_MAGs/ (Nextflow-based, separate workflow)
# See Kayla for DRAM2 submission instructions

# Step 2: Build genome_info.csv for dRep (see Notes section below)

# Step 3: Dereplicate pooled MAGs at 99% ANI
sbatch 09_cross_drep.sh

# Step 4: Map all samples to dereplicated MAG catalog
sbatch 10_bowtie2_map.sh

# Step 5: Between-sample abundance tables
sbatch 11_coverm.sh

# Step 6: SingleM appraise + MAG pipe (deferred steps from singleM_pipe.py)
# Re-run singleM_pipe.py after 09_cross_drep.sh completes — it will skip
# already-completed per-sample steps and run Steps 2 and 3 automatically
python3 singleM_pipe.py
```

### Monitoring jobs

```bash
squeue -u $USER
watch -n 5 squeue -u $USER
tail -f slurm_JOBID.out
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
sacct -u $USER --format=JobID,JobName,State,ExitCode,Elapsed --starttime=today
```

---

## Script overview

| Script | Phase | Resources | What it does |
|--------|-------|-----------|--------------|
| `00_config.sh` | — | — | Shared config; sourced by every script. Edit this first. |
| `00b_init_metadata.sh` | — | 1 cpu, 1gb, 5min | Initialize pipeline_metadata.tsv |
| `00c_raw_qc.sh` | 1 | 4 cpu, 50gb, 24hr | Raw read counts + FastQC + MultiQC |
| `01_trim.sh` | 1 | 10 cpu, 50gb, 24hr | sickle + bbduk adapter/poly-G removal |
| `02_assemble.sh` | 1 | 50 cpu, 450gb, 14 days | MEGAHIT per-sample assembly |
| `02b_coassemble.sh` | 2 | 50 cpu, 450gb, 14 days | Concatenate group reads + MEGAHIT co-assembly |
| `02c_coassemble_bin.sh` | 2 | 50 cpu, 450gb, 48hr | BBMap + MetaBAT + CheckM2 on co-assemblies |
| `02d_iterative_assembly.sh` | 2 | 50 cpu, 450gb, 14 days | Two rounds iterative assembly + binning |
| `02e_pool_mqhq_mags.sh` | 2 | 1 cpu, 50gb, 2hr | Pool all 938 MQ/HQ MAGs → ALL_MQHQ_MAGs/ |
| `03_filter_contigs.sh` | 1 | 4 cpu, 50gb, 12hr | pullseq >=2500bp + contig_stats.pl |
| `04_map_to_contigs.sh` | 1 | 50 cpu, 450gb, 14 days | BBMap → BAM → sort (no identity filter) |
| `05_metabat.sh` | 1 | 50 cpu, 200gb, 48hr | MetaBAT binning (default sensitivity) |
| `06_checkm2.sh` | 1 | 10 cpu, 100gb, 48hr | CheckM2 QC + copy MQ/HQ MAGs |
| `07_gtdbtk.sh` | 1 | 20 cpu, 200gb, 14 days | GTDB-Tk taxonomy on per-sample MAGs |
| `08_cleanup.sh` | 1 | 1 cpu, 10gb, 6hr | Delete intermediate mapping files |
| `singleM_pipe.py` | 1+3 | 20 cpu, varies | SingleM pipe (reads, Phase 1) + appraise (Phase 3) |
| `09_cross_drep.sh` | 3 | 50 cpu, 200gb, 48hr | dRep at 99% ANI + rename + concatenate |
| `10_bowtie2_map.sh` | 3 | 50 cpu, 450gb, 14 days | Build Bowtie2 index; map all samples |
| `11_coverm.sh` | 3 | 15 cpu, 100gb, 24hr | Between-sample abundance tables |

---

## MAG sources and counts (ApoE Carnitine pilot)

| Source | Script | MAGs |
|--------|--------|------|
| Per-sample assemblies (13 samples) | 06_checkm2.sh | 693 |
| Control co-assembly | 02c_coassemble_bin.sh | 83 |
| Carnitine co-assembly | 02c_coassemble_bin.sh | 102 |
| Control iterative round 1 | 02d_iterative_assembly.sh | 14 |
| Carnitine iterative round 1 | 02d_iterative_assembly.sh | 22 |
| Control iterative round 2 | 02d_iterative_assembly.sh | 10 |
| Carnitine iterative round 2 | 02d_iterative_assembly.sh | 14 |
| **Total pre-dRep pool** | 02e_pool_mqhq_mags.sh | **938** |

After dRep at 99% ANI, redundant MAGs are collapsed to primary cluster representatives.
The final non-redundant catalog lives in `cross_sample_bowtie_DB/all_samples_derep_MAGs.fa`.

---

## Co-assembly design

Co-assembly was performed for two treatment groups only:
- **Control** (5 samples): reads concatenated and assembled together
- **Carnitine** (5 samples): reads concatenated and assembled together

The antibiotic-only and E. limosum + antibiotic samples were NOT co-assembled because
they have only 3 samples each, which provides less benefit from read pooling.

Co-assembly recovers low-abundance organisms by pooling coverage across samples. A taxon
present at 1x per sample becomes effectively 5x in the co-assembly, bringing it above the
minimum threshold for assembly. The tradeoff is that strain-variable organisms may not
assemble cleanly when reads from multiple samples are mixed.

---

## Iterative assembly design

Iterative assembly follows the WCRC workshop protocol, adapted for co-assembly groups.

**Reference database for unmapped capture:** All MQ/HQ MAGs from per-sample and
co-assembly steps are concatenated into a pre-dRep reference FASTA. BBMap with
semiperfectmode=t maps the full concatenated group reads against this database.
Read pairs that do not recruit to any existing MAG (outu1/outu2) form the unmapped
pool for reassembly. This pre-dRep approach (following the workshop) ensures reads
mapping to any version of an organism are excluded, not just the dereplicated representative.

**Depth BAM:** After MEGAHIT assembles the unmapped reads, the FULL concatenated group
reads are mapped back to the new iterative contigs to generate coverage depth for MetaBAT.
Using the full read set gives MetaBAT a richer depth signal than the unmapped subset alone.
This follows KKA's approach.

**MetaBAT sensitivity:** Default used throughout for internal consistency across all three
assembly strategies. KKA uses --verysensitive; default chosen here so all MAG sources are
binned under identical parameters before pooling.

**Stopping criterion:** Two rounds were performed. Diminishing returns were clear:
- Iter 1: 36 new MQ/HQ MAGs (14 control + 22 carnitine)
- Iter 2: 24 new MQ/HQ MAGs (10 control + 14 carnitine)
The unmapped read pool barely decreased between iterations (~106M → ~104M for control),
indicating remaining reads come from organisms that resist assembly.

---

## SingleM

SingleM provides read-based community profiling independent of assembly and binning.
It uses conserved single-copy marker genes (not 16S) and GTDB taxonomy.

**Why SingleM in addition to CoverM:**
CoverM measures abundance of MAGs you successfully recovered. SingleM measures what was
in the community directly from reads — organisms that failed to assemble, assembled into
unbinned contigs, or failed CheckM2 still appear in SingleM profiles. The appraise step
directly compares these two views to quantify what fraction of the community your MAGs
captured and identify significant missed taxa.

**Per-sample steps (Phase 1, singleM_pipe.py Steps 1, 4, 5) — run on raw reads:**
- Step 1: `singlem pipe` on raw paired-end reads → reads OTU table + profile
- Step 4: `singlem prokaryotic_fraction` → fraction of reads that are prokaryotic
- Step 5: `singlem summarise` → phylum-level relative abundance table + Krona HTML

**Cross-sample steps (Phase 3, singleM_pipe.py Steps 2, 3) — deferred until after dRep:**
- Step 2: `singlem pipe` on cross-sample dereplicated MAG catalog (run ONCE)
- Step 3: `singlem appraise` — compare each sample's reads OTU table against MAG OTU
  table; outputs unrecovered OTUs (community members not in MAG catalog)

The script gracefully skips Steps 2 and 3 if the dRep MAG directory doesn't exist yet,
then runs them automatically when singleM_pipe.py is re-run after 09_cross_drep.sh.

**Output locations:**
- Per-sample: `{SAMPLE}/singlem_output/`
- Cross-sample: `singlem_cross_sample/`

---

## dRep and genome_info.csv

dRep must receive CheckM2 quality scores via --genomeInfo to use completeness and
contamination values from CheckM2 rather than running its own internal CheckM assessment.

Before running 09_cross_drep.sh, build genome_info.csv from all quality_report.tsv files:

```bash
echo "genome,completeness,contamination" > genome_info.csv

# Per-sample MAGs
for SAMPLE in $(cat scripts/apoe_sample_list.txt); do
    TSV="${PROJECT_DIR}/${SAMPLE}/megahit_out/checkm2_v1.1.0_fa/quality_report.tsv"
    [[ -f "${TSV}" ]] && awk -F'\t' -v s="${SAMPLE}" \
        'NR>1{print s"_"$1".fa,"$2","$3}' "${TSV}" >> genome_info.csv
done

# Co-assembly MAGs
for GROUP in control carnitine; do
    TSV="${PROJECT_DIR}/coassemblies/${GROUP}/checkm2_v1.1.0_fa/quality_report.tsv"
    [[ -f "${TSV}" ]] && awk -F'\t' -v g="${GROUP}" \
        'NR>1{print g"_"$1".fa,"$2","$3}' "${TSV}" >> genome_info.csv
done

# Iterative MAGs
for GROUP in control carnitine; do
    for ITER in 1 2; do
        TSV="${PROJECT_DIR}/coassemblies/${GROUP}/iter${ITER}/checkm2_v1.1.0_fa/quality_report.tsv"
        [[ -f "${TSV}" ]] && awk -F'\t' -v g="${GROUP}" -v i="${ITER}" \
            'NR>1{print g"_iter"i"_"$1".fa,"$2","$3}' "${TSV}" >> genome_info.csv
    done
done
```

Pass it to dRep in 09_cross_drep.sh with `--genomeInfo genome_info.csv`.

---

## Skip logic

| Script | Skips if... |
|--------|-------------|
| `00c_raw_qc.sh` | raw_reads in metadata AND FastQC HTMLs exist |
| `01_trim.sh` | Trimmed .fastq files already exist |
| `02_assemble.sh` | final.contigs.fa already exists |
| `02b_coassemble.sh` | final.contigs.fa exists in co-assembly dir |
| `02c_coassemble_bin.sh` | CheckM2 results exist AND MAGs/ has .fa files |
| `02d_iterative_assembly.sh` | Each sub-step checked per group per iteration |
| `02e_pool_mqhq_mags.sh` | ALL_MQHQ_MAGs/ exists with .fa files |
| `03_filter_contigs.sh` | Filtered scaffolds .fa already exists |
| `04_map_to_contigs.sh` | Sorted BAM already exists |
| `05_metabat.sh` | Bins directory exists with .fa files |
| `06_checkm2.sh` | quality_report.tsv exists AND MAGs/ has .fa files |
| `07_gtdbtk.sh` | gtdbtk.*.summary.tsv exists |
| `08_cleanup.sh` | Always safe to re-run |
| `singleM_pipe.py` | Each output file checked per sample per step |
| `09_cross_drep.sh` | dRep, rename, concatenation checked independently |
| `10_bowtie2_map.sh` | Bowtie2 index .bt2 files exist; POSSORT BAM per sample |
| `11_coverm.sh` | Each of three output files checked independently |

---

## Notes

- **Raw read format**: .fastq and .fastq.gz both handled automatically.
- **poly-G trimming**: bbduk in 01_trim.sh trims poly-G tails — essential for NextSeq
  1000/2000 two-color chemistry artifacts.
- **04_map_to_contigs.sh**: No identity filter applied before MetaBAT. MetaBAT receives
  the raw sorted BAM directly. The original reads_to_MAGs_pipeline.py applied a 99%
  reformat.sh filter; this was removed because unfiltered BAMs give MetaBAT better
  coverage signal.
- **MetaBAT sensitivity**: Default used throughout (per-sample, co-assembly, iterative).
  KKA uses --verysensitive; default chosen for internal consistency across all sources.
- **GTDB-Tk**: Run per-sample in 07_gtdbtk.sh. Co-assembly and iterative MAGs receive
  taxonomy through the post-dRep GTDB-Tk run in 09.1_cross_drep_gtdbtk.sh.
- **Separate dereplication**: ApoE_Carnitine and Manure Lagoons must be dereplicated
  separately — they are biologically independent environments. Use apoe_sample_list.txt
  for ApoE-specific cross-sample steps.
- **DRAM2**: Nextflow-based, run on ALL_MQHQ_MAGs/ before 09_cross_drep.sh.
  Contact Kayla for submission instructions.
- **Tools required**: fastqc, multiqc, sickle, bbduk, bbmap, megahit, pullseq, samtools,
  MetaBAT, contig_stats.pl, CheckM2 (conda: checkm2_v1.1.0), GTDB-Tk (conda: gtdbtk_v2.7.0),
  dRep (conda: drep_3.4.2), rename_bins_like_dram.py (/ORG-Data/scripts/), Bowtie2,
  CoverM, SingleM (conda: singlem_0.18.3), DRAM2 (Nextflow).
