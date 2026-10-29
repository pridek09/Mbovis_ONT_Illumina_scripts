# Illumina vs ONT pipeline: run order and matching Methods text

## Run order

Put this folder at `/home/pridek/Illumina_ONT/scripts/` and run `bash 06_run_all.sh` (analysis only) or
`bash 06_run_all.sh --full` (downsampling and both vSNP3 pipelines first).

| Step | Script | What it does | Writes |
|---|---|---|---|
| 0 | `00_downsample_to_15x.sh` | Rasusa, 15x, genome 4,349,904, seed 42, per-sample logs, seqkit depth table | `downsampled/` (the folders the vSNP3 wrappers read) and `downsampled/logs/realised_depth.tsv` |
| 1 | `vSNP3_pipeline.sh`, `vSNP3_ONT.sh` | vSNP3 step 1 and step 2 per platform, snp-dists | `analyses/vsnp3_supp/{illumina,ont}_15x/` |
| 1b | `01_record_versions.sh` | Every software version | `software_versions.txt` |
| 2 | `02_make_samples.py` | samples.csv; stops unless exactly 21 paired isolates, one VCF and one `_nodup.bam` each | `samples.csv` |
| 3 | `03_snp_metrics_and_concordance.py` | Tables 1-4, QUAL diagnostics | `concordance/` |
| 4 | `04_stats_and_figures.py` | Wilcoxon, Benjamini-Hochberg, Figures 1, 2, 4, post hoc correlations | `stats/` |
| 5 | `05_alignment_distances.py` | Alignment distances: three masking conditions (i-iii) on the raw vSNP3 output, plus the curated-matrix conditions (iv) when `Final_curated_SNPs_ONT_Illumina_v6.xlsx` is present; Figure 3; snp-dists cross-check | `distances/` |

Retired (do not use for reported numbers): `04_compare_platform_distances.py`, `04_pairwise_snp_distance.py`,
`callable_mask_and_compare.py`, `mask_ppe_and_compare.py`, `vSNP3_combined.sh`. Step 5 replaces them.

## Curated SNP set (Final_curated_SNPs_ONT_Illumina_v6.xlsx)

`05_alignment_distances.py --curated-xlsx ... --samples samples.csv` adds three conditions that use the `Cleaned_ONT` and
`Cleaned_Illumina` sheets instead of the raw vSNP3 tables: `curated` (each platform's own matrix: ONT 833 sites, Illumina 824),
`curated_shared_sites` (824) and `curated_complete_case` (sites unambiguous in all isolates on both platforms, 714).
No PE/PPE or callability masking is applied to them, because the curation already did that and the extra masks would remove
the clusters you validated and kept. Isolates are matched by name (ONT names with two ERR accessions are reduced to the first
accession and checked against `samples.csv`); N, IUPAC codes and `-` are never counted as differences. `--curated-only` runs just
this part without BAMs or FASTA. Steps 03 and 04 (Tables 1-4, Figures 1, 2, 4) still describe the UNCURATED vSNP3 calls.

Expected on the v6 file: 824 shared sites, 9 ONT-only, 0 discordant calls where both platforms are unambiguous, 99.35% identical
genotypes (including ambiguity codes), largest distance difference 9 SNPs, identical distances for all 210 pairs at the 714
complete-case sites, mean distance 0.64 on both platforms in the 67 close pairs.

## Decisions to make before the rerun

1. **ONT QUAL offset.** vSNP3 adds 100 to ONT QUAL. Step 3 prints the minimum QUAL of AC=2 SNP records per platform.
   If the ONT minimum is about 100 above Illumina, rerun step 3 with `--ont-qual-offset 100`. Report whichever you use.
2. **Group split.** `GROUP_SPLIT = 300` in step 4 reflects 169-174 vs 498-656 Illumina SNPs. Check it still separates your groups.
3. **Depth wording.** The realised depth in `realised_depth.tsv` (input) and Table 1 (BAM) must be reported separately.

## Methods text that matches these scripts

Replace the corresponding sections. Items in [brackets] come from the output files named.

**Downsampling.** Illumina and ONT reads from each of the 21 isolates were independently downsampled to a target depth of 15x
(total bases divided by the 4,349,904 bp reference length) with Rasusa v[`rasusa_run_parameters.txt`] using random seed 42 [11].
Isolates with less than 15x raw depth were retained at their full depth ([n] isolates; `realised_depth.tsv`). Depth in the
duplicate-removed alignments (Table 1) was lower than the Rasusa target.

**vSNP3 analysis.** The downsampled datasets were analysed independently with vSNP3 v[`software_versions.txt`] [17], using M. bovis
AF2122/97 (NC_002945.4) as the reference and the Mycobacterium_AF2122 dependency set. Step 1 used BWA-MEM v[ ] and FreeBayes v[ ]
for Illumina reads, and minimap2 v[ ] with bcftools v[ ] (the `--nanopore` option, which raises ONT QUAL values by 100) for ONT reads.
Duplicates were marked and removed with samtools v[ ]. Step 2 was run with the `-a` option on each platform separately.
All other options were left at their defaults.

**SNP definitions.** SNP calls were taken from the vSNP3 step 1 zero-coverage VCFs (`*_zc.vcf`). A high-confidence SNP was a record with
single-base REF and ALT (A, C, G or T), AC=2 and QUAL >= 20 [after subtracting 100 from ONT QUAL / ONT QUAL used as written: choose].
Records with AC=1 were ambiguous, and records with REF=N and ALT='.' were uncalled (zero coverage). At a position where one platform had a
high-confidence SNP, the other platform was classed in this order as shared, uncalled, ambiguous or reference, where reference means
none of the above (it includes positions with a lower-quality or non-SNP record). Quality SNPs (Figure 1f) were single-base records
with QUAL > 150 regardless of AC. Average depth, genome coverage and zero-coverage positions came from `samtools depth -a` on the
duplicate-removed BAM files. Ambiguous SNPs in Figure 1e are AC=1 records (counted per record).

**Allelic and position concordance.** Allelic concordance (Table 2) was the number of shared positions with the identical ALT base divided
by the number of shared positions. Position concordance (Table 3) was shared positions divided by validated positions. Validated positions
were the union of the two platforms' high-confidence SNP positions minus excluded positions. A position was excluded when the platform
without a high-confidence SNP was ambiguous or uncalled there. Concordance without exclusions was shared positions divided by all positions
with a high-confidence SNP on either platform. Values are per isolate and pooled over the 21 isolates (counts summed, then the ratio taken),
and the mean +/- SD (sample SD, n - 1) of the per-isolate percentages is also given.

**Recovery, precision and platform-specific SNPs.** Illumina was the comparator, not a truth set. Recovery (Table 4) was shared SNPs divided
by Illumina high-confidence SNPs, and precision was shared SNPs divided by ONT high-confidence SNPs. Precision describes agreement with
Illumina and does not measure ONT error. SNPs called by only one platform were classed by the other platform's call at that position
(reference, ambiguous or uncalled; Figure 4).

**Pairwise SNP distance.** Distances were calculated from the vSNP3 step 2 SNP alignment of each platform (Mbovis-All group; [863]
variable sites for Illumina and [941] for ONT, from `05` output). The distance between two isolates was the number of alignment positions
at which both isolates carried A, C, G or T and the bases differed, as in snp-dists; the unfiltered matrices were checked against snp-dists
v[ ] and agreed. The root sequence was excluded, leaving 210 pairs. Distances were compared pair by pair as ONT minus Illumina and
summarised as the percentage of pairs that agreed exactly or within 1, 2 or 5 SNPs, the median and interquartile range of the difference,
the Spearman correlation, and, for the [67] pairs with an Illumina distance of 10 or fewer, the mean distance on each platform (Figure 3).
Conditions were (i) the step 2 alignment without additional masking, (ii) after removing sites annotated PE, PPE or PE_PGRS, and
(iii) after PE/PPE masking plus a callability filter. In (iii) a site was removed from both alignments if, on either platform, the median
depth over isolates with reads at the site was below 5, no isolate had reads at the site, more than 25% of the pooled reads had mapping
quality below 20, or the median depth of the isolates carrying the alternative base was below 5. Ambiguity codes alone were not a criterion.
The callability filter and its thresholds were set after inspecting the data and are exploratory.
(iv) Curated set: the comparison was repeated on matrices from which PE/PPE sites, SNPs within 10 bp of another variant that were
absent from Illumina, and ONT-only SNPs explained by repeats, ambiguous calls or paralogs had been removed, and in which manually
reviewed SNPs (IGV, both platforms) had been kept or removed [state final counts: ONT 833, Illumina 824]. No further masking was
applied, and ambiguity codes were excluded pairwise, so they never counted as differences.

**Statistics.** Platform comparisons were paired by isolate (n = 21) with a two-sided Wilcoxon signed-rank test (SciPy v[ ]; isolates with
identical values on both platforms are dropped by the default zero method). The seven Table 1 p-values are reported unadjusted and were
adjusted with the Benjamini-Hochberg procedure over those seven tests [state result from `stats_summary.txt`]. Correlations are
exploratory and unadjusted. Pairwise distance comparisons are descriptive. Panels c-e of Figure 1 use a log scale [state if zeros present].

## Versions confirmed by your first run (`software_versions.txt`)

Rasusa 4.1.0; vSNP3 3.35; BWA 0.7.19; FreeBayes 1.3.9; minimap2 2.30; bcftools 1.23.1; samtools 1.23.1 (vSNP3 environment:
mapping, duplicate removal); snp-dists 1.2.0; seqkit 2.12.0; RAxML: conda package 8.2.13 but the binary prints 8.2.12, so cite
"RAxML v8.2.12 (as reported by the binary)" or check `raxmlHPC -v` yourself.
The analysis scripts (steps 03 and 05) ran in the base environment with **samtools 1.19.2**, so the Methods sentence about
depth, mapping quality and zero-coverage extraction must say v1.19.2, not v1.23.1. Python 3.13.11, SciPy 1.18.1, NumPy 2.4.3,
pandas 3.0.2, openpyxl 3.1.5 and Matplotlib 3.11.2 all match the existing Methods text.

## ONT QUAL offset (resolved by your first run)

ONT QUAL minimum was 120 against 21 for Illumina, i.e. the +100 offset is present. `06_run_all.sh` now passes
`--ont-qual-offset 100`. High-confidence SNP counts do not change (every ONT AC=2 record is already >= 120), but
"quality SNPs" (QUAL > 150) now means raw ONT QUAL > 250, so Table 1 and Figure 1f change. Rerun steps 03-04.
