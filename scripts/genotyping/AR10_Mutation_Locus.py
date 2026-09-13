#!/usr/bin/env python3
"""AR10_Mutation_Locus.py -- read-level quantification, statistics, and audit
visualization of the editing lesion itself (not a downstream DAR/DEG -- the
actual mutation) for either KS1/KS2 clone genotype, selected with --locus.

    conda activate viz_env
    python scripts/AR10_Mutation_Locus.py --locus KMT2D
    python scripts/AR10_Mutation_Locus.py --locus KDM6A
    python scripts/AR10_Mutation_Locus.py --locus KMT2D --igv-report

Both loci's coordinates were derived by hand from each gene's MANE Select
transcript CDS map, then CONFIRMED DIRECTLY against read evidence (not
trusted from the CDS arithmetic alone -- see per-locus notes below, both of
which caught real off-by-something/undercounting problems that pure
coordinate math would have silently missed).

KMT2D: chr12:49,052,381 (GRCh38, + strand), c.1301_1302insA / p.L434fs.
Sits on a homopolymer (chr12:49,052,379-381 = TTT). A 1bp T insertion
anywhere in/adjacent to that run produces an IDENTICAL resulting read
sequence regardless of which exact offset the aligner's CIGAR 'I' op lands
on -- classic indel-in-homopolymer placement ambiguity. A first version of
this script checked one fixed CIGAR offset and undercounted MUTATED ~3x.
A second bug (also found by manual read inspection, reads
LL00123:...1191:48113:11729 and LL00123:...1270:24066:6012): STAR sometimes
represents this exact insertion NOT as an explicit CIGAR 'I' op at all, but
as a silent 1bp register shift -- an all-M block that just happens to
mismatch at the 1-2 positions right after the homopolymer before the shift
becomes invisible again inside the following GGGGGG run. Neither bug is
visible to a classifier that trusts the CIGAR's M/I partition.

KDM6A: chrX:44,961,364 (GRCh38, + strand), c.306_318delinsTAAgaggatcc /
p.N103* -- a 13bp deletion (chrX:44,961,365-377) replaced by an 11bp
insert. This coordinate sits in the middle of exon 3 (44,961,284-392, well
clear of any splice boundary), so the same clipping pattern seen for
KMT2D's reads is caused by the edit itself, not genomic architecture:
divergence from reference triggers STAR to soft-clip immediately at the
edit rather than continue matching (WT reads at the identical anchor
continue matching cleanly for 40+ more bp with zero clipping, confirming
this). Soft-clip windows past the edit are short and often only 1
sequencing-error away from noise, so no fixed suffix requirement is used
here (see LOCI below) -- prefix + a long, distinctive middle sequence
(11-13bp) is unambiguous on its own without needing to also match what
comes after.

FIX APPLIED TO BOTH (the general lesson from the two KMT2D bugs above):
classification never depends on the CIGAR's placement of an insertion, nor
on where an M block "says" a position is. It instead locates the read's
own query index for one fixed, invariant reference anchor position (the
end of a locus-specific PREFIX, chosen to sit just before the variable
region) via a real M/EQ/X alignment, then reads the read's OWN raw sequence
forward from that anchor -- soft-clipped bases included, since they are
still part of the physically sequenced read, just not part of the
aligner's chosen CIGAR path. That raw chunk is compared directly against
two literal candidate strings (PREFIX+SANE_MIDDLE(+SUFFIX) vs
PREFIX+MUTATED_MIDDLE(+SUFFIX)) rather than re-interpreting the CIGAR at
all downstream of the anchor. For KMT2D specifically (homopolymer_char set)
an extra run-length scan additionally distinguishes SANE/MUTATED/
OTHER_INSERTION (run length 3 vs 4 vs >4) rather than just exact-string
matching, since the T-run's length -- not its exact CIGAR placement -- is
the thing that actually varies read-to-read.

Categories (identical semantics for both loci):
  SANE            -- window matches the unedited reference pattern
  MUTATED         -- window matches the edited (mutant-allele) pattern
  OTHER_INSERTION -- (KMT2D only) homopolymer run longer than the mutant
                     pattern -- an insertion bigger than expected
  MISMATCH        -- window covered but matches neither pattern cleanly --
                     sequencing error or a distinct variant nearby
  NOT_COVERING    -- read doesn't reach/fully cover the anchor -- excluded
  DUPLICATE       -- PCR/optical duplicate-flagged (same molecule seen
                     twice, e.g. both mates of a dup-flagged pair) --
                     excluded from the sane/mutated tally as not independent
                     evidence, but kept visible in its own bucket (with what
                     it would have counted as) rather than silently dropped

bam-readcount (KMT2D only -- meaningless for KDM6A's complex delins, which
it has no simple indel notation for) is kept as a REFERENCE point, not a
pass/fail cross-check: it anchors indels at one fixed position the same way
this script's first (buggy) version did, so it is expected to
systematically undercount MUTATED at KMT2D due to the homopolymer -- the
comparison is printed for transparency, not treated as an error.

Quantification ("Sane" / "Sane+Mutated" / "Total") and two normalizations:
  - % of total reads overlapping the site (includes noise categories)
  - % of informative reads only (Sane+Mutated denominator) -- the standard
    variant-calling VAF definition and the recommended one: other categories
    are noise (seq errors, splice edge cases) irrelevant to the sane/mutated
    question, and diluting the fraction by them understates the true
    allelic ratio.

Statistics:
  - PRIMARY/best: two-sided Fisher's exact test on the pooled 2x2 table
    (sane vs mutated) x (WT vs mutant genotype). Fisher's exact is used (not
    a chi-square/proportions z-test) because the WT-mutated cell is
    structurally zero -- exact methods are required, not a normal
    approximation, and Fisher's is exact for 2x2 tables of any cell size.
    NOTE this tests the mutant:sane RATIO within a genotype, NOT a direct
    sane-vs-sane comparison -- see AR10b_Sane_Transcript_DESeq2.R for that.
  - Secondary: two-sided exact binomial test (mutated vs sane+mutated, null
    p=0.5) per mutant-genotype replicate and pooled, testing deviation from
    the naive 50:50 heterozygous/clonal expectation.

RNA-seq only for quantification/stats (ATAC coverage at a single base is
typically 0-3 reads/sample -- too shallow); ATAC still appears in the
per-sample read-evidence audit plot for completeness, so "zero at this
base" is visibly confirmed rather than silently assumed.

Env: viz_env (pysam, matplotlib, scipy; bam-readcount and igv-reports were
added to viz_env for this analysis -- see envs/viz_env.yml).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
from scipy.stats import fisher_exact, binomtest

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_FASTA = REPO_ROOT / "data" / "reference" / "GRCh38.primary_assembly.genome.fa"
ATAC_DIR = REPO_ROOT / "data" / "ATAC_nfcore_output" / "bwa" / "merged_library"
RNA_DIR = REPO_ROOT / "data" / "RNA_nfcore_output" / "star_salmon"
TABLES_DIR = REPO_ROOT / "results" / "tables"

GENOTYPES = ["WT", "KMT2D_Het", "KDM6A_ko"]
REPLICATES = [1, 2, 3]
GENOTYPE_COLOR = {"WT": "#4D4D4D", "KDM6A_ko": "#1B9E77", "KMT2D_Het": "#7570B3"}
BASE_COLOR = {"A": "#4DAF4A", "C": "#377EB8", "G": "#B2662E", "T": "#E41A1C", "N": "#999999"}
CLASS_ORDER = ["SANE", "MUTATED", "OTHER_INSERTION", "MISMATCH", "NOT_COVERING", "DUPLICATE"]

# ---- per-locus configuration -----------------------------------------------
# prefix/sane_middle/mutated_middle/suffix were all confirmed against actual
# read sequence (see module docstring), not derived from CDS arithmetic
# alone -- both loci's naive coordinate math was off by 1-2bp vs what reads
# actually show.
LOCI = {
    "KMT2D": dict(
        chrom="chr12", gene="KMT2D", hgvsc="c.1301_1302insA", hgvsp="p.L434fs",
        mutant_genotype="KMT2D_Het",
        prefix_start_1based=49052373, prefix="TCCTCG",
        sane_middle="TTT", mutated_middle="TTTT", suffix="A",
        homopolymer_char="T",
        anchor_1based=49052381, window_flank=10,
        vcf_pos=49052381, vcf_ref="T", vcf_alt="TT",
        run_bam_readcount_crosscheck=True,
    ),
    "KDM6A": dict(
        chrom="chrX", gene="KDM6A", hgvsc="c.306_318delinsTAAgaggatcc", hgvsp="p.N103*",
        mutant_genotype="KDM6A_ko",
        prefix_start_1based=44961354, prefix="TAGGTCACTTC",
        sane_middle="AACCTCTTATTGG", mutated_middle="TAAGAGGATCC", suffix="",
        homopolymer_char=None,
        anchor_1based=44961364, window_flank=22,
        vcf_pos=44961364, vcf_ref="CAACCTCTTATTGG", vcf_alt="CTAAGAGGATCC",
        run_bam_readcount_crosscheck=False,
    ),
}

LOCUS = None  # populated by configure_locus() at the top of main()


def configure_locus(name):
    global LOCUS, CHROM, GENE, HGVSC, HGVSP, MUTANT_GENOTYPE, PLOT_GENOTYPES
    global PREFIX_START_0BASED, PREFIX, SANE_MIDDLE, MUTATED_MIDDLE, SUFFIX, HOMOPOLYMER_CHAR
    global ANCHOR_1BASED, ANCHOR_0BASED, WINDOW_FLANK, PLOTS_DIR, OUT_PREFIX

    LOCUS = LOCI[name]
    CHROM = LOCUS["chrom"]
    GENE, HGVSC, HGVSP = LOCUS["gene"], LOCUS["hgvsc"], LOCUS["hgvsp"]
    MUTANT_GENOTYPE = LOCUS["mutant_genotype"]
    PLOT_GENOTYPES = ["WT", MUTANT_GENOTYPE]

    PREFIX_START_0BASED = LOCUS["prefix_start_1based"] - 1
    PREFIX, SANE_MIDDLE, MUTATED_MIDDLE, SUFFIX = (
        LOCUS["prefix"], LOCUS["sane_middle"], LOCUS["mutated_middle"], LOCUS["suffix"])
    HOMOPOLYMER_CHAR = LOCUS["homopolymer_char"]

    ANCHOR_1BASED = LOCUS["anchor_1based"]
    ANCHOR_0BASED = ANCHOR_1BASED - 1
    WINDOW_FLANK = LOCUS["window_flank"]

    PLOTS_DIR = REPO_ROOT / "results" / "figures" / "genotype_validation" / f"AR10_{GENE}_Mutation_Locus"
    OUT_PREFIX = f"AR10_{GENE}"


def bam_path(genotype, rep, assay):
    if assay == "RNA":
        return RNA_DIR / f"{genotype}_{rep}_RNA.markdup.sorted.bam"
    return ATAC_DIR / f"{genotype}_{rep}_ATAC_REP1.mLb.clN.sorted.bam"


# BAM CIGAR operation codes (pysam/htslib numeric encoding)
CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_EQ, CIGAR_X = 0, 1, 2, 3, 4, 5, 7, 8
REF_CONSUMING = {CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X}
QUERY_CONSUMING = {CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X}


def classify_reads(bam_file):
    """One entry per overlapping read. Walks read.cigartuples directly
    (never pysam's get_aligned_pairs(), which represents trailing soft-clip
    bases the same way as true insertions -- see module docstring). Builds:
      - window_chars: a display-only dict (ref_0based -> char) for the wider
        +/-WINDOW_FLANK audit pileup plot.
      - the fixed-anchor structural classification: locate the read's own
        query index for PREFIX_START via a real M/EQ/X alignment, then
        compare the read's raw sequence from there against the sane/mutated
        candidate patterns -- see module docstring for why raw sequence,
        not re-parsing the CIGAR downstream of the anchor."""
    win_start = ANCHOR_0BASED - WINDOW_FLANK
    win_end = ANCHOR_0BASED + WINDOW_FLANK + 1
    max_run_check = 8  # extra homopolymer bases scanned for OTHER_INSERTION (KMT2D only)
    raw_chunk_len = len(PREFIX) + max(len(SANE_MIDDLE), len(MUTATED_MIDDLE)) + max_run_check + len(SUFFIX)
    out = []
    for read in bam_file.fetch(CHROM, win_start, win_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.cigartuples is None:
            continue

        ref_pos = read.reference_start
        query_pos = 0
        window_chars = {}
        anchor_query_idx = None       # query index where PREFIX_START_0BASED aligns
        variant_disqualified = False  # a D/N inside/near the window -- can't trust raw readthrough

        for op, length in read.cigartuples:
            if op in (CIGAR_M, CIGAR_EQ, CIGAR_X):
                for offset in range(length):
                    rp = ref_pos + offset
                    base = read.query_sequence[query_pos + offset]
                    if win_start <= rp < win_end:
                        window_chars[rp] = base
                    if rp == PREFIX_START_0BASED:
                        anchor_query_idx = query_pos + offset
            elif op == CIGAR_D:
                for offset in range(length):
                    rp = ref_pos + offset
                    if win_start <= rp < win_end:
                        window_chars[rp] = "-"
                # a deletion anywhere from the anchor through a few bp past the window can
                # break the raw-readthrough assumption below -- disqualify defensively
                if not (ref_pos + length <= PREFIX_START_0BASED or ref_pos > win_end + 5):
                    variant_disqualified = True
            elif op == CIGAR_N:
                if not (ref_pos + length <= PREFIX_START_0BASED or ref_pos > win_end + 5):
                    variant_disqualified = True  # spliced out within/just past the window
            # CIGAR_I, CIGAR_S, CIGAR_H: deliberately NOT specially handled here -- see docstring

            if op in REF_CONSUMING:
                ref_pos += length
            if op in QUERY_CONSUMING:
                query_pos += length

        if anchor_query_idx is None or variant_disqualified:
            underlying_cls, obs, extra = "NOT_COVERING", "", ""
        else:
            raw_chunk = read.query_sequence[anchor_query_idx:anchor_query_idx + raw_chunk_len]
            min_len = len(PREFIX) + min(len(SANE_MIDDLE), len(MUTATED_MIDDLE)) + len(SUFFIX)
            if len(raw_chunk) < min_len:
                underlying_cls, obs, extra = "NOT_COVERING", "", ""  # too close to read's end to resolve
            elif not raw_chunk.startswith(PREFIX):
                underlying_cls, obs, extra = "MISMATCH", raw_chunk[:min_len], ""
            elif HOMOPOLYMER_CHAR is not None:
                # run-length scan (KMT2D-style): distinguishes SANE/MUTATED/OTHER_INSERTION
                # by how many extra homopolymer_char bases follow the prefix, independent of
                # exactly where the aligner's CIGAR would have placed an insertion.
                rest = raw_chunk[len(PREFIX):]
                n_run = 0
                while n_run < len(rest) and rest[n_run] == HOMOPOLYMER_CHAR:
                    n_run += 1
                ref_run_len = len(SANE_MIDDLE)
                obs = PREFIX + rest[:n_run + 1]
                if n_run < len(rest) and (not SUFFIX or rest[n_run] == SUFFIX):
                    if n_run == ref_run_len:
                        underlying_cls, extra = "SANE", ""
                    elif n_run == ref_run_len + 1:
                        underlying_cls, extra = "MUTATED", HOMOPOLYMER_CHAR
                    elif n_run > ref_run_len:
                        underlying_cls, extra = "OTHER_INSERTION", HOMOPOLYMER_CHAR * (n_run - ref_run_len)
                    else:
                        underlying_cls, extra = "MISMATCH", ""
                else:
                    underlying_cls, extra = "MISMATCH", ""
            else:
                # direct literal-pattern comparison (KDM6A-style delins): no run-length
                # concept, no fixed suffix required (see module docstring) -- just does the
                # raw chunk start with prefix+sane_middle, or prefix+mutated_middle?
                rest = raw_chunk[len(PREFIX):]
                if rest.startswith(SANE_MIDDLE) and (not SUFFIX or rest[len(SANE_MIDDLE):].startswith(SUFFIX)):
                    underlying_cls, extra, obs = "SANE", "", PREFIX + SANE_MIDDLE
                elif rest.startswith(MUTATED_MIDDLE) and (not SUFFIX or rest[len(MUTATED_MIDDLE):].startswith(SUFFIX)):
                    underlying_cls, extra, obs = "MUTATED", MUTATED_MIDDLE, PREFIX + MUTATED_MIDDLE
                else:
                    underlying_cls, extra, obs = "MISMATCH", "", raw_chunk[:min_len]

        # PCR/optical duplicates: same original molecule seen twice -- excluded from the
        # sane/mutated tally as not independent evidence, matching bam-readcount's default
        # behavior, but kept fully visible in its own bucket rather than silently dropped.
        cls = "DUPLICATE" if read.is_duplicate else underlying_cls

        out.append({"read": read, "cls": cls, "underlying_cls": underlying_cls,
                     "inserted": extra, "variant_window_obs": obs, "window": window_chars})
    return out


_SITE_FILE = None


def bam_readcount_crosscheck(bam):
    """Independent count via bam-readcount (KMT2D only -- see module
    docstring for why this is skipped for KDM6A's complex delins), purely
    as a REFERENCE point, not a pass/fail cross-check: it anchors indels at
    one fixed position, so it's expected to undercount at a homopolymer."""
    global _SITE_FILE
    if _SITE_FILE is None:
        import tempfile
        fh = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        fh.write(f"{CHROM}\t{ANCHOR_1BASED}\t{ANCHOR_1BASED}\n")
        fh.close()
        _SITE_FILE = fh.name
    result = subprocess.run(
        ["bam-readcount", "-w1", "-i", "-f", str(REF_FASTA),
         "-l", _SITE_FILE, str(bam)],
        capture_output=True, text=True,
    )
    line = result.stdout.strip()
    if not line:
        return {"sane": 0, "mutated": 0}
    fields = line.split("\t")
    sane, mutated = 0, 0
    for allele_field in fields[4:]:
        allele, count = allele_field.split(":")[0], int(allele_field.split(":")[1])
        if allele == HOMOPOLYMER_CHAR:
            sane = count
        elif allele == f"+{HOMOPOLYMER_CHAR}":
            mutated += count
    return {"sane": sane, "mutated": mutated}


def quantify_all():
    rows, evidence = [], {}
    for genotype in GENOTYPES:
        for rep in REPLICATES:
            for assay in ("RNA", "ATAC"):
                bam = bam_path(genotype, rep, assay)
                sample = f"{genotype}_{rep}"
                with pysam.AlignmentFile(str(bam), "rb") as bf:
                    reads = classify_reads(bf)
                counts = {c: sum(1 for r in reads if r["cls"] == c) for c in CLASS_ORDER}
                sane, mutated = counts["SANE"], counts["MUTATED"]
                informative = sane + mutated
                # "Total" excludes PCR/optical duplicates (not independent evidence -- see
                # classify_reads); full duplicate-inclusive count kept as a separate audit column.
                total = sum(v for c, v in counts.items() if c != "DUPLICATE")
                total_incl_dup = total + counts["DUPLICATE"]

                if LOCUS["run_bam_readcount_crosscheck"]:
                    xcheck = bam_readcount_crosscheck(bam)
                    if mutated > 0 or xcheck["mutated"] > 0:
                        print(f"[AR10:{GENE}] {sample} {assay}: fixed-anchor sane={sane} mutated={mutated}  "
                              f"vs  bam-readcount (single-position-anchored) sane={xcheck['sane']} mutated={xcheck['mutated']}",
                              file=sys.stderr)
                else:
                    xcheck = {"sane": None, "mutated": None}

                rows.append({
                    "sample": sample, "genotype": genotype, "assay": assay,
                    "sane": sane, "mutated": mutated, "other_insertion": counts["OTHER_INSERTION"],
                    "mismatch": counts["MISMATCH"], "not_covering": counts["NOT_COVERING"],
                    "duplicate_excluded": counts["DUPLICATE"],
                    "sane_plus_mutated": informative, "total_reads_overlapping": total,
                    "total_reads_incl_duplicates": total_incl_dup,
                    "pct_mutated_of_total": round(100 * mutated / total, 1) if total else None,
                    "pct_sane_of_total": round(100 * sane / total, 1) if total else None,
                    "pct_mutated_of_informative": round(100 * mutated / informative, 1) if informative else None,
                    "pct_sane_of_informative": round(100 * sane / informative, 1) if informative else None,
                    "bam_readcount_sane": xcheck["sane"], "bam_readcount_mutated": xcheck["mutated"],
                })
                evidence[(genotype, rep, assay)] = reads
    return rows, evidence


def run_stats(rows):
    mut_rna = [r for r in rows if r["genotype"] == MUTANT_GENOTYPE and r["assay"] == "RNA"]
    wt_rna = [r for r in rows if r["genotype"] == "WT" and r["assay"] == "RNA"]

    mut_sane = sum(r["sane"] for r in mut_rna)
    mut_mut = sum(r["mutated"] for r in mut_rna)
    wt_sane = sum(r["sane"] for r in wt_rna)
    wt_mut = sum(r["mutated"] for r in wt_rna)

    stats_rows = []

    # PRIMARY: pooled Fisher's exact test, sane vs mutated, WT vs mutant genotype
    odds_ratio, p_fisher = fisher_exact([[wt_sane, wt_mut], [mut_sane, mut_mut]])
    stats_rows.append({
        "test": "Fisher exact (pooled, primary)",
        "comparison": f"WT vs {MUTANT_GENOTYPE}, sane vs mutated read counts",
        "n_wt_sane": wt_sane, "n_wt_mutated": wt_mut, "n_mut_sane": mut_sane, "n_mut_mutated": mut_mut,
        "odds_ratio": round(odds_ratio, 4) if odds_ratio not in (float("inf"),) else "inf",
        "p_value": p_fisher,
        "interpretation": f"sane-allele fraction significantly reduced in {MUTANT_GENOTYPE} vs WT" if p_fisher < 0.05 else "not significant",
    })

    # SECONDARY: per-replicate + pooled one-sample exact binomial test vs 0.5
    for r in mut_rna:
        n = r["sane_plus_mutated"]
        if n == 0:
            continue
        res = binomtest(r["mutated"], n, p=0.5, alternative="two-sided")
        stats_rows.append({
            "test": "Binomial exact (one-sample, secondary)", "comparison": f"{r['sample']} vs null p=0.5",
            "n_wt_sane": "", "n_wt_mutated": "", "n_mut_sane": r["sane"], "n_mut_mutated": r["mutated"],
            "odds_ratio": "", "p_value": res.pvalue,
            "interpretation": "mutated allele significantly below 50% (NMD-consistent)" if res.pvalue < 0.05 else "not significant",
        })
    res_pooled = binomtest(mut_mut, mut_sane + mut_mut, p=0.5, alternative="two-sided")
    stats_rows.append({
        "test": "Binomial exact (one-sample, pooled, secondary)", "comparison": f"{MUTANT_GENOTYPE} pooled vs null p=0.5",
        "n_wt_sane": "", "n_wt_mutated": "", "n_mut_sane": mut_sane, "n_mut_mutated": mut_mut,
        "odds_ratio": "", "p_value": res_pooled.pvalue,
        "interpretation": "mutated allele significantly below 50% (NMD-consistent)" if res_pooled.pvalue < 0.05 else "not significant",
    })
    return stats_rows


def write_csv(rows, path, fieldnames):
    import csv
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[AR10:{GENE}] Wrote {path} ({len(rows)} rows)")


def sig_stars(p):
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_quantification(rows, stats_rows):
    rna_rows = [r for r in rows if r["assay"] == "RNA" and r["genotype"] in PLOT_GENOTYPES]
    order = [f"{g}_{r}" for g in GENOTYPES for r in REPLICATES]
    rna_rows.sort(key=lambda r: order.index(r["sample"]))

    binom_p = {s["comparison"].split(" vs")[0]: s["p_value"]
               for s in stats_rows if s["test"].startswith("Binomial exact (one-sample, secondary)")}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, key, title in zip(
        axes,
        ("pct_mutated_of_total", "pct_mutated_of_informative"),
        ("% of total reads overlapping site\n(includes noise categories)",
         "% of informative (sane+mutated) reads\n(standard VAF definition -- recommended)"),
    ):
        xs = range(len(rna_rows))
        vals = [r[key] or 0 for r in rna_rows]
        colors = [GENOTYPE_COLOR[r["genotype"]] for r in rna_rows]
        ax.bar(xs, vals, color=colors)
        for x, r in zip(xs, rna_rows):
            ax.text(x, (r[key] or 0) + 0.8, f"{r[key]}%\n(n={r['total_reads_overlapping'] if key.endswith('total') else r['sane_plus_mutated']})",
                    ha="center", va="bottom", fontsize=7)
            if key == "pct_mutated_of_informative" and r["sample"] in binom_p:
                ax.text(x, (r[key] or 0) + 6.5, sig_stars(binom_p[r["sample"]]), ha="center", va="bottom",
                        fontsize=9, color=GENOTYPE_COLOR[MUTANT_GENOTYPE], fontweight="bold")
        ax.axhline(50, color="black", lw=0.8, ls="--", alpha=0.4)
        ax.set_xticks(xs)
        ax.set_xticklabels([r["sample"] for r in rna_rows], rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 60)
        ax.set_title(title, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].text(0.98, 0.95, "* vs null p=0.5 (exact binomial)\n**** p<0.0001", transform=axes[1].transAxes,
                 ha="right", va="top", fontsize=6.5, color="#555555")
    axes[0].set_ylabel("Mutated allele fraction (%)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=g) for g, c in GENOTYPE_COLOR.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.04), fontsize=8)
    fig.suptitle(f"{GENE} {HGVSC} ({HGVSP}) -- {CHROM}:{ANCHOR_1BASED:,} (GRCh38) -- RNA-seq allele quantification",
                 fontsize=10, y=1.1)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = PLOTS_DIR / f"{OUT_PREFIX}_mutation_allele_quant.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[AR10:{GENE}] Wrote {out}")
    plt.close(fig)


def load_size_factors():
    """DESeq2 size factors from the R02 fit, exported by
    AR10b_Sane_Transcript_DESeq2.R. Used to normalize the raw locus-anchor
    read counts onto the same depth-corrected scale as every other
    count/LFC in this project, rather than plotting raw counts that
    conflate real signal with per-sample sequencing-depth differences."""
    path = TABLES_DIR / "R02_size_factors.csv"
    if not path.exists():
        return None
    import csv as _csv
    with open(path) as fh:
        return {row["sample_rna"]: float(row["size_factor"]) for row in _csv.DictReader(fh)}


def plot_counts_breakdown(rows, stats_rows, include_other=True):
    """The literal Sane / Sane+Mutated / Total breakdown: a stacked bar per
    RNA sample (Sane at bottom, Mutated on top) plus, if include_other=True,
    a translucent Total-including-noise outline on top of that -- DESeq2-
    size-factor-normalized so bar heights are comparable across samples
    despite differing sequencing depth, with the Sane-WT-vs-Sane-mutant
    DESeq2 result annotated as a bracket (NOT Fisher's exact -- see module
    docstring for the distinction). include_other=False drops the hatched
    mismatch/not-covering extension entirely, showing only the informative
    Sane+Mutated stack -- a separate, differently-named output, not a
    replacement for the include_other=True version."""
    rna_rows = [r for r in rows if r["assay"] == "RNA" and r["genotype"] in PLOT_GENOTYPES]
    order = [f"{g}_{r}" for g in GENOTYPES for r in REPLICATES]
    rna_rows.sort(key=lambda r: order.index(r["sample"]))

    size_factors = load_size_factors()
    norm = (lambda r, v: v / size_factors[f"{r['sample']}_RNA"]) if size_factors else (lambda r, v: v)
    unit = "DESeq2-normalized read count" if size_factors else "Read count (raw -- run AR10b for normalization)"

    fig, ax = plt.subplots(figsize=(9, 5))
    xs = list(range(len(rna_rows)))
    sane_vals = [norm(r, r["sane"]) for r in rna_rows]
    mut_vals = [norm(r, r["mutated"]) for r in rna_rows]

    ax.bar(xs, sane_vals, color="#4D4D4D", label="Sane (ref allele)")
    ax.bar(xs, mut_vals, bottom=sane_vals, color="#E41A1C", label="Mutated (edited allele)")
    bottoms = [s + m for s, m in zip(sane_vals, mut_vals)]

    if include_other:
        other_vals = [norm(r, r["other_insertion"] + r["mismatch"] + r["not_covering"]) for r in rna_rows]
        ax.bar(xs, other_vals, bottom=bottoms, color="#CCCCCC", alpha=0.5, hatch="//",
               edgecolor="#999999", label="Other/excluded (mismatch, not-covering)\n= Total")
        top_vals = [norm(r, r["total_reads_overlapping"]) for r in rna_rows]
        top_label = "Total"
    else:
        top_vals = bottoms
        top_label = "Sane+Mutated"

    for x, r, sv, mv in zip(xs, rna_rows, sane_vals, mut_vals):
        ax.text(x, sv / 2, f"{sv:.1f}", ha="center", va="center", fontsize=7.5, color="white")
        if r["mutated"]:
            ax.text(x, sv + mv / 2, f"{mv:.1f}", ha="center", va="center",
                     fontsize=7.5, color="white", fontweight="bold")
    for x, tv in zip(xs, top_vals):
        ax.text(x, tv + 1.5, f"{top_label}={tv:.1f}", ha="center", va="bottom", fontsize=6.5, color="#555555")

    mut_idx = [i for i, r in enumerate(rna_rows) if r["genotype"] == MUTANT_GENOTYPE]
    wt_idx = [i for i, r in enumerate(rna_rows) if r["genotype"] == "WT"]
    y_bracket = max(top_vals) + 10
    x0, x1 = (min(wt_idx) + max(wt_idx)) / 2, (min(mut_idx) + max(mut_idx)) / 2
    ax.plot([x0, x0, x1, x1], [y_bracket - 1, y_bracket, y_bracket, y_bracket - 1], color="black", lw=1)
    deseq2_path = TABLES_DIR / f"{OUT_PREFIX}_sane_transcript_DESeq2.csv"
    if deseq2_path.exists():
        import csv as _csv
        with open(deseq2_path) as fh:
            d = next(_csv.DictReader(fh))
        lfc, padj = float(d["log2FoldChange"]), float(d["padj"])
        label = f"Sane WT vs Sane {MUTANT_GENOTYPE} (DESeq2): log2FC={lfc:.2f}, padj={padj:.2g} ({sig_stars(padj)})"
    else:
        fisher_p = next(s["p_value"] for s in stats_rows if s["test"].startswith("Fisher exact"))
        label = (f"Fisher's exact (mutant:sane ratio, NOT sane-vs-sane) p={fisher_p:.2g} "
                 f"({sig_stars(fisher_p)}) -- run AR10b for the sane-vs-sane test")
    ax.text((x0 + x1) / 2, y_bracket + 0.5, label, ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_xticks(xs)
    ax.set_xticklabels([r["sample"] for r in rna_rows], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(unit)
    ax.set_ylim(0, y_bracket + 6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)
    norm_note = "DESeq2 size-factor normalized" if size_factors else "raw counts, NOT normalized"
    scope_note = "Sane + Mutated only (no other/excluded reads shown)" if not include_other else "Sane / Mutated / Total"
    ax.set_title(f"{GENE} {HGVSC} ({HGVSP}) -- {CHROM}:{ANCHOR_1BASED:,} (GRCh38)\n"
                 f"{scope_note} read counts, RNA-seq -- {norm_note}",
                 fontsize=9.5)
    fig.tight_layout()
    suffix = "" if include_other else "_sane_mutated_only"
    for ext in ("png", "pdf"):
        out = PLOTS_DIR / f"{OUT_PREFIX}_read_count_breakdown{suffix}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[AR10:{GENE}] Wrote {out}")
    plt.close(fig)


def load_sane_deseq2_caption():
    """One-line caption of the Sane-WT-vs-Sane-mutant DESeq2 Wald test (see
    AR10b_Sane_Transcript_DESeq2.R), for display on the pileup pages --
    distinct from the Fisher's exact test, which tests the mutant:sane RATIO
    within a genotype, not a direct sane-vs-sane comparison."""
    path = TABLES_DIR / f"{OUT_PREFIX}_sane_transcript_DESeq2.csv"
    if not path.exists():
        return f"Sane WT vs Sane {MUTANT_GENOTYPE} (DESeq2): not yet run -- run AR10b_Sane_Transcript_DESeq2.R --locus {GENE}"
    import csv as _csv
    with open(path) as fh:
        d = next(_csv.DictReader(fh))
    lfc, padj = float(d["log2FoldChange"]), float(d["padj"])
    return f"Sane WT vs Sane {MUTANT_GENOTYPE} (DESeq2 Wald test, size-factor normalized): log2FC={lfc:.2f}, padj={padj:.2g} ({sig_stars(padj)})"


def plot_read_evidence(evidence, genotypes=None, show_duplicates=False):
    """One page per sample x assay: every classified read rendered as a
    monospace, base-colored row, grouped by classification bucket with a
    counted header, so the sane/mutated tally can be visually audited read
    by read rather than trusted as a black-box number. `genotypes` restricts
    which genotypes get a page (default: WT + the locus's mutant genotype --
    the third genotype is a different gene's edit, irrelevant here).
    `show_duplicates` (default False) omits the DUPLICATE bucket from the
    rendered page entirely -- those reads are already excluded from the
    sane/mutated tally either way; this just keeps the page itself clean."""
    genotypes = genotypes if genotypes is not None else PLOT_GENOTYPES
    class_order = [c for c in CLASS_ORDER if show_duplicates or c != "DUPLICATE"]
    caption = load_sane_deseq2_caption()
    win_start = ANCHOR_0BASED - WINDOW_FLANK
    win_end = ANCHOR_0BASED + WINDOW_FLANK + 1
    ref = pysam.FastaFile(str(REF_FASTA))
    ref_seq = ref.fetch(CHROM, win_start, win_end)
    ref.close()
    anchor_col = ANCHOR_0BASED - win_start  # 0-based column index of the anchor base

    out_pdf = PLOTS_DIR / f"{OUT_PREFIX}_read_evidence_pileup.pdf"
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf) as pdf:
        for genotype in genotypes:
            for rep in REPLICATES:
                for assay in ("RNA", "ATAC"):
                    reads = evidence[(genotype, rep, assay)]
                    title = (f"{genotype}_{rep} {assay} -- {CHROM}:{win_start+1:,}-{win_end:,} -- "
                             f"anchor {CHROM}:{ANCHOR_1BASED:,} (edit highlighted)\n{caption}")
                    _render_pileup_page(pdf, reads, ref_seq, win_start, anchor_col, class_order, title)
    print(f"[AR10:{GENE}] Wrote {out_pdf} (one page per sample x assay, {len(genotypes)*len(REPLICATES)*2} pages)")


def plot_read_evidence_by_genotype(evidence, genotypes=None, assay="RNA", show_duplicates=False):
    """Genotype-level companion to plot_read_evidence: pools reads across all
    3 replicates of each genotype into ONE page per genotype -- read-level
    evidence for the same 3-way comparison as the sane_transcript_by_genotype
    bar chart, not just a per-replicate view."""
    genotypes = genotypes if genotypes is not None else GENOTYPES
    class_order = [c for c in CLASS_ORDER if show_duplicates or c != "DUPLICATE"]
    caption = load_sane_deseq2_caption()
    win_start = ANCHOR_0BASED - WINDOW_FLANK
    win_end = ANCHOR_0BASED + WINDOW_FLANK + 1
    ref = pysam.FastaFile(str(REF_FASTA))
    ref_seq = ref.fetch(CHROM, win_start, win_end)
    ref.close()
    anchor_col = ANCHOR_0BASED - win_start

    out_pdf = PLOTS_DIR / f"{OUT_PREFIX}_read_evidence_pileup_by_genotype.pdf"
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf) as pdf:
        for genotype in genotypes:
            reads = [r for rep in REPLICATES for r in evidence[(genotype, rep, assay)]]
            title = (f"{genotype} {assay} (pooled across {len(REPLICATES)} replicates) -- "
                     f"{CHROM}:{win_start+1:,}-{win_end:,} -- anchor {CHROM}:{ANCHOR_1BASED:,} "
                     f"(edit highlighted)\n{caption}")
            _render_pileup_page(pdf, reads, ref_seq, win_start, anchor_col, class_order, title)
    print(f"[AR10:{GENE}] Wrote {out_pdf} (one page per genotype, {len(genotypes)} pages)")


def _render_pileup_page(pdf, reads, ref_seq, win_start, anchor_col, class_order, title):
    """Shared page renderer used by both plot_read_evidence (per-sample) and
    plot_read_evidence_by_genotype (pooled-per-genotype): every read in
    `reads`, grouped by classification bucket, one monospace color-coded row
    each, the edit highlighted starting at the anchor column."""
    def row_string(entry):
        return [entry["window"].get(win_start + i, " ") for i in range(len(ref_seq))]

    groups = {c: [r for r in reads if r["cls"] == c] for c in class_order}
    n_rows = sum(len(v) + (1 if v else 0) for v in groups.values()) + 3
    fig_h = max(2.2, 0.22 * n_rows + 1.2)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.set_xlim(-0.5, len(ref_seq) + 8.5)
    ax.set_ylim(0, n_rows)
    ax.axis("off")

    y = n_rows - 1
    for i, base in enumerate(ref_seq):
        ax.text(i, y, base, family="monospace", fontsize=9, ha="center", va="center",
                color=BASE_COLOR.get(base, "black"), fontweight="bold")
    ax.text(len(ref_seq) + 1, y, "reference", family="monospace", fontsize=8, va="center", color="#555555")
    ax.axvspan(anchor_col + 0.5, anchor_col + 1.5, ymin=0, ymax=1, color="#7570B3", alpha=0.08, lw=0)
    y -= 1.3

    for cls in class_order:
        grp = groups[cls]
        if not grp:
            continue
        ax.text(-0.5, y, f"{cls} (n={len(grp)})", family="monospace", fontsize=8.5,
                fontweight="bold", va="center", ha="left", color="#222222")
        y -= 1
        for entry in grp:
            chars = row_string(entry)
            for i, ch in enumerate(chars):
                color = BASE_COLOR.get(ch, "#BBBBBB") if ch not in (" ", "-") else "#CCCCCC"
                ax.text(i, y, ch if ch != " " else "·", family="monospace", fontsize=7.5,
                        ha="center", va="center", color=color)
            u_cls = entry["underlying_cls"]
            if u_cls == "MUTATED":
                box_color = "#E41A1C" if cls != "DUPLICATE" else "#F4A6A6"
                ax.text(anchor_col + 1, y, entry["inserted"][:4], family="monospace", fontsize=7.5,
                        ha="center", va="center", color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor=box_color, edgecolor="none"))
            elif u_cls == "OTHER_INSERTION":
                box_color = "#FF7F00" if cls != "DUPLICATE" else "#FFC98C"
                ax.text(anchor_col + 1, y, entry["inserted"][:4] or "?", family="monospace", fontsize=7.5,
                        ha="center", va="center", color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor=box_color, edgecolor="none"))
            read = entry["read"]
            tag = f" [dup; would be {u_cls}]" if cls == "DUPLICATE" else ""
            obs_note = f" win={entry['variant_window_obs']}" if entry["variant_window_obs"] else ""
            ax.text(len(ref_seq) + 1, y, f"{read.query_name[-14:]}{tag}{obs_note}", family="monospace",
                    fontsize=6, va="center", color="#888888")
            y -= 1
        y -= 0.3

    ax.set_title(title, fontsize=8.5, loc="left")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def build_igv_report():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    vcf_path = PLOTS_DIR / f"_{GENE}_variant.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n'
        '##INFO=<ID=HGVSC,Number=1,Type=String,Description="cDNA HGVS notation">\n'
        '##INFO=<ID=HGVSP,Number=1,Type=String,Description="Protein HGVS notation">\n'
        f"##contig=<ID={CHROM},length=133275309>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"{CHROM}\t{LOCUS['vcf_pos']}\t{GENE}_{HGVSP}\t{LOCUS['vcf_ref']}\t{LOCUS['vcf_alt']}\t.\t.\t"
        f"GENE={GENE};HGVSC={HGVSC};HGVSP={HGVSP}\n"
    )
    tracks = []
    for g in GENOTYPES:
        for rep in REPLICATES:
            tracks.append(str(bam_path(g, rep, "RNA")))
    for g in GENOTYPES:
        for rep in REPLICATES:
            tracks.append(str(bam_path(g, rep, "ATAC")))
    out_html = PLOTS_DIR / f"{OUT_PREFIX}_{HGVSP}_igv_report.html"
    subprocess.run(
        ["create_report", str(vcf_path), "--fasta", str(REF_FASTA),
         "--tracks", *tracks, "--info-columns", "GENE", "HGVSC", "HGVSP",
         "--flanking", "60", "--sort", "BASE", "--standalone",
         "--title", f"{GENE} {HGVSC} / {HGVSP} -- {CHROM}:{ANCHOR_1BASED:,} (GRCh38)",
         "--output", str(out_html)],
        check=True,
    )
    vcf_path.unlink()
    print(f"[AR10:{GENE}] Wrote {out_html}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", required=True, choices=sorted(LOCI), help="which mutation locus to analyze")
    ap.add_argument("--igv-report", action="store_true", help="also (re)build the standalone igv.js HTML report")
    args = ap.parse_args()

    configure_locus(args.locus)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[AR10:{GENE}] {GENE} {HGVSC} / {HGVSP} -- anchor {CHROM}:{ANCHOR_1BASED:,} (GRCh38), "
          f"sane='{PREFIX}{SANE_MIDDLE}{SUFFIX}' vs mutated='{PREFIX}{MUTATED_MIDDLE}{SUFFIX}'")

    rows, evidence = quantify_all()
    write_csv(rows, TABLES_DIR / f"{OUT_PREFIX}_mutation_allele_quant.csv", list(rows[0].keys()))

    stats_rows = run_stats(rows)
    write_csv(stats_rows, TABLES_DIR / f"{OUT_PREFIX}_mutation_allele_stats.csv", list(stats_rows[0].keys()))
    for s in stats_rows:
        print(f"[AR10:{GENE}] {s['test']:45s} {s['comparison']:45s} p={s['p_value']:.3g}  {s['interpretation']}")

    plot_quantification(rows, stats_rows)
    plot_counts_breakdown(rows, stats_rows)
    plot_counts_breakdown(rows, stats_rows, include_other=False)
    plot_read_evidence(evidence)
    plot_read_evidence_by_genotype(evidence)

    if args.igv_report:
        build_igv_report()


if __name__ == "__main__":
    main()
