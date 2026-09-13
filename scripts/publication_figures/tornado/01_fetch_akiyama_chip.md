# External ChIP-seq input: MLL4 / UTX / FLAG-control (Akiyama et al.)

The Figure4 tornado plot overlays public ChIP-seq for MLL4 (=KMT2D) and UTX
(=KDM6A) in human ES cells, from a published reanalysis:

> Akiyama et al., "Functional redundancy between UTY and UTX in regulating
> the localization of transcription factors involved in pluripotency."
> bioRxiv 2025.07.03.663017 / *Development* 2026. PMID 41906541.
> GEO SuperSeries **GSE301298**, ChIP-seq SubSeries **GSE301295**.

This is a heavier, multi-tool alignment pipeline (nf-core/chipseq or
equivalent BWA + deepTools workflow) and is deliberately not duplicated here
-- it is cited and re-run from raw reads, not vendored as code. Align with
default nf-core/chipseq parameters against the same GRCh38 primary assembly
used throughout this repo, single-end, following the paper's stated protocol
(50bp SE, HiSeq 2500).

## Samples needed

| sample_id | target | GSM | SRX | SRR run(s) | ENA fastq URL |
|---|---|---|---|---|---|
| MLL4_hES_rep_1 | MLL4 | GSM9080167 | SRX29493329 | SRR34327783 | `ftp.sra.ebi.ac.uk/vol1/fastq/SRR343/083/SRR34327783/SRR34327783.fastq.gz` |
| Flag_hES_ctrl (input) | Input | GSM9080147 | SRX29493309 | SRR34327805 | `ftp.sra.ebi.ac.uk/vol1/fastq/SRR343/005/SRR34327805/SRR34327805.fastq.gz` |
| UTX_Flag_hES_rep_1 | UTX | GSM9080141 | SRX29493303 | SRR34327812, SRR34327813 | `.../SRR343/012/SRR34327812/...` , `.../SRR343/013/SRR34327813/...` |
| UTX_Flag_hES_rep_2 | UTX | GSM9080142 | SRX29493304 | SRR34327810, SRR34327811 | `.../SRR343/010/SRR34327810/...` , `.../SRR343/011/SRR34327811/...` |
| UTX_Flag_hES_rep_3 | UTX | GSM9080143 | SRX29493305 | SRR34327809 | `.../SRR343/009/SRR34327809/...` |

Multi-run samples: download each SRR separately, concatenate the gzipped
fastqs (valid -- concatenated gzip streams read back as one continuous
`.fastq.gz`).

## What this repo needs from that alignment

Three bigWigs (fold-change-over-control, e.g. `bamCompare`/MACS `bdgcmp -m
FE` against `Flag_hES_ctrl`), placed at:

```
data/external/akiyama_chip/MLL4_hES_rep_1.bigWig
data/external/akiyama_chip/Flag_hES_ctrl.bigWig
data/external/akiyama_chip/UTX_hES_avg.bigWig   # mean of the 3 UTX replicate bigWigs
```

`UTX_hES_avg.bigWig` is the per-bin mean of the three `UTX_Flag_hES_rep_*`
tracks (e.g. `bigwigAverage` from deepTools), not any single replicate.

## Not reproducible from this repo alone

The ChIP alignment itself (adapter trimming, BWA alignment, duplicate
marking, peak/coverage calling) is out of scope for this script chain --
it is a full companion pipeline, cited above by accession rather than
reimplemented. Once the three bigWigs above exist at the paths shown,
`02_define_bound_highconf_canon.sh` onward will run unmodified.
