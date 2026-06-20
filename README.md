# WGDI Circos Plot

Generate publication-style Circos plots from WGDI synteny results.

This repository provides a Python wrapper around the mature Perl Circos
renderer. The script reads WGDI-format `gff`, WGDI `len`, and
`blockinfo.csv` / `block_information.csv`, maps WGDI gene-order coordinates
back to real chromosome base-pair coordinates, and renders block-level
synteny ribbons between two genomes.

The main plot uses one interval ribbon per WGDI synteny block. Anchor gene
pairs are exported for audit, but they are not drawn as thousands of
point-to-point links by default.

## Features

- Uses real chromosome lengths from WGDI `len` files.
- Preserves chromosome order from the `len` files.
- Maps `blockinfo` gene-order ranges through WGDI `gff` gene coordinates.
- Filters scaffold, contig, unplaced, random, and other small fragments by
  default.
- Renders with Perl Circos and writes all generated Circos input files.
- Produces audit tables for mapped blocks, drawn blocks, all mappable blocks,
  and anchor-level pairs.

## Installation

The script expects Python with `pandas` and a working Perl Circos installation.
A conda environment is the simplest setup:

```bash
conda create -n circos -c bioconda -c conda-forge \
  python pandas circos perl-gd perl-font-ttf perl-math-round \
  perl-config-general perl-clone perl-math-bezier
```

Check the installation:

```bash
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --help
conda run -n circos circos -v
```

By default the script renders through:

```bash
conda run -n circos circos -conf circos.conf
```

This is intentional because directly calling `circos` can fail in some conda
setups when Perl modules are not found.

## Quick Start

```bash
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py \
  --gff1 genomeA.gff \
  --len1 genomeA.len \
  --name1 GenomeA \
  --gff2 genomeB.gff \
  --len2 genomeB.len \
  --name2 GenomeB \
  --block-info genomeA_vs_genomeB_block_information.csv \
  --outdir results/genomeA_vs_genomeB \
  --label GenomeA_vs_GenomeB \
  --output-prefix genomeA_vs_genomeB_circos
```

The command renders PNG and SVG files by default. To only generate Circos
inputs, configuration, and audit tables:

```bash
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py \
  --gff1 genomeA.gff \
  --len1 genomeA.len \
  --name1 GenomeA \
  --gff2 genomeB.gff \
  --len2 genomeB.len \
  --name2 GenomeB \
  --block-info genomeA_vs_genomeB_block_information.csv \
  --outdir results/genomeA_vs_genomeB \
  --no-render
```

## Required Inputs

### WGDI `gff`

This is the WGDI-format gene table, not a standard GFF3 file. The script uses:

| Column | Meaning |
|---:|---|
| 1 | Chromosome or sequence name |
| 2 | Gene ID |
| 3 | Gene start in real bp coordinates |
| 4 | Gene end in real bp coordinates |
| 6 | WGDI gene order index |

### WGDI `len`

The script uses:

| Column | Meaning |
|---:|---|
| 1 | Chromosome or sequence name |
| 2 | Real chromosome length in bp |

Chromosome order in the Circos plot follows the order in each `len` file.

### WGDI `blockinfo.csv`

The block information file must contain:

| Column | Meaning |
|---|---|
| `chr1` | Sequence name in genome A |
| `chr2` | Sequence name in genome B |
| `block1` | Gene-order list or range for genome A |
| `block2` | Gene-order list or range for genome B |

Recommended additional columns include `id`, `start1`, `end1`, `start2`,
`end2`, and `length`. The script requires `blockinfo`; it does not try to
reconstruct blocks from `.collinearity.txt`.

## Important Options

| Option | Default | Description |
|---|---:|---|
| `--min-block-anchors` | `50` | Draw only synteny blocks with at least this many mapped anchors. All mapped blocks are still written to audit files. |
| `--max-links` | `0` | After `--min-block-anchors`, keep only the largest N block ribbons. `0` means no limit. |
| `--exclude-regex` | `scaffold|contig|unplaced|unlocalized|random|chrUn|tig` | Remove likely non-chromosomal fragments. |
| `--keep-regex` | none | If set, retain only sequence names matching this regex before applying `--exclude-regex`. |
| `--chrom-spacing` / `--spacing` | `0.006r` | Gap between adjacent chromosome sectors. Smaller values make the circle more compact. |
| `--bezier-radius` | `0.08r` | Bezier control radius for ribbons. Smaller values reduce broad central ribbon bands. |
| `--link-opacity` | `0.22` | Ribbon opacity. Lower values reduce overplotting. |
| `--no-render` | off | Write files and reports without running Circos. |

Filtering order is fixed: blocks are first filtered by
`--min-block-anchors`; `--max-links` is then applied to the remaining blocks.
Largest blocks are ranked by mapped anchor count, then mapped bp span, then
original order in `blockinfo`.

## Output

Each output directory contains:

| File | Description |
|---|---|
| `*.png`, `*.svg` | Rendered Circos figures, unless `--no-render` is used. |
| `circos.conf` | Main Circos configuration file. |
| `karyotype.txt` | Circos sectors using real chromosome lengths. |
| `links.txt` | Block interval ribbons drawn in the main plot. |
| `links.all_blocks.txt` | All mappable block interval links. |
| `block_links.tsv` | Table of blocks drawn in the main plot. |
| `mapped_blocks.tsv` | Mapped block audit table with WGDI order and bp ranges. |
| `all_block_links.tsv` | All mappable block intervals for audit. |
| `anchor_links.tsv` | Anchor gene pair audit table; not drawn by default. |
| `validation_report.md` | Mapping and rendering summary. |

## Notes on Ribbon Width

The plot uses `ribbon = yes` in Circos, so each synteny block is drawn as a
filled interval-to-interval ribbon. Ribbons can look wider or darker near the
center because Bezier boundaries converge and transparent ribbons overlap.
This is a rendering geometry effect, not a dynamic change in link thickness.

For a cleaner overview, increase `--min-block-anchors`, set `--max-links`, or
lower `--link-opacity`. The default `--bezier-radius 0.08r` is chosen to keep
major block structure visible while reducing broad central bands.

## Batch Use

The script intentionally handles one comparison per command. For multiple
comparisons, use a shell loop or call the script repeatedly with different
input paths and output directories.

## License

No license has been specified yet.
