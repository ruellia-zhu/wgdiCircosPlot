#!/usr/bin/env python3
"""Convert WGDI block_information.csv results to Perl Circos inputs.

The plotting backend is the mature Perl Circos package.  This script only
parses WGDI files, maps gene-order coordinates back to real chromosome bp
coordinates, writes Circos input/config files, and optionally runs Circos.
"""

from __future__ import annotations

import argparse
import bisect
import colorsys
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - exercised only on bad envs
    raise SystemExit(
        "Missing Python dependency: pandas. Run this script in the circos environment, "
        "or install dependencies with `conda install -n circos -c conda-forge pandas`."
    ) from exc


DEFAULT_PALETTE = [
    (122, 127, 0),
    (192, 0, 0),
    (255, 30, 30),
    (255, 0, 184),
    (255, 208, 0),
    (255, 152, 0),
    (247, 199, 199),
    (166, 124, 0),
    (108, 118, 0),
    (70, 130, 180),
    (0, 150, 150),
    (120, 80, 180),
    (80, 170, 90),
    (210, 90, 40),
]

DEFAULT_EXCLUDE_RE = r"scaffold|contig|unplaced|unlocalized|random|chrUn|tig"


@dataclass(frozen=True)
class Gene:
    chrom: str
    gene_id: str
    start: int
    end: int
    order: int


@dataclass(frozen=True)
class LensData:
    path: Path
    lengths: Dict[str, int]
    gene_counts: Dict[str, Optional[int]]
    order: List[str]
    skipped: List[str]


@dataclass
class LinkRecord:
    block_id: str
    block_rank: int
    source: str
    chr1: str
    start1: int
    end1: int
    gene_order1: Optional[int]
    gene_id1: str
    chr2: str
    start2: int
    end2: int
    gene_order2: Optional[int]
    gene_id2: str
    color: str


@dataclass
class BlockLink:
    block_id: str
    block_rank: int
    source: str
    chr1: str
    start1: int
    end1: int
    chr2: str
    start2: int
    end2: int
    original_anchor_count: int
    mapped_anchor_count: int
    color: str


@dataclass
class BlockRecord:
    block_id: str
    source: str
    chr1: str
    chr2: str
    original_start1: Optional[int]
    original_end1: Optional[int]
    original_start2: Optional[int]
    original_end2: Optional[int]
    original_anchor_count: int
    mapped_anchor_count: int
    failed_anchor_count: int
    mapped_start1: Optional[int]
    mapped_end1: Optional[int]
    mapped_start2: Optional[int]
    mapped_end2: Optional[int]
    pvalue: Optional[float] = None
    ks_median: Optional[float] = None
    note: str = ""


@dataclass(frozen=True)
class DatasetSpec:
    label: str
    gff1: Path
    len1: Path
    name1: str
    gff2: Path
    len2: Path
    name2: str
    block_info: Path
    output_prefix: str


def safe_id(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "x"


def sector_id(prefix: str, chrom: str) -> str:
    return f"{safe_id(prefix)}_{safe_id(chrom)}"


def parse_rgb_line(line: str) -> Optional[Tuple[int, int, int]]:
    s = line.strip()
    if not s or s.startswith("#") and not re.fullmatch(r"#[0-9A-Fa-f]{6}", s):
        return None
    if "=" in s:
        s = s.split("=", 1)[1].strip()
    if s.startswith("#") and re.fullmatch(r"#[0-9A-Fa-f]{6}", s):
        return tuple(int(s[i : i + 2], 16) for i in (1, 3, 5))
    parts = re.split(r"[,\s]+", re.sub(r"\s+#.*$", "", s))
    if len(parts) >= 3:
        try:
            rgb = tuple(int(float(x)) for x in parts[:3])
        except ValueError:
            return None
        if all(0 <= x <= 255 for x in rgb):
            return rgb
    return None


def generated_palette(n: int) -> List[Tuple[int, int, int]]:
    colors = list(DEFAULT_PALETTE)
    h = 0.08
    while len(colors) < n:
        h = (h + 0.61803398875) % 1.0
        r, g, b = colorsys.hsv_to_rgb(h, 0.70, 0.82)
        colors.append((int(r * 255), int(g * 255), int(b * 255)))
    return colors[:n]


def load_palette(path: Optional[Path], n: int) -> List[Tuple[int, int, int]]:
    colors: List[Tuple[int, int, int]] = []
    if path:
        with path.open() as handle:
            for line in handle:
                rgb = parse_rgb_line(line)
                if rgb is not None:
                    colors.append(rgb)
    if len(colors) < n:
        colors.extend(generated_palette(n)[len(colors) :])
    return colors[:n]


def read_lens(path: Path, exclude_re: Optional[re.Pattern], keep_re: Optional[re.Pattern]) -> LensData:
    lengths: Dict[str, int] = {}
    gene_counts: Dict[str, Optional[int]] = {}
    order: List[str] = []
    skipped: List[str] = []

    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                skipped.append(f"line {line_no}: expected at least 2 columns")
                continue
            chrom = parts[0]
            if keep_re and not keep_re.search(chrom):
                skipped.append(chrom)
                continue
            if exclude_re and exclude_re.search(chrom):
                skipped.append(chrom)
                continue
            try:
                length = int(float(parts[1]))
            except ValueError:
                skipped.append(f"{chrom}: invalid length {parts[1]}")
                continue
            if length <= 0:
                skipped.append(f"{chrom}: non-positive length {length}")
                continue
            count = None
            if len(parts) >= 3:
                try:
                    count = int(float(parts[2]))
                except ValueError:
                    count = None
            lengths[chrom] = length
            gene_counts[chrom] = count
            order.append(chrom)

    if not order:
        raise ValueError(f"No chromosomes retained from {path}")
    return LensData(path=path, lengths=lengths, gene_counts=gene_counts, order=order, skipped=skipped)


def read_wgdi_gff(path: Path, allowed_chroms: Iterable[str]) -> Tuple[Dict[Tuple[str, int], Gene], Dict[str, Gene], Dict[str, List[Gene]], List[str]]:
    allowed = set(allowed_chroms)
    by_order: Dict[Tuple[str, int], Gene] = {}
    by_id: Dict[str, Gene] = {}
    by_chrom: Dict[str, List[Gene]] = {chrom: [] for chrom in allowed}
    warnings: List[str] = []

    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                warnings.append(f"{path}:{line_no}: expected WGDI gff with at least 6 columns")
                continue
            chrom, gene_id = parts[0], parts[1]
            if chrom not in allowed:
                continue
            try:
                start = int(float(parts[2]))
                end = int(float(parts[3]))
                order = int(float(parts[5]))
            except ValueError:
                warnings.append(f"{path}:{line_no}: invalid coordinate/order")
                continue
            if start > end:
                start, end = end, start
            gene = Gene(chrom=chrom, gene_id=gene_id, start=start, end=end, order=order)
            key = (chrom, order)
            if key in by_order:
                warnings.append(f"{path}:{line_no}: duplicate order {chrom}:{order}; keeping later record")
            by_order[key] = gene
            by_id[gene_id] = gene
            by_chrom.setdefault(chrom, []).append(gene)

    for genes in by_chrom.values():
        genes.sort(key=lambda g: (g.start, g.end, g.order))
    return by_order, by_id, by_chrom, warnings


def parse_order_list(value: object) -> List[int]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    orders: List[int] = []
    for token in re.split(r"[_;,\s]+", text):
        if not token:
            continue
        try:
            orders.append(int(float(token)))
        except ValueError:
            continue
    return orders


def optional_int(value: object) -> Optional[int]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def optional_float(value: object) -> Optional[float]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def clamp_interval(start: int, end: int, chrom_len: int) -> Tuple[int, int]:
    start = max(0, min(chrom_len, int(start)))
    end = max(0, min(chrom_len, int(end)))
    if start > end:
        start, end = end, start
    if start == end:
        end = min(chrom_len, start + 1)
        if start == end:
            start = max(0, end - 1)
    return start, end


def make_link(
    block_id: str,
    block_rank: int,
    source: str,
    gene1: Gene,
    gene2: Gene,
    lens1: LensData,
    lens2: LensData,
    name2: str,
    chrom_colors: Dict[str, str],
) -> LinkRecord:
    s1, e1 = clamp_interval(gene1.start, gene1.end, lens1.lengths[gene1.chrom])
    s2, e2 = clamp_interval(gene2.start, gene2.end, lens2.lengths[gene2.chrom])
    target_sid = sector_id(name2, gene2.chrom)
    color = chrom_colors.get(target_sid, "c01")
    return LinkRecord(
        block_id=str(block_id),
        block_rank=block_rank,
        source=source,
        chr1=gene1.chrom,
        start1=s1,
        end1=e1,
        gene_order1=gene1.order,
        gene_id1=gene1.gene_id,
        chr2=gene2.chrom,
        start2=s2,
        end2=e2,
        gene_order2=gene2.order,
        gene_id2=gene2.gene_id,
        color=color,
    )


def summarize_block(block_id: str, source: str, chr1: str, chr2: str, original_start1: Optional[int], original_end1: Optional[int], original_start2: Optional[int], original_end2: Optional[int], original_anchor_count: int, links: Sequence[LinkRecord], failed: int, pvalue: Optional[float] = None, ks_median: Optional[float] = None, note: str = "") -> BlockRecord:
    if links:
        mapped_start1 = min(link.start1 for link in links)
        mapped_end1 = max(link.end1 for link in links)
        mapped_start2 = min(link.start2 for link in links)
        mapped_end2 = max(link.end2 for link in links)
    else:
        mapped_start1 = mapped_end1 = mapped_start2 = mapped_end2 = None
    return BlockRecord(
        block_id=str(block_id),
        source=source,
        chr1=chr1,
        chr2=chr2,
        original_start1=original_start1,
        original_end1=original_end1,
        original_start2=original_start2,
        original_end2=original_end2,
        original_anchor_count=original_anchor_count,
        mapped_anchor_count=len(links),
        failed_anchor_count=failed,
        mapped_start1=mapped_start1,
        mapped_end1=mapped_end1,
        mapped_start2=mapped_start2,
        mapped_end2=mapped_end2,
        pvalue=pvalue,
        ks_median=ks_median,
        note=note,
    )


def parse_blockinfo_links(
    path: Path,
    by_order1: Dict[Tuple[str, int], Gene],
    by_order2: Dict[Tuple[str, int], Gene],
    lens1: LensData,
    lens2: LensData,
    name2: str,
    chrom_colors: Dict[str, str],
) -> Tuple[List[LinkRecord], List[BlockRecord], List[str]]:
    df = pd.read_csv(path)
    required = {"chr1", "chr2", "block1", "block2"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

    links: List[LinkRecord] = []
    blocks: List[BlockRecord] = []
    warnings: List[str] = []

    for rank, (_, row) in enumerate(df.iterrows(), 1):
        block_id = str(row.get("id", rank))
        chr1 = str(row["chr1"])
        chr2 = str(row["chr2"])
        orders1 = parse_order_list(row.get("block1"))
        orders2 = parse_order_list(row.get("block2"))
        original_anchor_count = optional_int(row.get("length")) or max(len(orders1), len(orders2))
        note = ""

        if chr1 not in lens1.lengths or chr2 not in lens2.lengths:
            note = "skipped_non_retained_chromosome"
            blocks.append(
                summarize_block(
                    block_id,
                    "block_information",
                    chr1,
                    chr2,
                    optional_int(row.get("start1")),
                    optional_int(row.get("end1")),
                    optional_int(row.get("start2")),
                    optional_int(row.get("end2")),
                    original_anchor_count,
                    [],
                    original_anchor_count,
                    optional_float(row.get("pvalue")),
                    optional_float(row.get("ks_median")),
                    note,
                )
            )
            continue

        if len(orders1) != len(orders2):
            warnings.append(f"{path.name} block {block_id}: block1/block2 anchor counts differ ({len(orders1)} vs {len(orders2)})")
        pair_count = min(len(orders1), len(orders2))
        block_links: List[LinkRecord] = []
        failed = abs(len(orders1) - len(orders2))

        for order1, order2 in zip(orders1[:pair_count], orders2[:pair_count]):
            gene1 = by_order1.get((chr1, order1))
            gene2 = by_order2.get((chr2, order2))
            if gene1 is None or gene2 is None:
                failed += 1
                continue
            block_links.append(make_link(block_id, rank, "block_information", gene1, gene2, lens1, lens2, name2, chrom_colors))

        links.extend(block_links)
        blocks.append(
            summarize_block(
                block_id,
                "block_information",
                chr1,
                chr2,
                optional_int(row.get("start1")),
                optional_int(row.get("end1")),
                optional_int(row.get("start2")),
                optional_int(row.get("end2")),
                original_anchor_count,
                block_links,
                failed,
                optional_float(row.get("pvalue")),
                optional_float(row.get("ks_median")),
                note,
            )
        )

    return links, blocks, warnings


def color_name(index: int, width: int = 2) -> str:
    return f"c{index:0{width}d}"


def assign_colors(order1: Sequence[str], order2: Sequence[str], name1: str, name2: str, palette_path: Optional[Path], repeat_per_genome: bool) -> Tuple[Dict[str, str], List[Tuple[str, Tuple[int, int, int]]]]:
    if repeat_per_genome:
        n_colors = max(len(order1), len(order2), 1)
    else:
        n_colors = max(len(order1) + len(order2), 1)
    palette = load_palette(palette_path, n_colors)
    width = max(2, len(str(n_colors)))
    color_defs = [(color_name(i + 1, width), rgb) for i, rgb in enumerate(palette)]
    chrom_colors: Dict[str, str] = {}

    for idx, chrom in enumerate(order1):
        cname = color_defs[idx % len(color_defs)][0]
        chrom_colors[sector_id(name1, chrom)] = cname
    for idx, chrom in enumerate(order2):
        color_idx = idx if repeat_per_genome else len(order1) + idx
        cname = color_defs[color_idx % len(color_defs)][0]
        chrom_colors[sector_id(name2, chrom)] = cname
    return chrom_colors, color_defs


def write_colors(path: Path, color_defs: Sequence[Tuple[str, Tuple[int, int, int]]], opacity: float) -> None:
    with path.open("w") as out:
        for cname, (r, g, b) in color_defs:
            out.write(f"{cname} = {r},{g},{b}\n")
            out.write(f"{cname}_t = {r},{g},{b},{opacity}\n")


def write_colors_fonts_patterns(path: Path) -> None:
    path.write_text(
        """<colors>
<<include etc/colors.conf>>
<<include colors.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<patterns>
<<include etc/patterns.conf>>
</patterns>
"""
    )


def write_karyotype(path: Path, lens1: LensData, name1: str, lens2: LensData, name2: str, chrom_colors: Dict[str, str]) -> List[str]:
    order_ids: List[str] = []
    with path.open("w") as out:
        for chrom in lens1.order:
            sid = sector_id(name1, chrom)
            order_ids.append(sid)
            out.write(f"chr - {sid} {chrom} 0 {lens1.lengths[chrom]} {chrom_colors[sid]}\n")
        for chrom in lens2.order:
            sid = sector_id(name2, chrom)
            order_ids.append(sid)
            out.write(f"chr - {sid} {chrom} 0 {lens2.lengths[chrom]} {chrom_colors[sid]}\n")
    return order_ids


def write_chrom_bands(path: Path, lens1: LensData, name1: str, lens2: LensData, name2: str, chrom_colors: Dict[str, str]) -> None:
    with path.open("w") as out:
        for lens, name in ((lens1, name1), (lens2, name2)):
            for chrom in lens.order:
                sid = sector_id(name, chrom)
                out.write(f"{sid} 0 {lens.lengths[chrom]} fill_color={chrom_colors[sid]}\n")


def write_chromosomes_tsv(path: Path, lens1: LensData, name1: str, lens2: LensData, name2: str, chrom_colors: Dict[str, str]) -> None:
    with path.open("w") as out:
        out.write("genome\tsector_id\tchromosome\tlength_bp\tgene_count\tcolor\n")
        for lens, name in ((lens1, name1), (lens2, name2)):
            for chrom in lens.order:
                sid = sector_id(name, chrom)
                count = lens.gene_counts.get(chrom)
                out.write(f"{name}\t{sid}\t{chrom}\t{lens.lengths[chrom]}\t{'' if count is None else count}\t{chrom_colors[sid]}\n")


def density_rows(lens: LensData, name: str, by_chrom: Dict[str, List[Gene]], chrom_colors: Dict[str, str], window: int, step: int) -> List[Tuple[str, int, int, int, str]]:
    rows: List[Tuple[str, int, int, int, str]] = []
    for chrom in lens.order:
        genes = by_chrom.get(chrom, [])
        starts = sorted(g.start for g in genes)
        ends = sorted(g.end for g in genes)
        sid = sector_id(name, chrom)
        color = chrom_colors[sid]
        pos = 0
        chrom_len = lens.lengths[chrom]
        while pos < chrom_len:
            end = min(chrom_len, pos + window)
            n_start_before_end = bisect.bisect_right(starts, end)
            n_end_before_start = bisect.bisect_left(ends, pos)
            count = max(0, n_start_before_end - n_end_before_start)
            rows.append((sid, pos, max(pos + 1, end), count, color))
            pos += step
    return rows


def write_gene_density(path: Path, rows: Sequence[Tuple[str, int, int, int, str]]) -> int:
    values = [row[3] for row in rows]
    with path.open("w") as out:
        for sid, start, end, count, color in rows:
            out.write(f"{sid} {start} {end} {count} fill_color={color}\n")
    if not values:
        return 1
    values_sorted = sorted(values)
    idx = min(len(values_sorted) - 1, max(0, math.ceil(len(values_sorted) * 0.99) - 1))
    return max(1, int(math.ceil(values_sorted[idx] * 1.15)))


def block_links_from_blocks(
    blocks: Sequence[BlockRecord],
    name2: str,
    chrom_colors: Dict[str, str],
    min_block_anchors: int,
) -> List[BlockLink]:
    links: List[BlockLink] = []
    for rank, block in enumerate(blocks, 1):
        if block.mapped_anchor_count < min_block_anchors:
            continue
        if (
            block.mapped_start1 is None
            or block.mapped_end1 is None
            or block.mapped_start2 is None
            or block.mapped_end2 is None
        ):
            continue
        target_sid = sector_id(name2, block.chr2)
        links.append(
            BlockLink(
                block_id=block.block_id,
                block_rank=rank,
                source=block.source,
                chr1=block.chr1,
                start1=block.mapped_start1,
                end1=block.mapped_end1,
                chr2=block.chr2,
                start2=block.mapped_start2,
                end2=block.mapped_end2,
                original_anchor_count=block.original_anchor_count,
                mapped_anchor_count=block.mapped_anchor_count,
                color=chrom_colors.get(target_sid, "c01"),
            )
        )
    return links


def block_link_span(link: BlockLink) -> int:
    return max(link.end1 - link.start1, link.end2 - link.start2)


def limit_block_links(links: Sequence[BlockLink], max_links: int) -> Tuple[List[BlockLink], int]:
    if max_links <= 0 or len(links) <= max_links:
        return list(links), 0
    ranked = sorted(
        enumerate(links),
        key=lambda item: (item[1].mapped_anchor_count, block_link_span(item[1]), -item[0]),
        reverse=True,
    )
    retained_indexes = {idx for idx, _ in ranked[:max_links]}
    retained = [link for idx, link in enumerate(links) if idx in retained_indexes]
    return retained, len(links) - len(retained)


def block_link_keys(links: Sequence[BlockLink]) -> set:
    return {(link.block_id, link.source, link.chr1, link.chr2) for link in links}


def write_links(path: Path, links: Sequence[BlockLink], name1: str, name2: str) -> None:
    with path.open("w") as out:
        for link in links:
            out.write(
                f"{sector_id(name1, link.chr1)} {link.start1} {link.end1} "
                f"{sector_id(name2, link.chr2)} {link.start2} {link.end2} "
                f"color={link.color}_t\n"
            )


def write_block_links_tsv(path: Path, links: Sequence[BlockLink]) -> None:
    with path.open("w") as out:
        out.write(
            "block_id\tblock_rank\tsource\tchr1\tstart1\tend1\tlength1_bp\t"
            "chr2\tstart2\tend2\tlength2_bp\toriginal_anchor_count\tmapped_anchor_count\tcolor\n"
        )
        for link in links:
            out.write(
                f"{link.block_id}\t{link.block_rank}\t{link.source}\t{link.chr1}\t{link.start1}\t{link.end1}\t"
                f"{link.end1 - link.start1}\t{link.chr2}\t{link.start2}\t{link.end2}\t{link.end2 - link.start2}\t"
                f"{link.original_anchor_count}\t{link.mapped_anchor_count}\t{link.color}\n"
            )


def write_anchor_links_tsv(path: Path, links: Sequence[LinkRecord]) -> None:
    with path.open("w") as out:
        out.write(
            "block_id\tblock_rank\tsource\tchr1\tstart1\tend1\tgene_order1\tgene_id1\t"
            "chr2\tstart2\tend2\tgene_order2\tgene_id2\tcolor\n"
        )
        for link in links:
            out.write(
                f"{link.block_id}\t{link.block_rank}\t{link.source}\t{link.chr1}\t{link.start1}\t{link.end1}\t"
                f"{'' if link.gene_order1 is None else link.gene_order1}\t{link.gene_id1}\t"
                f"{link.chr2}\t{link.start2}\t{link.end2}\t{'' if link.gene_order2 is None else link.gene_order2}\t"
                f"{link.gene_id2}\t{link.color}\n"
            )


def block_orientation(block: BlockRecord) -> str:
    if (
        block.original_start1 is None
        or block.original_end1 is None
        or block.original_start2 is None
        or block.original_end2 is None
    ):
        return "unknown"
    dir1 = 1 if block.original_end1 >= block.original_start1 else -1
    dir2 = 1 if block.original_end2 >= block.original_start2 else -1
    return "same" if dir1 == dir2 else "inverted"


def write_blocks_tsv(path: Path, blocks: Sequence[BlockRecord]) -> None:
    with path.open("w") as out:
        out.write(
            "block_id\tsource\tchr1\tchr2\toriginal_start1\toriginal_end1\toriginal_start2\toriginal_end2\t"
            "orientation\toriginal_anchor_count\tmapped_anchor_count\tfailed_anchor_count\tmapped_start1\tmapped_end1\t"
            "mapped_length1_bp\tmapped_start2\tmapped_end2\tmapped_length2_bp\tpvalue\tks_median\tnote\n"
        )
        for block in blocks:
            mapped_length1 = None
            mapped_length2 = None
            if block.mapped_start1 is not None and block.mapped_end1 is not None:
                mapped_length1 = block.mapped_end1 - block.mapped_start1
            if block.mapped_start2 is not None and block.mapped_end2 is not None:
                mapped_length2 = block.mapped_end2 - block.mapped_start2
            vals = [
                block.block_id,
                block.source,
                block.chr1,
                block.chr2,
                block.original_start1,
                block.original_end1,
                block.original_start2,
                block.original_end2,
                block_orientation(block),
                block.original_anchor_count,
                block.mapped_anchor_count,
                block.failed_anchor_count,
                block.mapped_start1,
                block.mapped_end1,
                mapped_length1,
                block.mapped_start2,
                block.mapped_end2,
                mapped_length2,
                block.pvalue,
                block.ks_median,
                block.note,
            ]
            out.write("\t".join("" if v is None else str(v) for v in vals) + "\n")


def write_pair_summary(path: Path, links: Sequence[BlockLink]) -> None:
    summary: Dict[Tuple[str, str], int] = {}
    for link in links:
        summary[(link.chr1, link.chr2)] = summary.get((link.chr1, link.chr2), 0) + 1
    with path.open("w") as out:
        out.write("chr1\tchr2\tlinks\n")
        for (chr1, chr2), count in sorted(summary.items(), key=lambda x: (-x[1], x[0][0], x[0][1])):
            out.write(f"{chr1}\t{chr2}\t{count}\n")


def find_housekeeping() -> Optional[Path]:
    candidates: List[Path] = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "etc" / "housekeeping.conf")
    circos_path = shutil.which("circos")
    if circos_path:
        candidates.append(Path(circos_path).resolve().parent.parent / "etc" / "housekeeping.conf")
    candidates.extend(
        [
            Path("/media/desk16/tl5024/miniconda3/envs/circos/etc/housekeeping.conf"),
            Path("/media/desk16/tl5024/miniconda3/envs/circos/bin/../etc/housekeeping.conf"),
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def write_housekeeping(path: Path, max_links: int) -> None:
    source = find_housekeeping()
    if source:
        text = source.read_text()
        if re.search(r"^\s*max_links\s*=", text, flags=re.M):
            text = re.sub(r"^\s*max_links\s*=.*", f"max_links = {max_links}", text, flags=re.M)
        else:
            text += f"\nmax_links = {max_links}\n"
        path.write_text(text)
    else:
        path.write_text(f"max_links = {max_links}\n")


def write_circos_conf(path: Path, order_ids: Sequence[str], output_prefix: str, args: argparse.Namespace, hist_max: int) -> None:
    conf = f"""
<<include colors_fonts_patterns.custom.conf>>

karyotype = karyotype.txt

chromosomes_units = 1000000
chromosomes_display_default = yes
chromosomes_order = {','.join(order_ids)}

angle_offset = {args.angle_offset}

<image>
dir = .
file = {output_prefix}.png
png = yes
svg = yes
radius = {args.image_radius}
background = white
</image>

<ideogram>
<spacing>
default = {args.spacing}
</spacing>

radius           = {args.ideogram_radius}
thickness        = 1p
fill             = no
stroke_thickness = 0p

show_label       = yes
label_font       = default
label_radius     = {args.label_radius}
label_size       = {args.label_size}
label_parallel   = yes
</ideogram>

<highlights>
<highlight>
file             = chrom_bands.txt
r0               = {args.band_r0}
r1               = {args.band_r1}
stroke_color     = white
stroke_thickness = {args.band_stroke}
</highlight>
</highlights>

<plots>
<plot>
type             = histogram
file             = gene_density.txt
r0               = {args.hist_r0}
r1               = {args.hist_r1}
min              = 0
max              = {hist_max}
fill_under       = yes
stroke_thickness = 0p
orientation      = out
</plot>
</plots>

<links>
<link>
file             = links.txt
radius           = {args.link_radius}
bezier_radius    = {args.bezier_radius}
thickness        = {args.link_thickness}
ribbon           = yes
flat             = yes
</link>
</links>

<ticks>
show_ticks = no
show_tick_labels = no
</ticks>

<<include housekeeping.conf>>
"""
    path.write_text(conf.lstrip())


def validate_coordinates(links: Sequence[BlockLink], lens1: LensData, lens2: LensData) -> List[str]:
    errors: List[str] = []
    for idx, link in enumerate(links, 1):
        len1 = lens1.lengths.get(link.chr1)
        len2 = lens2.lengths.get(link.chr2)
        if len1 is None or len2 is None:
            errors.append(f"link {idx}: chromosome missing from retained lens")
            continue
        if not (0 <= link.start1 <= link.end1 <= len1):
            errors.append(f"link {idx}: chr1 coordinate out of range")
        if not (0 <= link.start2 <= link.end2 <= len2):
            errors.append(f"link {idx}: chr2 coordinate out of range")
        if len(errors) >= 20:
            errors.append("additional coordinate errors suppressed")
            break
    return errors


def validate_block_link_coverage(links: Sequence[BlockLink], blocks: Sequence[BlockRecord], exact: bool) -> List[str]:
    if not exact:
        return []
    mapped_blocks = [
        block
        for block in blocks
        if block.mapped_anchor_count > 0
        and block.mapped_start1 is not None
        and block.mapped_end1 is not None
        and block.mapped_start2 is not None
        and block.mapped_end2 is not None
    ]
    errors: List[str] = []
    if len(links) != len(mapped_blocks):
        errors.append(f"block ribbon link count {len(links)} != mapped block count {len(mapped_blocks)}")
    seen = {
        (
            link.block_id,
            link.source,
            link.chr1,
            link.start1,
            link.end1,
            link.chr2,
            link.start2,
            link.end2,
        )
        for link in links
    }
    for block in mapped_blocks[:]:
        key = (
            block.block_id,
            block.source,
            block.chr1,
            block.mapped_start1,
            block.mapped_end1,
            block.chr2,
            block.mapped_start2,
            block.mapped_end2,
        )
        if key not in seen:
            errors.append(f"mapped block {block.block_id} missing exact ribbon link")
            if len(errors) >= 20:
                errors.append("additional block-link coverage errors suppressed")
                break
    return errors


def main_block_warnings(blocks: Sequence[BlockRecord], top_n: int = 20) -> List[str]:
    ranked = sorted(blocks, key=lambda b: b.original_anchor_count, reverse=True)[:top_n]
    return [
        f"top block {block.block_id} ({block.chr1}-{block.chr2}, anchors={block.original_anchor_count}) has no mapped anchors"
        for block in ranked
        if block.mapped_anchor_count == 0
    ]


def main_block_retention_warnings(blocks: Sequence[BlockRecord], links: Sequence[BlockLink], top_n: int = 20) -> List[str]:
    retained = block_link_keys(links)
    ranked = sorted(blocks, key=lambda b: b.original_anchor_count, reverse=True)[:top_n]
    return [
        f"top block {block.block_id} ({block.chr1}-{block.chr2}, anchors={block.original_anchor_count}) is not retained in drawn ribbons"
        for block in ranked
        if (block.block_id, block.source, block.chr1, block.chr2) not in retained
    ]


def write_validation_report(
    path: Path,
    spec: DatasetSpec,
    lens1: LensData,
    lens2: LensData,
    block_links: Sequence[BlockLink],
    all_block_links: Sequence[BlockLink],
    anchor_links: Sequence[LinkRecord],
    blocks: Sequence[BlockRecord],
    warnings: Sequence[str],
    coordinate_errors: Sequence[str],
    block_link_errors: Sequence[str],
    min_block_anchors: int,
    threshold_block_link_count: int,
    max_links: int,
    max_links_removed: int,
    render_returncode: Optional[int],
) -> None:
    mapped_blocks = sum(1 for block in blocks if block.mapped_anchor_count > 0)
    failed_anchors = sum(block.failed_anchor_count for block in blocks)
    text = [
        f"# Validation Report: {spec.label}",
        "",
        "## Inputs",
        "- source: block_information",
        f"- gff1: {spec.gff1}",
        f"- len1: {spec.len1}",
        f"- gff2: {spec.gff2}",
        f"- len2: {spec.len2}",
        f"- block_info: {spec.block_info}",
        "",
        "## Counts",
        f"- genome1 chromosomes: {len(lens1.order)}",
        f"- genome2 chromosomes: {len(lens2.order)}",
        f"- total Circos sectors: {len(lens1.order) + len(lens2.order)}",
        f"- blocks total: {len(blocks)}",
        f"- blocks mapped: {mapped_blocks}",
        f"- min_block_anchors: {min_block_anchors}",
        f"- block ribbons after min_block_anchors: {threshold_block_link_count}",
        f"- max_links: {max_links}",
        f"- block ribbons removed by max_links: {max_links_removed}",
        f"- block ribbon links written: {len(block_links)}",
        f"- all mapped block links audited: {len(all_block_links)}",
        f"- anchor pairs audited: {len(anchor_links)}",
        f"- failed anchors: {failed_anchors}",
        f"- drawn model: one interval ribbon per mapped synteny block",
        "",
        "## Coordinate Validation",
        f"- coordinate errors: {len(coordinate_errors)}",
    ]
    if coordinate_errors:
        text.extend(f"  - {msg}" for msg in coordinate_errors[:20])
    text.extend(["", "## Block Ribbon Validation", f"- block-link coverage errors: {len(block_link_errors)}"])
    text.extend(f"  - {msg}" for msg in block_link_errors[:20])
    text.extend(["", "## Warnings", f"- warnings: {len(warnings)}"])
    text.extend(f"  - {msg}" for msg in warnings[:50])
    main_warnings = main_block_warnings(blocks)
    text.extend(["", "## Main Block Check", f"- top block mapping warnings: {len(main_warnings)}"])
    text.extend(f"  - {msg}" for msg in main_warnings)
    if render_returncode is not None:
        text.extend(["", "## Render", f"- circos return code: {render_returncode}"])
    path.write_text("\n".join(text) + "\n")


def validate_spec_paths(spec: DatasetSpec) -> None:
    for path in [spec.gff1, spec.len1, spec.gff2, spec.len2]:
        if not path.exists():
            raise FileNotFoundError(path)
    if spec.block_info is None:
        raise FileNotFoundError(f"No block_information/blockinfo CSV found for {spec.label}")
    if not spec.block_info.exists():
        raise FileNotFoundError(spec.block_info)


def render_circos(outdir: Path, args: argparse.Namespace) -> int:
    cmd = [args.conda_exe, "run", "-n", args.conda_env, "circos", "-conf", "circos.conf"]
    if args.no_conda_run_render:
        cmd = ["circos", "-conf", "circos.conf"]
    proc = subprocess.run(cmd, cwd=outdir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (outdir / "circos.render.log").write_text(proc.stdout)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
    return proc.returncode


def check_dependencies(args: argparse.Namespace) -> None:
    """Fail early if the expected circos environment is not usable."""
    messages = [f"pandas {pd.__version__}"]

    if args.no_conda_run_render:
        circos_exe = shutil.which("circos")
        if circos_exe is None:
            raise SystemExit("Missing executable: circos. Install Perl Circos in the active environment.")
        cmd = [circos_exe, "-v"]
    else:
        conda_exe = Path(args.conda_exe)
        if not conda_exe.exists():
            raise SystemExit(f"Missing conda executable: {conda_exe}")
        cmd = [str(conda_exe), "run", "-n", args.conda_env, "circos", "-v"]

    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise SystemExit(
            "Circos dependency check failed. Use the circos conda environment and install missing "
            "packages there. Command output:\n" + proc.stdout
        )
    version_line = proc.stdout.strip().splitlines()[0] if proc.stdout.strip() else "circos ok"
    messages.append(version_line)
    if not args.quiet_deps:
        print("Dependency check OK: " + "; ".join(messages), file=sys.stderr)


def generate_one(spec: DatasetSpec, outdir: Path, args: argparse.Namespace) -> Dict[str, object]:
    validate_spec_paths(spec)
    outdir.mkdir(parents=True, exist_ok=True)

    exclude_re = re.compile(args.exclude_regex, re.I) if args.exclude_regex else None
    keep_re = re.compile(args.keep_regex) if args.keep_regex else None

    lens1 = read_lens(spec.len1, exclude_re, keep_re)
    lens2 = read_lens(spec.len2, exclude_re, keep_re)
    chrom_colors, color_defs = assign_colors(lens1.order, lens2.order, spec.name1, spec.name2, args.palette, not args.sequential_colors)

    by_order1, _by_id1, by_chrom1, gff_warnings1 = read_wgdi_gff(spec.gff1, lens1.order)
    by_order2, _by_id2, by_chrom2, gff_warnings2 = read_wgdi_gff(spec.gff2, lens2.order)

    warnings: List[str] = []
    warnings.extend(gff_warnings1)
    warnings.extend(gff_warnings2)
    if lens1.skipped:
        warnings.append(f"{spec.len1}: skipped {len(lens1.skipped)} non-retained entries")
    if lens2.skipped:
        warnings.append(f"{spec.len2}: skipped {len(lens2.skipped)} non-retained entries")

    anchor_links, blocks, parse_warnings = parse_blockinfo_links(
        spec.block_info, by_order1, by_order2, lens1, lens2, spec.name2, chrom_colors
    )
    warnings.extend(parse_warnings)

    all_block_links = block_links_from_blocks(blocks, spec.name2, chrom_colors, 1)
    threshold_block_links = block_links_from_blocks(blocks, spec.name2, chrom_colors, args.min_block_anchors)
    block_links, max_links_removed = limit_block_links(threshold_block_links, args.max_links)
    exact_block_coverage = args.min_block_anchors <= 1 and not args.max_links
    if args.min_block_anchors > 1:
        warnings.append(f"summary mode: retained blocks with at least {args.min_block_anchors} mapped anchors")
    if max_links_removed:
        warnings.append(
            f"summary mode: kept the {len(block_links)} largest block ribbons after --max-links {args.max_links}"
        )

    order_ids = write_karyotype(outdir / "karyotype.txt", lens1, spec.name1, lens2, spec.name2, chrom_colors)
    write_chrom_bands(outdir / "chrom_bands.txt", lens1, spec.name1, lens2, spec.name2, chrom_colors)
    write_chromosomes_tsv(outdir / "chromosomes.tsv", lens1, spec.name1, lens2, spec.name2, chrom_colors)
    write_colors(outdir / "colors.conf", color_defs, args.link_opacity)
    write_colors_fonts_patterns(outdir / "colors_fonts_patterns.custom.conf")
    density = density_rows(lens1, spec.name1, by_chrom1, chrom_colors, args.gene_density_window, args.gene_density_step)
    density.extend(density_rows(lens2, spec.name2, by_chrom2, chrom_colors, args.gene_density_window, args.gene_density_step))
    hist_max = write_gene_density(outdir / "gene_density.txt", density)
    write_links(outdir / "links.full.txt", block_links, spec.name1, spec.name2)
    write_links(outdir / "links.txt", block_links, spec.name1, spec.name2)
    write_links(outdir / "links.all_blocks.txt", all_block_links, spec.name1, spec.name2)
    write_block_links_tsv(outdir / "block_links.tsv", block_links)
    write_block_links_tsv(outdir / "all_block_links.tsv", all_block_links)
    write_anchor_links_tsv(outdir / "anchor_links.tsv", anchor_links)
    write_blocks_tsv(outdir / "mapped_blocks.tsv", blocks)
    write_pair_summary(outdir / "links.pair_summary.tsv", block_links)
    write_housekeeping(outdir / "housekeeping.conf", max(len(block_links) * 2, 100000))
    write_circos_conf(outdir / "circos.conf", order_ids, spec.output_prefix, args, hist_max)

    coordinate_errors = validate_coordinates(block_links, lens1, lens2)
    block_link_errors = validate_block_link_coverage(block_links, blocks, exact_block_coverage)
    warnings.extend(main_block_warnings(blocks))
    warnings.extend(main_block_retention_warnings(blocks, block_links))
    render_returncode: Optional[int] = None
    if not args.no_render:
        render_returncode = render_circos(outdir, args)

    write_validation_report(
        outdir / "validation_report.md",
        spec,
        lens1,
        lens2,
        block_links,
        all_block_links,
        anchor_links,
        blocks,
        warnings,
        coordinate_errors,
        block_link_errors,
        args.min_block_anchors,
        len(threshold_block_links),
        args.max_links,
        max_links_removed,
        render_returncode,
    )

    return {
        "label": spec.label,
        "outdir": outdir,
        "sectors": len(order_ids),
        "block_links": len(block_links),
        "all_block_links": len(all_block_links),
        "anchor_links": len(anchor_links),
        "blocks": len(blocks),
        "mapped_blocks": sum(1 for block in blocks if block.mapped_anchor_count > 0),
        "coordinate_errors": len(coordinate_errors),
        "block_link_errors": len(block_link_errors),
        "warnings": len(warnings),
        "render_returncode": render_returncode,
    }


def single_spec_from_args(args: argparse.Namespace) -> DatasetSpec:
    required = ["gff1", "len1", "name1", "gff2", "len2", "name2", "block_info"]
    missing = [arg for arg in required if getattr(args, arg) is None]
    if missing:
        raise SystemExit(f"Missing required single-plot arguments: {', '.join('--' + x.replace('_', '-') for x in missing)}")
    name = args.label or f"{args.name1}_vs_{args.name2}"
    return DatasetSpec(
        label=name,
        gff1=args.gff1,
        len1=args.len1,
        name1=args.name1,
        gff2=args.gff2,
        len2=args.len2,
        name2=args.name2,
        block_info=args.block_info,
        output_prefix=args.output_prefix or f"{safe_id(args.name1)}_vs_{safe_id(args.name2)}_circos",
    )


def write_batch_summary(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w") as out:
        out.write(
            "label\toutdir\tsectors\tblock_links\tall_block_links\tanchor_links\tblocks\tmapped_blocks\t"
            "coordinate_errors\tblock_link_errors\twarnings\trender_returncode\n"
        )
        for row in rows:
            out.write(
                f"{row['label']}\t{row['outdir']}\t{row['sectors']}\t{row['block_links']}\t"
                f"{row['all_block_links']}\t{row['anchor_links']}\t{row['blocks']}\t{row['mapped_blocks']}\t"
                f"{row['coordinate_errors']}\t{row['block_link_errors']}\t{row['warnings']}\t"
                f"{'' if row['render_returncode'] is None else row['render_returncode']}\n"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Map WGDI-format gff/len plus blockinfo CSV gene-order coordinates "
            "to real bp coordinates and render a block-level Circos plot."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--outdir", type=Path, default=Path("circos_plot/results"), help="Output directory.")
    parser.add_argument("--no-render", action="store_true", help="Only write Circos input/config and audit files; do not run Perl Circos.")

    parser.add_argument("--gff1", type=Path, help="Genome A WGDI-format gff file.")
    parser.add_argument("--len1", type=Path, help="Genome A WGDI len file.")
    parser.add_argument("--name1", help="Display/prefix name for genome A sectors.")
    parser.add_argument("--gff2", type=Path, help="Genome B WGDI-format gff file.")
    parser.add_argument("--len2", type=Path, help="Genome B WGDI len file.")
    parser.add_argument("--name2", help="Display/prefix name for genome B sectors.")
    parser.add_argument("--block-info", type=Path, help="WGDI blockinfo/block_information CSV file.")
    parser.add_argument("--label", help="Dataset label used in reports.")
    parser.add_argument("--output-prefix", help="Output image prefix inside --outdir.")

    parser.add_argument("--palette", type=Path, help="Optional RGB/hex palette file.")
    parser.add_argument("--sequential-colors", action="store_true", help="Use unique sequential colors across both genome halves instead of repeating per genome.")
    parser.add_argument("--exclude-regex", default=DEFAULT_EXCLUDE_RE, help="Regex for non-chromosomal sequences to remove.")
    parser.add_argument("--keep-regex", help="Optional regex; if set, retain only matching sequence names before applying exclude regex.")
    parser.add_argument("--min-block-anchors", type=int, default=50, help="Draw only blocks with at least this many mapped anchors; all mapped blocks are still written to audit files.")
    parser.add_argument("--max-links", type=int, default=0, help="After --min-block-anchors, keep only the largest N block ribbons. 0 keeps all threshold-passing blocks.")

    parser.add_argument("--gene-density-window", type=int, default=1_000_000)
    parser.add_argument("--gene-density-step", type=int, default=1_000_000)
    parser.add_argument("--link-opacity", type=float, default=0.22)
    parser.add_argument("--link-radius", default="0.850r")
    parser.add_argument("--bezier-radius", default="0.08r", help="Bezier control radius for block ribbons. Smaller values keep ribbons from forming broad central bands.")
    parser.add_argument("--link-thickness", default="1p")
    parser.add_argument("--image-radius", default="1600p")
    parser.add_argument("--ideogram-radius", default="0.885r")
    parser.add_argument("--band-r0", default="0.865r")
    parser.add_argument("--band-r1", default="0.918r")
    parser.add_argument("--band-stroke", default="3p")
    parser.add_argument("--hist-r0", default="0.928r")
    parser.add_argument("--hist-r1", default="0.995r")
    parser.add_argument("--label-radius", default="1.068r")
    parser.add_argument("--label-size", default="22p")
    parser.add_argument("--chrom-spacing", "--spacing", dest="spacing", default="0.006r", help="Circos spacing between adjacent chromosome sectors, e.g. 0.004r for tighter or 0.012r for wider gaps.")
    parser.add_argument("--angle-offset", default="-90")

    parser.add_argument("--conda-exe", default="/media/desk16/tl5024/miniconda3/bin/conda")
    parser.add_argument("--conda-env", default="circos")
    parser.add_argument("--no-conda-run-render", action="store_true", help="Render by calling `circos` directly instead of `conda run -n circos circos`.")
    parser.add_argument("--quiet-deps", action="store_true", help="Do not print successful dependency check details.")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.min_block_anchors < 1:
        raise SystemExit("--min-block-anchors must be >= 1")
    if args.max_links < 0:
        raise SystemExit("--max-links must be >= 0")
    if not (0 < args.link_opacity <= 1):
        raise SystemExit("--link-opacity must be > 0 and <= 1")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)
    check_dependencies(args)
    spec = single_spec_from_args(args)
    outdir = args.outdir if args.outdir.is_absolute() else Path.cwd() / args.outdir
    rows = [generate_one(spec, outdir, args)]
    write_batch_summary(outdir / "VALIDATION_SUMMARY.tsv", rows)

    for row in rows:
        print(
            f"{row['label']}: sectors={row['sectors']} block_links={row['block_links']} "
            f"all_block_links={row['all_block_links']} anchor_pairs={row['anchor_links']} "
            f"blocks={row['mapped_blocks']}/{row['blocks']} "
            f"coord_errors={row['coordinate_errors']} block_link_errors={row['block_link_errors']} "
            f"render={row['render_returncode']}"
        )
    if any(row["coordinate_errors"] for row in rows):
        return 2
    if any(row["block_link_errors"] for row in rows):
        return 2
    if any(row["render_returncode"] not in (None, 0) for row in rows):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
