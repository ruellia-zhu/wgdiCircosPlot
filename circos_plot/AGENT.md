# AGENT.md - WGDI blockinfo Circos 迭代规范

本目录以 `wgdi_blockinfo_to_circos.py` 为入口，使用 Perl Circos 渲染。Python 只做 WGDI 输入解析、基因序号到真实 bp 坐标映射、Circos 文件生成和验收报告；不要在 Python 中手写绘图引擎。

主图必须是共线性块区间到区间的 ribbon：一个 WGDI block 对应 `links.txt` 中一条 interval ribbon。不要把所有 anchor gene pair 画成点对点细线。

## 固定事实约束

- 用户必须显式提供 `--gff1 --len1 --name1 --gff2 --len2 --name2 --block-info --outdir`。
- 不提供任何内置批量入口。批量绘图用 shell 循环逐个调用脚本。
- 不解析 `.collinearity.txt` fallback；缺少 blockinfo 直接失败。
- WGDI `.len` 第 2 列是真实染色体 bp 长度；染色体顺序必须严格等于 `.len` 文件顺序。
- WGDI `.gff` 第 6 列是基因序号；`blockinfo.csv` 的 `block1/block2` 必须通过该序号映射回真实 bp 坐标。
- 默认主图绘制 `mapped_anchor_count >= 50` 的主要共线性块。全量可映射 block 必须保存在 `all_block_links.tsv` 和 `links.all_blocks.txt`。
- 默认染色体间距为 `--chrom-spacing 0.006r`；`--spacing` 是兼容别名。
- 默认 ribbon 曲率为 `--bezier-radius 0.08r`，用于减轻区间 ribbon 在中心变粗的视觉效果。
- `anchor_links.tsv` 只作审计，不参与绘图。
- 不要删除主染色体，不要只保留有 link 的染色体，不要用视觉效果覆盖事实错误。

## block 筛选规则

主图筛选顺序固定：

1. 解析和映射所有可映射 block，写入全量审计文件。
2. 应用 `--min-block-anchors`，默认 `50`。
3. 如果设置 `--max-links N` 且阈值后仍超过 N，则保留最大的 N 条 block ribbon。

最大 block 的排序规则：`mapped_anchor_count` 多者优先；相同则映射 bp span 更长者优先；再相同则保留 blockinfo 中更靠前者。

不要为了让中心更干净而把主图改成 anchor pair 细线或伪坐标。link 中心看似变粗通常是 interval ribbon 的 Bezier 几何和透明叠加造成的；优先调整 `--bezier-radius`、`--min-block-anchors` 或透明度。

## 当前 6 组输入

- `01_guohuai_vs_guohuai`: `Sjap.gff` + `Sjap.len` + `Sjap-Sjap_block_information.csv` -> `SjapL/SjapR`
- `02_jiyehuai_vs_jiyehuai`: `S.jap.o.gff` + `S.jap.o.len` + `Sjapo-Sjapo_blcok_information.csv` -> `SjapoL/SjapoR`
- `03_xianghuai_vs_xianghuai`: `Cwilsonii.gff` + `Cwilsonii.len` + `Cwilsonii-Cwilsonii_block_information.csv` -> `CwilsoniiL/CwilsoniiR`
- `04_cjxh_vs_cjxh`: `Pplatycarpum.gff` + `Pplatycarpum.len` + `Pplatycarpum-Pplatycarpum_block_information.csv` -> `PplatycarpumL/PplatycarpumR`
- `05_guohuai_vs_jiyehuai`: `jiyehuai.gff` + `jiyehuai.len` vs `Sjap.gff` + `Sjap.len` + `Sjapo-Sjap.blockinfo.csv` -> `Sjapo/Sjap`
- `06_guohuai_vs_Medicago`: `Sjap.gff` + `Sjap.len` vs `Mtrun.gff` + `Mtrun.len` + `Sjap-Mtrun_block_information.csv` -> `Sjap/Mtrun`

## 环境和渲染

脚本默认渲染。只生成配置时使用 `--no-render`。

渲染通过：

```bash
conda run -n circos circos -conf circos.conf
```

不要直接调用 `/media/desk16/tl5024/miniconda3/envs/circos/bin/circos`，直接调用可能找不到 Perl 模块。

## 验收标准

每次改动后至少检查：

```bash
OUT=circos_plot/results_v3_spacing006_bezier008
for d in "$OUT"/*
do
  printf '%s ' "$(basename "$d")"
  awk '$1=="chr"{n++} END{printf "sectors=%d ", n+0}' "$d/karyotype.txt"
  printf 'drawn_block_links='
  wc -l < "$d/links.txt"
  printf 'all_block_links='
  tail -n +2 "$d/all_block_links.tsv" | wc -l
done
```

预期 sector 数：

- `01-05`: 28
- `06`: 22

检查坐标错误：

```bash
cut -f1-12 "$OUT"/*/VALIDATION_SUMMARY.tsv
```

`coordinate_errors` 必须为 0。`validation_report.md` 中排名前 20 的主要 block 不能出现 unmapped 警告，除非输入文件本身缺失对应 gene order，并且报告中已说明。

## 禁止复用的旧流程

正式结果不要使用旧的筛选脚本：

```text
code/filter_circos_links_by_pair.py
code/prune_circos_to_linked_chromosomes.py
```

这些脚本会过滤染色体对或删除扇区，容易造成染色体数量和物种结构错误。
