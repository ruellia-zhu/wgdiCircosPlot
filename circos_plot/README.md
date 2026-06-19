# WGDI blockinfo Circos 绘图脚本

`wgdi_blockinfo_to_circos.py` 用 WGDI 格式的 `gff`、`len` 和 `blockinfo/block_information.csv` 生成基因组共线性 Circos 图。绘图后端是成熟的 Perl Circos；Python 脚本只负责解析 WGDI 文件、把 block 中的基因序号映射成真实染色体 bp 坐标、写出 Circos 配置并默认渲染。

主图绘制的是 **block-to-block 区间 ribbon**：一个 WGDI 共线性 block 在 `links.txt` 中对应一条染色体区间到染色体区间的 ribbon。脚本不会把所有 anchor gene pair 画成点对点细线；anchor pair 只写入 `anchor_links.tsv` 作为审计。

## 环境

建议使用已有的 `circos` conda 环境：

```bash
cd /media/desk16/tl5024/work/wgdi
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --help
conda run -n circos circos -v
```

需要的依赖：

- Python: `pandas`
- Perl Circos: `circos`
- 常见 Circos Perl 依赖：`perl-gd`、`perl-font-ttf`、`perl-math-round`、`perl-config-general`、`perl-clone`、`perl-math-bezier`

如缺依赖，只安装到 `circos` 环境：

```bash
conda install -n circos -c bioconda -c conda-forge \
  circos pandas perl-gd perl-font-ttf perl-math-round \
  perl-config-general perl-clone perl-math-bezier
```

默认渲染时脚本调用 `conda run -n circos circos -conf circos.conf`。不要直接调用环境里的 `bin/circos`，直接调用有时找不到 Perl 模块。

## 单图用法

用户只需要提供基因组 A、基因组 B 的 WGDI 文件、名称、blockinfo 和输出路径：

```bash
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py \
  --gff1 06_guohuai_vs_Medicago/Sjap.gff \
  --len1 06_guohuai_vs_Medicago/Sjap.len \
  --name1 Sjap \
  --gff2 06_guohuai_vs_Medicago/Mtrun.gff \
  --len2 06_guohuai_vs_Medicago/Mtrun.len \
  --name2 Mtrun \
  --block-info 06_guohuai_vs_Medicago/Sjap-Mtrun_block_information.csv \
  --outdir circos_plot/results_v3_spacing006_bezier008/06_Sjap_vs_Mtrun \
  --label 06_Sjap_vs_Mtrun \
  --output-prefix Sjap_vs_Mtrun_circos
```

脚本默认会渲染 PNG/SVG。如果只想生成配置和审计文件，不运行 Circos，加：

```bash
--no-render
```

## 输入格式

- `--gff1/--gff2`：WGDI 格式 gff，不是普通 GFF3。脚本使用第 1 列染色体、第 2 列基因 ID、第 3/4 列真实 bp 坐标、第 6 列 WGDI 基因序号。
- `--len1/--len2`：WGDI len 文件。脚本使用第 1 列染色体名、第 2 列真实染色体长度；染色体顺序严格按 len 文件顺序。
- `--block-info`：WGDI blockinfo 或 `block_information.csv`。必须包含 `chr1`、`chr2`、`block1`、`block2`；建议包含 `id`、`start1`、`end1`、`start2`、`end2`、`length`。

脚本不解析 `.collinearity.txt` fallback；没有 blockinfo 会直接报错。

## 关键参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--min-block-anchors` | `50` | 主图只绘制 mapped anchor 数不少于该值的 block。默认取 50 是为了得到更宏观、更简洁的共线性结构。所有可映射 block 仍写入 `all_block_links.tsv`。 |
| `--max-links` | `0` | 在 `--min-block-anchors` 过滤之后，如果仍然太多，可只保留最大的 N 条 block ribbon。`0` 表示不限制。 |
| `--exclude-regex` | `scaffold|contig|unplaced|unlocalized|random|chrUn|tig` | 默认去掉 contig/scaffold 等碎片序列。 |
| `--keep-regex` | 无 | 如果设置，先只保留匹配该正则的序列，再应用 `--exclude-regex`。 |
| `--no-render` | 关闭 | 只写输入文件和报告，不运行 Circos。 |
| `--link-opacity` | `0.22` | ribbon 透明度。 |
| `--chrom-spacing` / `--spacing` | `0.006r` | 相邻染色体 sector 的 Circos 间距。更小更紧凑，例如 `0.004r`；更大更疏，例如 `0.012r`。 |
| `--bezier-radius` | `0.08r` | block ribbon 的 Bezier 控制半径。较小值可减少中心宽带状叠加；较大值会让 ribbon 更容易形成宽的中心带。 |
| `--link-radius` | `0.850r` | ribbon 起止端所在半径。通常不需要改。 |

`--min-block-anchors` 和 `--max-links` 的顺序固定：先按 anchor 数过滤，再按 `--max-links` 保留最大的 block。最大的定义为：mapped anchor 数更多优先；如相同，映射 bp span 更长优先；再相同，保留原始 blockinfo 中更靠前的 block。

如果看到 link 在图中心像是变粗，通常不是脚本动态改变线宽，而是 `ribbon = yes` 的区间 ribbon 在 Bezier 曲线靠近圆心时发生张开和透明叠加。默认 `--bezier-radius 0.08r` 会减轻中心宽带，同时仍保留 block 区间宽度信息。

## 输出文件

每个输出目录包含：

- `*_circos.png`、`*_circos.svg`：渲染结果。
- `karyotype.txt`：Circos 染色体长度和顺序，长度为真实 bp 坐标。
- `links.txt` / `links.full.txt`：主图实际绘制的 block-level interval ribbons。
- `links.all_blocks.txt` / `all_block_links.tsv`：全量可映射 block interval links，用于审计。
- `block_links.tsv`：主图绘制的 block links 表格。
- `mapped_blocks.tsv`：block 级审计表，含原始 WGDI 基因序号范围和映射后真实 bp 范围。
- `anchor_links.tsv`：anchor gene pair 级审计表，不参与绘图。
- `chromosomes.tsv`、`gene_density.txt`、`chrom_bands.txt`、`circos.conf`、`colors.conf`、`housekeeping.conf`：Circos 输入。
- `validation_report.md`、`VALIDATION_SUMMARY.tsv`：验收报告。

## 当前 6 组数据批量示例

脚本没有内置批量入口。需要批量绘制时，在 shell 中逐个调用即可。下面示例把结果另存到 `circos_plot/results_v3_spacing006_bezier008`：

```bash
OUT=circos_plot/results_v3_spacing006_bezier008

conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 01_guohuai_vs_guohuai/Sjap.gff --len1 01_guohuai_vs_guohuai/Sjap.len --name1 SjapL --gff2 01_guohuai_vs_guohuai/Sjap.gff --len2 01_guohuai_vs_guohuai/Sjap.len --name2 SjapR --block-info 01_guohuai_vs_guohuai/Sjap-Sjap_block_information.csv --outdir "$OUT/01_Sjap_vs_Sjap" --label 01_Sjap_vs_Sjap --output-prefix SjapL_vs_SjapR_circos
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 02_jiyehuai_vs_jiyehuai/S.jap.o.gff --len1 02_jiyehuai_vs_jiyehuai/S.jap.o.len --name1 SjapoL --gff2 02_jiyehuai_vs_jiyehuai/S.jap.o.gff --len2 02_jiyehuai_vs_jiyehuai/S.jap.o.len --name2 SjapoR --block-info 02_jiyehuai_vs_jiyehuai/Sjapo-Sjapo_blcok_information.csv --outdir "$OUT/02_Sjapo_vs_Sjapo" --label 02_Sjapo_vs_Sjapo --output-prefix SjapoL_vs_SjapoR_circos
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 03_xianghuai_vs_xianghuai/Cwilsonii.gff --len1 03_xianghuai_vs_xianghuai/Cwilsonii.len --name1 CwilsoniiL --gff2 03_xianghuai_vs_xianghuai/Cwilsonii.gff --len2 03_xianghuai_vs_xianghuai/Cwilsonii.len --name2 CwilsoniiR --block-info 03_xianghuai_vs_xianghuai/Cwilsonii-Cwilsonii_block_information.csv --outdir "$OUT/03_Cwilsonii_vs_Cwilsonii" --label 03_Cwilsonii_vs_Cwilsonii --output-prefix CwilsoniiL_vs_CwilsoniiR_circos
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 04_cjxh_vs_cjxh/Pplatycarpum.gff --len1 04_cjxh_vs_cjxh/Pplatycarpum.len --name1 PplatycarpumL --gff2 04_cjxh_vs_cjxh/Pplatycarpum.gff --len2 04_cjxh_vs_cjxh/Pplatycarpum.len --name2 PplatycarpumR --block-info 04_cjxh_vs_cjxh/Pplatycarpum-Pplatycarpum_block_information.csv --outdir "$OUT/04_Pplatycarpum_vs_Pplatycarpum" --label 04_Pplatycarpum_vs_Pplatycarpum --output-prefix PplatycarpumL_vs_PplatycarpumR_circos
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 05_guohuai_vs_jiyehuai/jiyehuai.gff --len1 05_guohuai_vs_jiyehuai/jiyehuai.len --name1 Sjapo --gff2 05_guohuai_vs_jiyehuai/Sjap.gff --len2 05_guohuai_vs_jiyehuai/Sjap.len --name2 Sjap --block-info 05_guohuai_vs_jiyehuai/Sjapo-Sjap.blockinfo.csv --outdir "$OUT/05_Sjapo_vs_Sjap" --label 05_Sjapo_vs_Sjap --output-prefix Sjapo_vs_Sjap_circos
conda run -n circos python circos_plot/wgdi_blockinfo_to_circos.py --gff1 06_guohuai_vs_Medicago/Sjap.gff --len1 06_guohuai_vs_Medicago/Sjap.len --name1 Sjap --gff2 06_guohuai_vs_Medicago/Mtrun.gff --len2 06_guohuai_vs_Medicago/Mtrun.len --name2 Mtrun --block-info 06_guohuai_vs_Medicago/Sjap-Mtrun_block_information.csv --outdir "$OUT/06_Sjap_vs_Mtrun" --label 06_Sjap_vs_Mtrun --output-prefix Sjap_vs_Mtrun_circos
```

## 验收

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

预期：`01-05` 为 28 个 sector，`06` 为 22 个 sector；`validation_report.md` 中 `coordinate errors` 必须为 0。默认 `--min-block-anchors 50` 下，主图 link 数应约为 80-120 条，明显比上一版更宏观。
