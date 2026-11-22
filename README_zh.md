![DEGAP](./DEGAP.png)

# DEGAP 2.0 — 基因组组装路径动态延伸

[English](README.md) | **中文**

**DEGAP 2.0** 利用 HiFi 和/或 Oxford Nanopore（ONT）长读段数据，完成基因组组装缺口填补、contig 连接以及染色体末端延伸。

## DEGAP 2.0 新增特性

相较于 DEGAP 1.0（[v1 历史指南（英文）](v1/README.md)），DEGAP 2.0 带来四项核心升级：

1. **更高效的候选 reads 筛选** — DEGAP 1.0 已提供比对后的测序深度过滤；DEGAP 2.0 新增可选的 `--kmer_filter`（默认关闭），可在比对*之前*收窄候选 reads 范围，从而在大数据集上减少比对计算量。
2. **双向缺口延伸** — 将单侧缺口填补升级为顺序进行的左右交替延伸：每次生成的新累计版本都会与对侧全部历史版本（包括原始两侧序列）比对检查闭合；若一侧正常停止，其历史版本会冻结并继续作为对侧的匹配目标。
3. **TelSeeker 染色体末端延伸** — 新增针对指定染色体末端向端粒 reads 进行单端延伸的 `telseeker` 模式，以及面向单条 contig 末端的 `telseeker_ctg` 变体，作为内部缺口填补与 contig 连接的补充。
4. **更清晰的离线交互式报告** — 生成自包含的离线交互式 HTML 报告（并非实时监控面板），用于说明任务结果、延伸历史、所选连接与证据以及停止原因；在本地浏览器中直接打开生成的 `.html` 文件即可查看（详见[报告阅读指南](HOWTO_zh.md#查看-html-报告)；`TelSeeker Ctg` 目前通过 `result/status.tsv` 输出结果而非 HTML 报告）。

## 安装

DEGAP 2.0 运行于 Linux。仓库提供的 [environment.yml](environment.yml) 固定 Python 3.10，并列出所需依赖包与外部工具：

- **Python 依赖包**：Python 3.10、Biopython、Pysam、NumPy、Matplotlib、SciPy
- **外部工具**：Minimap2 (2.28)、Hifiasm (0.25.0)、SAMtools (1.22.1)、Seqkit (2.8.0)、MUMmer4、GNU Parallel、Ripgrep (15.1.0)

推荐使用 [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) 创建运行环境：

```bash
micromamba env create -f environment.yml
micromamba activate degap_env
```

当 `environment.yml` 更新后，可同步更新已有环境：

```bash
micromamba env update -f environment.yml --prune
```

如需不使用 Micromamba 手动安装，请参见 [HOWTO_zh.md — 安装](HOWTO_zh.md#安装)。

## 模式选择

| 模式 | 命令入口 | 用途 | 核心输入 |
| --- | --- | --- | --- |
| **GapFiller** | `python bin/DEGAP.py --mode gapfiller` | 基于双向延伸家族填补已知左右两侧序列（flank）之间的单个缺口 | `--seqleft`、`--seqright`、`--hifi`/`--ont`、`-o/--out` |
| **CtgLinker** | `python bin/DEGAP.py --mode ctglinker` | 使用主家族及正向/反向候选家族，从 contig 边缘连接 contig 并填补缺口 | `--ctgseq`、`--hifi`/`--ont`、`-o/--out` |
| **TelSeeker** | `python bin/DEGAP.py --mode telseeker` | 将指定的染色体末端（`<chr>.L` / `<chr>.R`）向端粒 reads 进行单端延伸 | `--genome`、`-m/--motif`、`-e/--target_ends`、`--hifi`/`--ont`、`-o/--out` |
| **TelSeeker Ctg** | `python bin/DEGAP.py --mode telseeker_ctg` | 对一个或多个单条记录的 contig 末端（`L` 或 `R`）进行单端延伸，并共享 reads 与端粒缓存 | `--ctg <fasta> <L\|R>`、`-m/--motif`、`--hifi`/`--ont`、`-o/--out` |
| **AutoGapfiller** | `python bin/AutoGapfiller.py` | 自动检测基因组中的 N 缺口，并为每个缺口运行一个双向缺口填补任务 | `--genome`、`--hifi`/`--ont`、`-o/--output`（可选 `--setup_only`） |

运行 TelSeeker 前，可先使用 `python bin/TelSeekerCheck.py --genome <genome.fa> --motif <MOTIF> --out <dir>` 生成染色体末端人工核查图与 2 kb 末端 FASTA 文件。

## 快速开始

填补已知左右两侧序列（flank）之间的单个缺口（`gapfiller`）：

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./left_flank.fa \
    --seqright ./right_flank.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./gap_out \
    --MaximumExtensionRound 25 \
    -t 20
```

任务完成后：

- 在本地浏览器中直接打开 `./gap_out/family-report.html`，查看闭合状态、延伸家族历史、候选连接、连接位点指标和停止原因。
- 使用 `./gap_out/<name>.final.fa` 前，请先在 `./gap_out/family-result.json`（或报告顶部状态标签）中确认 `status: "closed"`；进程以退出码 0 结束并不等同于缺口已闭合，也不能直接证明序列在生物学上完全正确。
- 凡涉及纯 ONT 组装或混合模式回退到纯 ONT 组装的场景，`--ont` 必须提供**原始 FASTQ/FASTQ.gz** 文件（`--hifi` 同时支持 FASTA 与 FASTQ）。

其余模式用法、共享预算规则、断点恢复策略与完整参数，请参阅 [HOWTO_zh.md](HOWTO_zh.md)。

## 文档导航

- **实用用户指南**：[HOWTO_zh.md](HOWTO_zh.md) — 全部五种模式的完整操作参考、参数表、HTML 报告解读、输出目录结构与常见问题排查。
- **设计说明**：
  - [docs/design/bidirectional-extension.md](docs/design/bidirectional-extension.md) — 双向历史延伸家族、共享预算与断点恢复规则。
  - [docs/design/junction-validation.md](docs/design/junction-validation.md)（英文）— 全序列回退比对、原始两侧序列（flank）核心保护、有 reads 支持的有界末端修订及锚点回退重试机制。
  - [docs/design/telomere-target-search.md](docs/design/telomere-target-search.md)（英文）— TelSeeker 共享目标索引、分批候选验证与局部接头佐证策略。
  - [docs/design/report-presentation.md](docs/design/report-presentation.md)（英文）— 离线 HTML 报告呈现与结构对齐说明。
- **阶段性验证记录（均为英文）**：
  - [SSH 180 扩展真实数据验证矩阵](docs/validation/2026-11-22-ssh180-expanded-matrix.md) — 涵盖固定缺口、多工作流模式、共享限额及负样本的 39 组工作负载。
  - [SSH 180：两个遮蔽 HiFi 缺口闭合验证](docs/validation/2026-11-22-ssh180-new-dataset.md) — 纯缺口序列实测对比与 ONT 适用边界。
  - [SSH 180 连接验证与末端修订测试](docs/validation/2026-11-22-junction-repair.md) — 连接检查与有界末端修订验证记录。
  - [SSH 180 TelSeeker 目标检索验证](docs/validation/2026-11-22-telomere-target-search.md) — 共享目标索引、分批比对验证与局部支持度检查。
  - 早期验证记录：[docs/validation/2026-11-21-ssh180-real-ont.md](docs/validation/2026-11-21-ssh180-real-ont.md)、[docs/validation/2026-11-21-ont-quality-repair.md](docs/validation/2026-11-21-ont-quality-repair.md)。
- **历史 v1 指南**：[v1/README.md](v1/README.md)（英文）— 旧版 DEGAP 1.0 文档（仅作历史参考）。

## 引用

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. *Briefings in Bioinformatics*, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

## 联系方式

如有意见、问题反馈或使用咨询，请联系 <jzhang@mail.hzau.edu.cn>。
