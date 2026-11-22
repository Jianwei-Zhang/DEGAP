![DEGAP](./DEGAP.png)

# DEGAP 2.0 — 用户指南

[English](HOWTO.md) | **中文**

**DEGAP 2.0** 的完整操作参考手册，涵盖运行、断点恢复、报告核查与常见问题排查。

## 目录

- [安装](#安装)
- [输入 Reads 与 ONT 质量契约](#输入-reads-与-ont-质量契约)
- [延伸模型、预算、错误处理与断点恢复](#延伸模型预算错误处理与断点恢复)
- [运行模式与示例](#运行模式与示例)
  - [1. GapFiller](#1-gapfiller)
  - [2. CtgLinker](#2-ctglinker)
  - [3. TelSeeker 与 TelSeekerCheck](#3-telseeker-与-telseekercheck)
  - [4. TelSeeker Ctg](#4-telseeker-ctg)
  - [5. AutoGapfiller](#5-autogapfiller)
- [参数参考](#参数参考)
- [查看 HTML 报告](#查看-html-报告)
- [输出目录结构](#输出目录结构)
- [常见问题排查](#常见问题排查)
- [相关文档](#相关文档)
- [引用与联系方式](#引用与联系方式)

---

## 安装

推荐使用 Linux 下的 Micromamba 环境以及 [environment.yml](environment.yml) 中固定的版本配置，具体命令见 [README_zh.md — 安装](README_zh.md#安装)。

### 手动安装（不使用 Micromamba）

安装 Python 依赖包（推荐 Python 3.10）：

```bash
pip install biopython pysam numpy matplotlib scipy
```

确保以下外部工具已加入 `PATH`（版本建议与 [environment.yml](environment.yml) 保持一致）：

- **Minimap2**：<https://github.com/lh3/minimap2>
- **Hifiasm**：<https://github.com/chhylp123/hifiasm>
- **SAMtools**：<http://www.htslib.org/>
- **Seqkit**：<https://github.com/shenwei356/seqkit>
- **MUMmer4**：<https://github.com/mummer4/mummer>
- **GNU Parallel**：<https://www.gnu.org/software/parallel/>
- **Ripgrep**：<https://github.com/BurntSushi/ripgrep>

---

## 输入 Reads 与 ONT 质量契约

必须至少提供 `--hifi` 或 `--ont` 之一。两个参数均支持传入多个文件（例如 `--hifi HIFI1.fq.gz HIFI2.fq.gz`）。

| 输入模式 | 支持的文件格式 | 局部组装要求 |
| --- | --- | --- |
| **纯 HiFi**（`--hifi`） | `.fa`、`.fasta`、`.fna`、`.fq`、`.fastq`（支持 `.gz`） | HiFi 的 FASTA 和 FASTQ 均受支持。 |
| **纯 ONT**（`--ont`） | **必须提供原始 `.fq` / `.fastq`（支持 `.gz`）** | `hifiasm --ont` 需要真实的 Phred+33 质量值；仅提供 FASTA 的 ONT 输入会在组装阶段被明确拒绝。 |
| **混合模式**（`--hifi` + `--ont`） | HiFi 支持 FASTA/FASTQ；ONT 推荐使用 FASTQ | 仅当本轮仍有可用 HiFi reads 时，ONT FASTA 才可作为补充的 `--ul` 输入使用。一旦某轮回退至纯 ONT 组装，就必须具备原始 ONT FASTQ 质量值。 |

ONT 质量恢复机制说明（`bin/ONTAssemblyReads.py`）：

- 尽管预处理阶段会将 reads 规范化为工作 FASTA，且命令行解析允许传入 FASTA 路径，但纯 ONT 局部组装**绝不会伪造合成质量值**。
- 原始 ONT FASTQ/FASTQ.gz 文件（单个或多个）会按需延迟索引为压缩 SQLite 查询库（`ont-quality-cache/`），并在单次搜索的各延伸家族与轮次间共享。
- 筛选出的 reads 会按所选方向还原为 FASTQ（反向互补的 read 会同步反转质量值），并执行严格的 read ID、序列内容与唯一裁剪位置校验；遇到缺失 ID、重复 ID、序列不一致、裁剪位置歧义或格式损坏的 FASTQ 输入时，程序会明确报错退出。
- 原始 ONT 文件（包括其中的质量值）属于断点校验契约的一部分；即使仅修改了质量值，也会导致已有断点无法复用。

---

## 延伸模型、预算、错误处理与断点恢复

### 延伸家族搜索模型

- **GapFiller / AutoGapfiller（双向、顺序交替）**：左侧（`L0, L1, ...`）与右侧（`R0, R1, ...`）两条延伸链按轮次左右交替延伸（并非两侧同时并行执行）。每个版本均为累计序列。先执行第 0 轮直接闭合检查（`L0` 对 `R0`）；此后每产生一个新版本，都会与对侧**全部**历史版本（包含 `L0`/`R0`）进行闭合比对。若一侧因无可用延伸 reads/contig 或零增长等原因正常停止，其历史版本会冻结并继续作为对侧的匹配目标，另一侧仍可继续延伸。
- **CtgLinker（主家族 + 正向/反向候选家族）**：每次外端连接搜索都会为当前活跃的 contig/scaffold 建立主家族，并为其余所有可用 contig/scaffold 分别建立正向（`+`）和反向互补（`-`）候选家族。主家族与各候选家族在该次端点搜索内交替推进，并共享同一份轮数与延伸长度预算。同一条原始输入 contig 不会被重复拼入。
- **TelSeeker / TelSeeker Ctg（单端延伸）**：复用相同的持久化单轮延伸引擎，但始终保持**单端**模式——仅维护一条染色体侧延伸家族，不虚构对侧家族。连接判定以匹配共享目标索引中的静态端粒 reads 为准，而非与初始 seed 发生重叠。
- **确定性择优规则**：完成当前轮次的比对前沿后，所有合格连接按以下优先级依次排序：
  1. 比对一致性（identity，从高到低）
  2. 有效重叠长度（从长到短）
  3. 末端悬垂长度（overhang，从短到长）
  4. 固定的家族/版本/坐标顺序

  最终仅采用并输出唯一最优的合格连接。需注意，在重复区域中，比对一致性最高并不直接等同于生物学正确，仍需结合导出的证据与比对结果核查。

完整的调度与状态机规范见 [docs/design/bidirectional-extension.md](docs/design/bidirectional-extension.md)。

### 共享预算

- `--MaximumExtensionRound` 统计**单次搜索内所有延伸家族累计尝试的延伸轮数**（GapFiller/AutoGapfiller 为左右两侧之和；CtgLinker 为单次端点搜索中主家族与全部候选家族之和）。正常失败的延伸尝试同样计入轮数。
- `--MaximumExtensionLength` 累加该次搜索中所有延伸家族的新增碱基数。原始两侧序列（flank）和历史快照本身不计入，最终闭合重叠部分也不扣减；即使最终选用了较早的历史版本完成闭合，已消耗的轮数与新增碱基数也不会退还。
- 长度上限在启动下一轮尝试**之前**检查，因此最后一轮完整执行的延伸可能会超出 `--MaximumExtensionLength`；该轮结果仍会完整保留并参与闭合检查，之后再停止。
- 设置 `--MaximumExtensionRound 0` 会跳过所有延伸尝试，但在支持第 0 轮直接检查的模式中仍允许直接闭合检查。负数上限会被入口直接拒绝。两个参数默认均为 `None`（不设上限）。

### 连接验证与两侧序列保护

- 若末端窗口比对未找到合格连接，程序会自动回退到全序列比对；满足条件的小段共线比对块可组合为一条经验证的连接，且不会为了强行闭合短重叠案例而静默降低最小有效重叠门槛。
- 原始两侧序列（flank）的核心区域受保护；有界末端修订**并非**无条件抛光。仅允许在面向缺口的每侧末端最多 2 kb 范围内修订，且必须同时满足原始跨越 reads 支持（包含针对 seed 继承错误的原始 reads 共识候选）、双链支持以及原始边界的唯一锚定定位。未获充分支持的差异仅作记录而不应用。
- 纯缺口长度（gap-only length）、两侧序列长度变化（flank length delta）与原始真实重叠裁剪（true-overlap trimming）作为独立指标分别报告。
- 初始端点的锚点向内回退重试属于一次额外调度的尝试，计入共享轮数预算。

详细默认阈值、歧义处理、AutoGapfiller 最小替换片段及证据文件结构，请参见 [docs/design/junction-validation.md](docs/design/junction-validation.md)（英文）。

### 运行状态与工具错误

请将**程序正常结束（退出码 0）**、**缺口闭合/末端连接成功**与**序列通过生物学验证**区分为三种独立状态：

- 即使 `hifiasm` 返回退出码 `0`，只要其标准输出/标准错误中出现显式错误信息，DEGAP 仍会按工具错误终止。
- 缺少主 contig GFA 文件（`.p_ctg.gfa`）按工具错误处理，不会被误判为生物学不可延伸。
- 若 `hifiasm` 正常完成并写出合法的空 GFA（0 条序列记录）且未产生可用延伸 contig，这属于正常停止（`no_extension_contigs`），与工具错误严格区分。
- 工具错误绝不会写入已完成的组装阶段回执，也不会发布为正常的 `exhausted` 结果；此前已完成的阶段（如 `mapping/` 和 `reads/`）在修复错误后仍可直接复用。

### 断点恢复、兼容性、`--flag` 与 `--remove`

- **断点恢复**：保持命令与输入文件完全不变后重新运行（或附加 `--resume_auto`），即可基于内容校验的延伸家族断点（`families.json`、`ctg-families.json` 及各阶段回执）自动续跑。中断时未提交的临时阶段会归档到 `.incomplete.*` 目录中，不会被直接信任。
- **不支持的旧版恢复参数**：`DEGAP.py` 的全部四种模式均拒绝 `--resume ROUND`（不再支持任意轮次手动回退）。
- **旧版目录与算法契约变更**：旧版单侧输出目录，以及在旧算法契约下生成的目录，均不会被静默迁移或覆盖。请保留原有证据目录不动，并使用**新的输出目录**重新运行。
- **`--flag`**：在 `gapfiller` 模式中已废弃，仅作兼容保留。传入该参数会打印提示信息并直接忽略，新命令中无需再写。
- **`--remove`**：旧版清理级别参数（`1`、`2`、`3`；默认 `2`）。在延伸家族工作流中，因为断点恢复依赖这些文件，预处理后的 reads、分块、索引、历史版本及已验证阶段在 `--remove 1` 和 `--remove 2` 下都会**始终保留**。确认不再需要历史记录与续跑后，请手动清理已结束的工作目录。

---

## 运行模式与示例

### 1. GapFiller

填补已知左右两侧序列（flank）之间的单个缺口。

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./gap_left.fa \
    --seqright ./gap_right.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_gapfiller \
    --MaximumExtensionRound 25 \
    -t 20
```

纯 HiFi 模式，同时启用可选的比对后深度过滤（`--filterDepthHifi`）与比对前 k-mer 预筛选（`--kmer_filter`）：

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./gap_left.fa \
    --seqright ./gap_right.fa \
    --hifi ./hifi_reads.fa \
    -o ./output_gapfiller_hifi \
    --filterDepthHifi 0.3 \
    --kmer_filter \
    -t 20
```

请检查 `family-result.json` 或 `family-report.html`：只有 `status: "closed"` 表示缺口已闭合为单条序列；否则 `<name>.final.fa` 中导出的是彼此独立的未闭合延伸家族序列。

### 2. CtgLinker

从 contig 边缘出发连接 contig 并填补缺口。

```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./contigs.fa \
    --hifi ./HIFI1.fq.gz ./HIFI2.fq.gz \
    --ont ./ONT1.fq.gz ./ONT2.fq.gz \
    --out ./output_ctglinker \
    --MaximumExtensionRound 25 \
    -t 20
```

纯 ONT 示例（必须使用原始 FASTQ/FASTQ.gz）：

```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./contigs.fa \
    --ont ./ont_reads.fq.gz \
    --out ./output_ctglinker_ont \
    -t 20
```

`DG.scaffold.agp` 基于 `DG.scaffold.components.fa` 精确重建 `DG.scaffold.fa`（其中包含实际延伸快照以及缺口/修订溯源组件，而不仅是未改动的原始 contig）。

### 3. TelSeeker 与 TelSeekerCheck

#### 步骤 A（可选）：生成染色体末端人工核查文件

```bash
python bin/TelSeekerCheck.py \
    --genome ./genome.fa \
    --motif TTAGGG \
    --out ./telomere_check
```

`TelSeekerCheck.py` 会在 `<out>/genome.telomere.check/` 中输出人工核查文件：

- `genome.telomere.check.left.2kb.fa` 与 `genome.telomere.check.right.2kb.fa`（默认提取长度 `--extract-length 2000`）
- `all_chromosomes_combined.png` 与 `<chromosome>_telomere_motif.png`（默认 1 kb 绘图窗口，`--plot-window-size 1000`）

该脚本**不会**自动输出分类清单（不生成 `genome.telomere.check.csv`、`need_extension_chr_end.txt` 或 `uncertain_chr_end.txt`）。请结合图示与末端序列人工判读，然后列出需要延伸的染色体末端（每行一个 `<chromosome_id>.L` 或 `<chromosome_id>.R`，允许空行和 `#` 注释）：

```text
# 根据 TelSeekerCheck 核查图选定的目标染色体末端
Chr01.L
Chr03.R
```

#### 步骤 B：运行 TelSeeker

既可直接传入目标末端（`-e Chr01.L Chr03.R`），也可通过文本文件传入（`-e target_ends.txt`）：

```bash
python bin/DEGAP.py --mode telseeker \
    --genome ./genome.fa \
    --motif TTAGGG \
    -e ./target_ends.txt \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_telseeker \
    --work 4 \
    --tel-n 100 \
    --tel-r 0.6 \
    --tel-mm 0 \
    --MaximumExtensionRound 25 \
    -t 20
```

常见的大写 A/C/G/T 端粒重复基序：`TTAGGG`（脊椎动物）、`TTTAGGG`（植物）、`TTAGG`（昆虫）。

- **末端端粒 reads 识别**：第 1 步在每条 read 的两个物理末端分别扫描 `--tel-n * len(--motif)` bp，使用基序的所有循环移位及反向互补循环移位进行匹配。当 `hits / --tel-n >= --tel-r` 时该 read 末端通过筛选，每个基序长度单元允许 `--tel-mm`（`0` 或 `1`）个错配。以 `TTAGGG` 为例，`--tel-n 100` 扫描 600 bp，`--tel-r 0.6` 要求至少命中 60 次。若第 1 步筛选到的端粒 reads 过少，可适当降低 `--tel-r`、设置 `--tel-mm 1` 或减小 `--tel-n`。
- **单端共享目标检索**：TelSeeker 为每个染色体末端仅维护一条染色体侧延伸家族，以共享的静态 SQLite 目标索引检索端粒 reads。所有检索到的候选记录均按限定批次规模进行多查询比对验证（并非 top-N 截断），经确定性排序后，仅当在非重复锚点及接头向外最多 **200 bp** 范围内获得**至少 2 个不同 read ID** 的佐证时才接受连接。支持度计数属于局部接头的下界，**不代表**整条新增端粒尾部的共识纠错或抛光。最终拼接时完整保留当前染色体侧序列，仅在其目标端接入所选 read 的外侧尾部序列（`L` 端前接、`R` 端后接）。具体阈值、缓存与采样边界见 [docs/design/telomere-target-search.md](docs/design/telomere-target-search.md)（英文）。
- **人工复核**：即使已生成 linker 并整合至 `part3.integration.results/final.genome.fa`，也不能直接保证恢复了完整的生物学端粒；请务必结合 `part3.integration.results/final.genome.telomere.check/` 中的末端序列与图示进行人工确认。

### 4. TelSeeker Ctg

对一个或多个独立的单记录 contig FASTA 末端（`L` 或 `R`）进行单端延伸，并共享预处理后的 reads 与 `part1.telo.reads/`。

```bash
python bin/DEGAP.py --mode telseeker_ctg \
    --ctg ./chr1_partial.fa L \
    --ctg ./chr2_partial.fa R \
    --motif TTAGGG \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_telseeker_ctg \
    --work 4 \
    --tel-n 100 \
    --tel-r 0.6 \
    --tel-mm 0 \
    --MaximumExtensionRound 25 \
    --MaximumExtensionLength 1000000 \
    -t 20
```

`telseeker_ctg` 使用约束：

- 每个 `--ctg <fasta> <L|R>` 指定的 FASTA 文件必须且只能包含**一条**序列记录。
- 仅当 `--hifi`、`--ont`、`--motif`、`--tel-n`、`--tel-r` 与 `--tel-mm` 与 `telseeker_ctg.input/cache_manifest.json` 完全一致时，才允许复用同一个 `--out` 目录追加新的 `--ctg` 任务；否则 DEGAP 会报缓存不匹配错误并停止。
- `telseeker_ctg` 不支持 `--filterDepthHifi` 与 `--filterDepthOnt`。
- 当前公开命令行入口**不会**生成 HTML 报告；请查看汇总表 `result/status.tsv` 以及各任务的 `telseeker_ctg.jobs/<safe_ctg_id>.<L|R>/result/`。缺少匹配端粒 reads 的任务会明确以 `no_tel_reads` 或 `no_tel_reads_for_target_end` 失败。

### 5. AutoGapfiller

自动检测基因组草图中的全部 N 缺口（长度 `>= 20 bp`，间距 `< 11,000 bp` 的相邻缺口自动合并），为**每个缺口创建一个双向延伸任务**，并行执行任务、整合已闭合缺口，并生成离线 HTML 报告。注意 `AutoGapfiller.py` 使用的输出目录参数为 `--output` / `-o`（不是 `--out`）。

```bash
python bin/AutoGapfiller.py \
    --genome ./genome_with_gaps.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --output ./output_autogap \
    --work 4 \
    --MaximumExtensionRound 25 \
    --MaximumExtensionLength 500000 \
    --kmer_filter \
    -t 20
```

添加 `--setup_only` 可仅执行步骤 1–3.1（reads 预处理、缺口检测及任务目录 / `run.sh` 生成），暂不启动实际缺口填补：

```bash
python bin/AutoGapfiller.py \
    --genome ./genome_with_gaps.fa \
    --hifi ./hifi_reads.fq.gz \
    --output ./output_autogap \
    --setup_only
```

---

## 参数参考

### `bin/DEGAP.py`

| 参数 | 说明 | 默认值 | 适用模式 |
| --- | --- | --- | --- |
| `--mode` | 运行模式：`gapfiller`、`ctglinker`、`telseeker` 或 `telseeker_ctg` | 必填 | 全部 |
| `-o, --out` | 输出目录 | 必填 | 全部 |
| `--hifi` | 一个或多个 HiFi reads 文件（FASTA/FASTQ，可选 `.gz`） | `None` | 全部（与 `--ont` 至少提供一个） |
| `--ont` | 一个或多个 ONT reads 文件（纯 ONT 组装及混合模式回退纯 ONT 时必须提供原始 FASTQ/FASTQ.gz） | `None` | 全部（与 `--hifi` 至少提供一个） |
| `-t, --thread` | 线程数 | `20` | 全部 |
| `--edge` | 边界与错拼检查的边缘容差（bp） | `500` | 全部 |
| `--MaximumExtensionRound` | 单次搜索的共享最大尝试轮数（`0` 仅允许第 0 轮直接闭合检查） | `None` | 全部 |
| `--MaximumExtensionLength` | 单次搜索的共享最大新增碱基数（bp，启动下一轮前检查） | `None` | 全部 |
| `--kmer_filter` | 启用比对前可选 k-mer 候选 reads 预筛选（`store_true`） | `False` | 全部 |
| `-ks, --kmer_size` | 启用 `--kmer_filter` 时的 k-mer 长度 | `41` | 全部 |
| `-kn, --kmer_num` | 启用 `--kmer_filter` 时的 k-mer 数量 | `20` | 全部 |
| `--filterDepthHifi` | HiFi 比对后深度过滤比例（如 `0.3`）；`gapfiller` 实现，`ctglinker`/`telseeker` 跳过，`telseeker_ctg` 拒绝 | `None` | `gapfiller` |
| `--filterDepthOnt` | ONT 比对后深度过滤比例；各模式支持情况与 `--filterDepthHifi` 相同 | `None` | `gapfiller` |
| `--remove` | 旧版清理级别（`1`、`2`、`3`）；延伸家族历史、预处理 reads、索引和已验证阶段仍会保留 | `2` | 全部 |
| `--resume_auto` | 从已验证的延伸家族断点恢复运行（与直接重跑相同命令等价） | `False` | 全部 |
| `--resume` | 已废弃的旧版单侧指定轮次恢复（传入会报错拒绝） | `None` | 不支持 |
| `--seqleft`, `--seqright` | 左右两侧序列（flank）FASTA 文件 | 必填 | `gapfiller` |
| `--flag` | 已废弃；传入会打印提示并忽略 | `None` | `gapfiller` |
| `--ctgseq` | 输入 contig FASTA 文件 | 必填 | `ctglinker` |
| `--genome` | 输入基因组/染色体 FASTA 文件（`--chr` 已废弃） | 必填 | `telseeker` |
| `-e, --target_ends` | 目标染色体末端列表（`Chr01.L Chr01.R`）或每行一个 `<chr>.<L\|R>` 的文本文件 | 必填 | `telseeker` |
| `--ctg` | 可重复传入的 `<contig.fa> <L\|R>` 单记录 contig 末端任务 | 必填 | `telseeker_ctg` |
| `-m, --motif` | 大写 A/C/G/T 端粒重复基序 | 必填 | `telseeker`, `telseeker_ctg` |
| `-tn, --tel-n` | 每条 read 末端扫描的基序长度单元数 | `100` | `telseeker`, `telseeker_ctg` |
| `-tr, --tel-r` | 末端命中比例下限（`hits / tel-n`），取值范围 `(0, 1]` | `0.6` | `telseeker`, `telseeker_ctg` |
| `-tmm, --tel-mm` | 每个基序长度单元允许的错配数（`0` 或 `1`） | `0` | `telseeker`, `telseeker_ctg` |
| `-w, --work` | 并行处理的染色体/contig 末端任务数 | `1` | `telseeker`, `telseeker_ctg` |

### `bin/AutoGapfiller.py`

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `-g, --genome` | 含 N 缺口的基因组草图 FASTA | 必填 |
| `-o, --output` | 输出目录（注意参数名为 `--output`，不是 `--out`） | 必填 |
| `--hifi` | 一个或多个 HiFi reads 文件（FASTA/FASTQ，可选 `.gz`） | `None` |
| `--ont` | 一个或多个 ONT reads 文件（纯 ONT 组装及混合模式回退纯 ONT 时必须提供原始 FASTQ/FASTQ.gz） | `None` |
| `-w, --work` | 并行运行的缺口填补任务数 | `1` |
| `--setup_only` | 仅执行预处理、缺口检测与任务脚本准备，不运行缺口填补 | `False` |
| `-t, --thread` | 传递给每个 DEGAP 子任务的线程数 | `20` |
| `--kmer_filter` | 启用比对前可选 k-mer 预筛选（`store_true`） | `False` |
| `-ks, --kmer_size` | k-mer 长度 | `41` |
| `-kn, --kmer_num` | k-mer 数量 | `20` |
| `--MaximumExtensionLength` | 每个缺口的共享最大新增碱基数（bp） | `None` |
| `--MaximumExtensionRound` | 每个缺口的共享最大尝试轮数 | `None` |
| `--filterDepthHifi`, `--filterDepthOnt` | 第 1 步预处理阶段可选的比对后深度过滤阈值 | `None` |
| `--edge` | 边缘容差（bp） | `500` |
| `--remove` | 旧版清理级别（`1`、`2`、`3`）；延伸家族恢复所需文件仍会保留 | `2` |

### `bin/TelSeekerCheck.py`

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `--genome` | 输入基因组 FASTA | 必填 |
| `--out` | 输出目录（生成 `<out>/genome.telomere.check/`） | 必填 |
| `-m, --motif` | 端粒重复基序 | `TTAGGG` |
| `--extract-length` | 提取至左右末端 FASTA 的末端序列长度（bp） | `2000` |
| `--plot-window-size` | 全基因组基序分布图的统计窗口大小（bp） | `1000` |
| `--check-length` | 末端检查区间长度（bp） | `100000` |

> **过滤机制说明**：DEGAP 1.0 已具备比对后深度过滤（旧版在 `--reads` 上使用 `--filterDepth` / `--filterDepthOnly`）。在 DEGAP 2.0 中，按测序平台区分的 `--filterDepthHifi` 与 `--filterDepthOnt`（用于 `DEGAP.py --mode gapfiller` 与 `AutoGapfiller.py` 第 1 步）在比对后通过 `selectRawReads.py` 执行过滤，保留覆盖深度落在 `(threshold * mode_depth, (2 - threshold) * mode_depth)` 之外的 reads、靠近序列末端的 reads、未比对 reads 以及比对覆盖比例 `< 99%` 的 reads，并写出 `Genome.inputCtg.usefullReads.fasta`（混合模式同时写出 `Genome.inputCtg.usefullReads.ont.fasta`）。与之相对，`--kmer_filter`（默认 `False`，可通过 `-ks`/`-kn` 调节）是可选的比对前 k-mer 预筛选，用于在映射*之前*缩小候选 reads 集合；其耗时与灵敏度表现取决于具体数据集。

---

## 查看 HTML 报告

DEGAP 2.0 生成的 HTML 报告均为**离线、自包含的运行后汇总报告**（无需联网或启动服务，**并非**实时监控面板）。请在保持输出目录结构完整的前提下，直接在本地浏览器中打开入口 `.html` 文件，以便文中的相对路径产物链接能够正常跳转。

HTML 报告仅在任务执行到报告发布阶段后才会写出；如果运行在早期因输入校验、契约冲突或外部工具报错而中止，可能不会生成任何 HTML 文件。同样，一份完整或显示已闭合的 HTML 报告仅记录计算证据，本身并不等同于生物学正确性证明。

### 各模式报告与结果入口路径

| 模式 | 主入口路径 | 关联详情报告 / 结果文件 |
| --- | --- | --- |
| **GapFiller** | `<out>/family-report.html` | 直接链接 `family-result.json`、`families.json`、`<name>.final.fa`、`process.agp`、`connection-diagnostics.json`，闭合时还包含 `junction-evidence.json`。 |
| **CtgLinker** | `<out>/ctglinker-family-report.html` | 链接各次端点搜索的 `project/endpoint-XXXXXX/family-report.html`，以及 `DG.scaffold.fa`、`DG.scaffold.agp`、`DG.scaffold.components.fa` 和 `ctglinker-result.json`。 |
| **AutoGapfiller** | `<output>/05.visualization/Global.report.html` | 链接复制到同目录的逐缺口报告 `<output>/05.visualization/<chrom>.<gap_id>.bidirectional.html`（产物链接已重定向至 `03.gap_filling/`）。使用 `--setup_only` 时不生成报告。 |
| **TelSeeker** | `<out>/visual.report/Global.report.html` | 链接各染色体末端的 `<out>/visual.report/<Chr>.<L\|R>.report.html` 报告及目标检索证据文件。 |
| **TelSeeker Ctg** | *当前公开 CLI 不生成 HTML 报告* | 结果入口为 `<out>/result/status.tsv`，逐任务结果位于 `telseeker_ctg.jobs/<safe_ctg_id>.<L\|R>/result/status.tsv`、`extended.ctg.fa` 或 `failed.reason.txt`，以及可用的 `extension/telomere-result.json`。 |

### 报告重点核查项

1. **区分缺口闭合/末端连接与程序正常退出**
   - 退出码为 0 不能作为闭合或准确性的判断依据；在批量任务模式中，还应逐项核查每个子任务的状态。
   - **GapFiller**：重点查看顶部 `Closed / 已闭合` 状态标签（对应 `family-result.json` 中的 `status: "closed"`）。凡显示 `Unclosed` 时，左右延伸家族序列均分别导出，并未拼接闭合。
   - **AutoGapfiller**：在任务表和总览指标中查看 **Filled**（`Filled / 是否填补`）列，而不是只看是否执行过延伸；点击具体缺口行可进一步打开对应的 `.bidirectional.html` 家族报告。
   - **CtgLinker**：核查 **Endpoint searches** 表格中每次搜索的状态（`closed` 或未闭合）以及 **Scaffolds** 成员构成表；未闭合的端点搜索历史会保留链接供复核，但不会被强行拼入 scaffold。
   - **TelSeeker**：`connected` 仅表示存在通过局部佐证的端粒 read linker 并已整合至 `final.genome.fa`；最终端粒是否完整仍需在 `part3.integration.results/final.genome.telomere.check/` 中人工复核。
2. **区分延伸增长量与纯缺口/连接序列长度**
   - **Added bases / Total extension** 统计各轮次与各延伸家族在扣减重叠前的累计新增碱基数，切勿直接当作缺口长度解读。
   - 在已闭合的 **GapFiller / AutoGapfiller** 结果中，需分别查看独立的连接位点指标：**Gap-only length**（投影到原始两侧序列边界之间恢复出的纯缺口碱基数）、**Flank length delta**（`<= 2 kb` 末端有界修订带来的净长度变化，如有）以及 **Original overlap trim**（原始两侧序列真实重叠时的去重裁剪长度）。在 AutoGapfiller 总览中，**Net Length Change** 反映全基因组整合后的净长度变化，染色体轨迹图按原始缺口坐标绘制。
   - 在 **CtgLinker** 中，**Scaffold length** 表示连接后的完整 scaffold 长度，而逐端点的 `family-report.html` 将延伸总量与连接位点指标分开展示。
   - 在 **TelSeeker** 中，**rounds** 统计染色体侧延伸尝试次数，**extension length** 累加染色体侧局部组装的新增碱基数（不含最终接入的端粒 read 尾部），而总览报告中的 **Linked Sequence Length** 则是完整 `linker.fa` 的长度（即当前染色体侧序列加上所选 read 接入目标端的外侧尾部——`L` 端前接、`R` 端后接）。
3. **所选连接与证据支持范围**
   - **GapFiller / CtgLinker**：在 **Eligible connections** 中核查最终胜出的历史版本组合（例如 `left/v3 + right/v2`）及其一致性、重叠长度和悬垂长度，并展开查看快速/全序列连接诊断表以及已应用/未应用的末端修订记录。
   - **TelSeeker**：在目标检索摘要中查看候选池规模、检索数、比对验证数、合格数、胜出 read ID 及支持度下界。报告中的佐证数量仅覆盖局部接头（不代表整条新增尾部的共识纠错或抛光），具体锚点与接头跨度规则见 [docs/design/telomere-target-search.md](docs/design/telomere-target-search.md)（英文）。
4. **停止原因与未闭合任务解读**
   - 在 **GapFiller / CtgLinker** 中，查看各延伸家族的 **Stop reason**（如 `no_extension_reads`、`no_extension_contigs` 或共享轮数/长度上限）及计入预算的 **anchor-backoff retries**。一侧正常停止只会冻结该侧历史版本，不会阻止对侧继续延伸并与之闭合。
   - 在 **TelSeeker** 中，`no_sampled_nonrepeat_anchor`、`no_eligible_alignment` 和 `insufficient_read_support` 均为**单次检查**（如 `direct.check` 或某一延伸轮次检查）的未命中/拒绝原因，不必然代表整个末端任务终止。单次检查为负后，染色体侧延伸在预算内仍可继续推进，最终任务停止原因可能是无可用延伸数据（如 `no_extension_reads`、`no_extension_contigs`）或达到轮数/长度上限。直接检查未命中或局部支持不足，本身并不意味着该末端在生物学上不存在端粒。
5. **证据与产物文件链接**
   - 可通过各模式对应的产物链接和可展开证据表，直接查阅机器可读记录（`family-result.json`、`families.json`、`process.agp`、`connection-diagnostics.json`、`junction-evidence.json`、`ctglinker-result.json`、`telomere-result.json`、`connection-support.json`、`candidates.tsv`）。原始比对、局部组装与阶段回执按模式保存在各自子目录中——例如延伸步骤位于 `attempts/`，成对重叠检查位于 `family-overlaps/` 或 `ctg-overlaps/`，TelSeeker 检索与比对位于 `direct.check/` 或 `telomere-alignments/`——具体布局见[输出目录结构](#输出目录结构)。

报告视觉呈现与中英文切换细节可参考 [docs/design/report-presentation.md](docs/design/report-presentation.md)（英文）。

---

## 输出目录结构

### GapFiller

```text
output_directory/
├── <name>.final.fa                 # 闭合后的单条序列，或未闭合时分开导出的延伸家族序列
├── family-result.json              # 经校验的状态（"closed" 或未闭合）、胜出组合与校验和
├── family-report.html              # 离线双语报告，展示双侧家族与合格候选连接
├── family-components.fa            # 供 process.agp 引用的累计序列组件
├── families.json                   # 断点文件、延伸家族历史及共享预算状态
├── process.agp                     # 从 family-components.fa 重建 <name>.final.fa 的 AGP
├── process.log                     # 人类可读的搜索过程与结果日志
├── process.summary                 # 各家族尝试轮数、版本数、新增碱基与停止原因汇总
├── connection-diagnostics.json     # 快速/全序列连接检查及版本对拒绝原因
├── junction-evidence.json          # 仅经验证闭合时存在：原始边界定位与末端修订支持证据
├── junction-integration.fa         # 仅经验证闭合时存在：供 AutoGapfiller 使用的最小替换片段
├── family-overlaps/                # 版本对 MUMmer 比对产物与回执
├── ont-quality-cache/              # 共享 SQLite 质量索引（纯 ONT 组装时延迟创建）
├── processed_reads/                # 规范化后的工作 FASTA reads
├── hifi_reads.idx / ont_reads.idx  # BioPython reads 随机访问索引（按输入存在）
├── HiFi.reads.stat / ONT.reads.stat
├── hifi_reads_part/ / ont_reads_part/
└── attempts/<family_index>/<attempt>/
    ├── result.json                 # 经校验的单次尝试回执
    ├── mapping/                    # 可复用的 minimap2/samtools BAM 与回执
    ├── reads/                      # 筛选出的候选 reads 与阶段回执
    ├── assembler/                  # 可复用的局部组装（hifiasm）与阶段回执
    └── assembly/                   # contig 选择、累计序列与阶段回执
```

### CtgLinker

```text
output_directory/
├── DG.scaffold.fa                  # 最终 scaffold 序列（保留未连接的 contig）
├── DG.scaffold.log                 # Scaffold N50 及汇总统计
├── DG.scaffold.agp                 # 引用 DG.scaffold.components.fa 的重建坐标
├── DG.scaffold.components.fa       # AGP 实际使用的原始输入与延伸快照序列
├── DG.project.fa                   # 已连接的组装快照
├── Genome.inputCtg.fa              # 经 N 分割与 seed 长度过滤后的输入 contig
├── Genome.inputCtg.log             # 输入 contig 统计
├── used.txt / unplaced.txt         # 原始 contig 归属追踪记录
├── ctg-families.json               # 项目断点及待处理的端点搜索状态
├── ctglinker-result.json           # 经校验的结果清单与产物校验和
├── ctglinker-family-report.html    # 离线双语总报告，链接各端点搜索报告
├── ctg-sources/                    # 不可变的原始输入组件及 .rc.fa 反向互补视图
├── ctg-assemblies/                 # 组合生成的 DG*.fa scaffold 快照
├── ctg-overlaps/                   # MUMmer 重叠比对缓存
├── shared-extension-stages/        # 基于内容寻址的可复用单轮延伸阶段缓存
└── project/
    └── endpoint-000001/            # 单次共享预算的端点连接搜索
        ├── families.json
        ├── attempts/
        ├── family-result.json
        └── family-report.html
```

### TelSeeker

```text
output_directory/
├── genome.telomere.check/
│   └── need_extension_chr_end.txt          # 从 -e/--target_ends 读取的目标染色体末端
├── part1.telo.reads/
│   ├── Global.left.telo.reads.fa
│   ├── Global.right.telo.reads.fa
│   ├── left.telo.reads.fa
│   ├── right.telo.reads.fa
│   └── part1.log
├── part2.chr.end.job/
│   ├── .telomere-target-cache/             # 共享的磁盘持久化目标索引
│   └── <Chr>.<L|R>/                        # 例如 Chr01.L、Chr01.R
│       ├── seed.fa
│       ├── telo.reads.fa
│       ├── direct.check/
│       ├── linker.fa / linker.info         # 连接成功时生成
│       ├── telomere-alignments/             # 缓存的批量比对结果（第 0 轮检查同样复用）
│       └── extension/
│           ├── families.json               # 单条真实染色体侧延伸家族断点
│           ├── attempts/0/<attempt>/       # 各轮尝试（含 reads/assembly/telomere 回执）
│           ├── telomere-result.json        # 供报告读取的经校验结果清单
│           ├── extension.log
│           └── extension.summary
├── part3.integration.results/
│   ├── check_part2_jobs.csv                # 各末端连接汇总表
│   ├── final.genome.fa                     # 整合已连接 linker 后的基因组
│   └── final.genome.telomere.check/        # 供人工复核的最终 2 kb 末端 FASTA 与基序图
└── visual.report/
    ├── Global.report.html                  # 全局总览报告
    └── <Chr>.<L|R>.report.html             # 单个染色体末端详情报告
```

### TelSeeker Ctg

```text
output_directory/
├── telseeker_ctg.input/
│   ├── cache_manifest.json                 # 共享 reads/基序/端粒参数的兼容性契约
│   ├── ctgs.fa
│   └── ctg_tasks.tsv
├── part1.telo.reads/                       # 共享端粒 reads 池
├── telseeker_ctg.jobs/
│   └── <safe_ctg_id>.<L|R>/
│       ├── task-contract.json
│       ├── direct.check/
│       ├── extension/                      # 单端延伸家族断点与 telomere-result.json
│       └── result/
│           ├── status.tsv
│           ├── receipt.json
│           ├── extended.ctg.fa             # 连接成功时输出
│           └── failed.reason.txt           # 失败时输出原因
└── result/
    └── status.tsv                          # 全部请求端点的汇总状态表
```

### AutoGapfiller

```text
output_directory/
├── autogapfiller_summary.json              # 流程总汇总（使用 --setup_only 时为 autogapfiller_setup_summary.json）
├── 01.reads_preprocessor/
│   ├── processed_reads/
│   │   ├── hifi_reads.fa                   # 提供 --hifi 时存在
│   │   └── ont_reads.fa                    # 提供 --ont 时存在
│   ├── HiFi.reads.stat / ONT.reads.stat
│   ├── hifi_reads.idx / ont_reads.idx
│   ├── hifi_reads_part/ / ont_reads_part/
│   └── preprocessing_config.json
├── 02.gap_detection/
│   ├── gapbase.log                         # 全部原始 N 区域记录
│   ├── gap.log                             # 过滤（>=20 bp）并合并（间距 <11 kb）后的缺口清单
│   ├── gapDataBase/                        # 原始两侧序列
│   ├── gapData/                            # 实际用于子任务的两侧序列
│   └── gap_detection_config.json
├── 03.gap_filling/
│   ├── job_manager_config.json
│   ├── gap_filling_results.json
│   └── <chrom>.<gap_id>.bidirectional/     # 每个缺口独立的双向 GapFiller 工作目录
│       ├── run.sh
│       ├── seqleft.fasta / seqright.fasta
│       ├── families.json
│       ├── family-result.json
│       ├── family-report.html
│       ├── <chrom>.<gap_id>.bidirectional.final.fa
│       └── attempts/
├── 04.genome_integration/
│   ├── genome.filled.fasta                 # 整合填补结果后的基因组
│   ├── integration.log
│   └── integration_config.json
└── 05.visualization/
    ├── Global.report.html                  # 交互式总览报告
    └── <chrom>.<gap_id>.bidirectional.html # 逐缺口家族报告（产物链接已重定向）
```

---

## 常见问题排查

1. **纯 ONT 组装因 FASTA / 质量值错误而失败**
   - 请提供带有真实 Phred+33 质量值的原始 ONT `.fq` / `.fastq`（或 `.gz`）文件。
   - 检查多个原始 ONT FASTQ 文件间是否存在重复的 read ID，以及是否存在截断或格式损坏的记录。
   - 如果此前在 ONT 质量修复提交（`c0fd8e3`）之前的旧版本中运行过任务，请保留旧目录不动，在**新输出目录**中重新运行。

2. **区分 `hifiasm` 显式报错 / 缺失 `.p_ctg.gfa` 与正常未闭合（`exhausted`）**
   - 检查 `attempts/<family>/<attempt>/assembler/hifiasm/hifiasm.stderr.log` 与 `hifiasm.stdout.log`。
   - 出现显式 `hifiasm` 错误行（即使退出码为 0）或缺少主 contig GFA 文件时，任务会按工具错误终止，而不会发布为 `exhausted`。
   - 若 `hifiasm` 正常完成并写出合法的空 GFA（0 条序列记录）且未产生可用延伸 contig，该延伸家族会以 `no_extension_contigs` 正常停止，整次搜索可能在未闭合状态下以 `exhausted` 结束。

3. **重跑时提示旧版输出目录或输入变更错误**
   - 旧版单侧输出目录（含 `checkpoint.info`、`process/`、旧式 `round*` 目录，或缺少家族清单的旧版 CtgLinker/TelSeeker 目录）会被拒绝而非覆盖，请改用新输出目录。
   - 修改输入 reads（包含 ONT 质量值）、两侧序列、运行参数或外部工具二进制文件都会使基于内容校验的断点失效；变更条件后请使用新输出目录。

4. **K-mer 预筛选相关问题**
   - 启用 `--kmer_filter` 需要环境中 `PATH` 包含 `seqkit`、GNU `parallel` 与 `rg`（ripgrep）。
   - 如果开启 `--kmer_filter` 后某个缺口未能延伸，请在新的输出目录中关闭 `--kmer_filter` 重试。

5. **AutoGapfiller 未检测到任何缺口**
   - 确认输入基因组 FASTA 中包含长度至少为 **20 bp** 的连续 `N` 区域（`AutoGapfillerPart2.py` 会过滤短于 20 bp 的 N 区域，并合并间距小于 11,000 bp 的相邻缺口）。
   - 查看 `02.gap_detection/gapbase.log` 与 `02.gap_detection/gap.log`。

6. **AutoGapfiller 子任务失败**
   - 检查 `03.gap_filling/<chrom>.<gap_id>.bidirectional/run.stdout.log`、`run.stderr.log`、`process.log` 以及 `03.gap_filling/gap_filling_results.json`。
   - 若触发内存或 CPU 限制，请适当调小 `--work` 或 `-t`。

7. **TelSeeker 未找到端粒 reads 或末端未能连接**
   - 确认 `-m/--motif` 使用符合目标物种的大写 `A/C/G/T` 基序（如 `TTAGGG`、`TTTAGGG`、`TTAGG`）。
   - 确认 `-e/--target_ends` 中的染色体名称与 `--genome` FASTA 标题一致，且格式为 `<chr>.L` 或 `<chr>.R`。
   - 如果第 1 步生成的 `part1.telo.reads/Global.left.telo.reads.fa` / `Global.right.telo.reads.fa` 为空，请在新的输出目录中放宽筛选条件重试，例如 `--tel-r 0.5 --tel-mm 1`（或减小 `--tel-n`）。
   - 遇到 `insufficient_read_support`：说明存在候选比对，但佐证同一连接的独立目标 read ID 不足 2 个。若仅有的两条 reads 分别指向不同末端或路径，或同一 read ID 重复出现，均不构成独立的双 read 佐证；请检查 `connection-support.json`，切勿直接将其等同于“该处不存在端粒”。

8. **TelSeeker Ctg 报错 `no_tel_reads`、`no_tel_reads_for_target_end` 或缓存不匹配**
   - 查看 `result/status.tsv` 与 `telseeker_ctg.jobs/<safe_ctg_id>.<L|R>/result/status.tsv`。
   - `no_tel_reads`：第 1 步未筛选到任何全局端粒 reads。
   - `no_tel_reads_for_target_end`：存在端粒 reads，但没有对应所请求 `L` 或 `R` 端的端粒 reads。
   - 复用 `--out` 时出现缓存不匹配：`--hifi`、`--ont`、`--motif`、`--tel-n`、`--tel-r` 与 `--tel-mm` 必须与 `telseeker_ctg.input/cache_manifest.json` 完全一致，否则请指定新的 `--out`。

9. **磁盘空间占用**
   - 为支持延伸家族断点恢复，无论 `--remove 1` 还是 `--remove 2`，预处理 reads、分块、索引、历史版本及已验证阶段回执均会保留。请提前规划存储空间，并在确认不再需要续跑后手动归档或清理已完成的工作目录。

---

## 相关文档

- **概览与快速开始**：[README_zh.md](README_zh.md)
- **设计说明**：[bidirectional-extension.md](docs/design/bidirectional-extension.md)、[junction-validation.md](docs/design/junction-validation.md)（英文）、[telomere-target-search.md](docs/design/telomere-target-search.md)（英文）、[report-presentation.md](docs/design/report-presentation.md)（英文）
- **阶段性验证记录（均为英文）**：[2026-11-22-ssh180-expanded-matrix.md](docs/validation/2026-11-22-ssh180-expanded-matrix.md)、[2026-11-22-ssh180-new-dataset.md](docs/validation/2026-11-22-ssh180-new-dataset.md)、[2026-11-22-junction-repair.md](docs/validation/2026-11-22-junction-repair.md)、[2026-11-22-telomere-target-search.md](docs/validation/2026-11-22-telomere-target-search.md)、[2026-11-21-ssh180-real-ont.md](docs/validation/2026-11-21-ssh180-real-ont.md)、[2026-11-21-ont-quality-repair.md](docs/validation/2026-11-21-ont-quality-repair.md)
- **历史 v1 指南**：[v1/README.md](v1/README.md)（英文）

---

## 引用与联系方式

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. *Briefings in Bioinformatics*, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

如有问题反馈或建议，请联系 <jzhang@mail.hzau.edu.cn>。
