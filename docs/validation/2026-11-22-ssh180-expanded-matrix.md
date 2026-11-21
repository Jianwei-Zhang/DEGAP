# SSH 180：扩展测试矩阵与阶段结果

日期：2026-11-22。分支：`feat/bidirectional-gap-extension`。测试源码：`c85007b036455f9a02ffeaec271723f21ccdcdc6`。

## 当前结论

在[两个 20 kb 用例](2026-11-22-ssh180-new-dataset.md)之后，新增 **29 个独立命令用例**，覆盖 GapFiller、AutoGapfiller、CtgLinker、k-mer 开关、预算边界和单侧停止后继续。29 个用例的流程及恢复检查均通过；这不等于 29 个缺口全部闭合。

- HiFi：9 个固定位置/长度的单 gap 全部闭合；3 个 k-mer 开启对照也闭合且序列与对应关闭版本相同。
- AutoGapfiller：HiFi 的 3 个 gap 全部填补；低覆盖 ONT 的 3 个 gap 全部未闭合。
- CtgLinker：3 条 contig（中间一条反向输入）连接成一条 scaffold，原始成员无丢失、无重复。
- 共享预算与冻结侧：3 个边界用例按预期停止；2 个缺少单侧支持的用例保留停止侧，另一侧继续尝试，没有立即整体失败。
- ONT：官方 Salmonella demo 的 9 个固定 gap 均未闭合。相对比较参考的 aligned-span 深度约 4.29×，本轮所有 local-assembly primary GFA 均没有序列段；不能把流程正常退出当作生物学成功。

本轮 HiFi 仍是同一批 reads 组装基线上的人工缺口，不是天然 gap 或独立真值准确率评估。ONT demo 和 LT2 比较参考的菌株一致性未证实，二者差异也不能直接计为测序错误。成功、正常未闭合、工具失败和验证脚本失败分开记录。

## 设计与可复查位置

运行目录（独立于旧测试及其他项目任务）：

```text
/share/home/xbzhang26/cgat-reads-debug/degap-expanded-20261122/
```

- 通过 Donau ARM 队列运行计算；普通用例各 4 CPU、最多 3,600 秒。AutoGapfiller `--work 1`，防止 worker 叠加超出资源。
- 沿用 site hifiasm `0.25.0-r726`、minimap2 `2.28-r1209`、SAMtools `1.21`、MUMmer `4.0.0rc1`、SeqKit `2.8.2`；每个 job 清理继承的 module 状态后重新加载。Python 环境只读复用，不更新共享软件。
- GapFiller 固定 60 kb flank。除边界用例外，预算为 **全部家族合计 12 次尝试、300,000 bp 新增序列**。不传 `--flag`，不调 depth 过滤参数。
- 主矩阵在基线的 1/4、1/2、3/4 位置，交叉 5/20/50 kb 长度。gap 起点为 `floor(reference_length * numerator / denominator) - gap_length / 2`，终点为起点加长度；坐标采用 0-based、右端不含。
- ONT 另固定 1/4、1/2、3/4 × 1/5/20 kb。坐标在查看填补结果前确定，不根据结果移动 gap、修改门槛或删掉失败例。
- 读取一次已有基线和原始 reads；已完成的 mapping、assembly 和结果不因增加评估而重跑。

远端 `evidence/plan.json`、`evidence/ont-selection-plan.json`、`evidence/freeze-selection-plan.json` 保存各轮预先设计，`evidence/<case>.inputs.json` 保存实际坐标和输入 hash。`scripts/prepare.py`、`scripts/prepare_ont_data.py`、`scripts/prepare_freeze.py` 和 `scripts/run_case.py` 保存准备、执行与恢复核验命令。

## HiFi：固定 gap 矩阵

原始示例 `chr11-2M.fa.gz` 共 5,626 条 / 84,042,542 bp；基线 `ptg000001l` 长 2,217,142 bp。输入及基线 hash 见[前一份记录](2026-11-22-ssh180-new-dataset.md)。以下 9 例均为 `closed`。

| Case | Donau job | 总尝试 | 选中版本 | Gap 内实际长度 | Gap-only 一致性 |
| --- | --- | ---: | --- | ---: | ---: |
| `hifi-q1-5k` | `43009632` | 1 | L1 + R0 | 5,000 bp | 100.0000% |
| `hifi-q1-20k` | `43009633` | 4 | L2 + R2 | 19,990 bp | 99.8150% |
| `hifi-q1-50k` | `43009634` | 5 | L3 + R2 | 50,028 bp | 99.7702% |
| `hifi-mid-5k` | `43009635` | 1 | L1 + R0 | 5,000 bp | 99.9800% |
| `hifi-mid-20k` | `43009636` | 4 | L0 + R2 | 19,984 bp | 99.9050% |
| `hifi-mid-50k` | `43009637` | 4 | L2 + R2 | 49,995 bp | 99.9780% |
| `hifi-q3-5k` | `43009638` | 1 | L1 + R0 | 5,000 bp | 99.8400% |
| `hifi-q3-20k` | `43009639` | 2 | L1 + R1 | 19,958 bp | 99.7400% |
| `hifi-q3-50k` | `43009640` | 5 | L3 + R1 | 49,980 bp | 99.9020% |

所有 gap-only 比对都覆盖基线 gap 的全部区间，不能以较长、较准确的 flank 稀释 gap 内差异。闭合 overlap 的 identity 只是搜索判定指标，不等于这里的填入序列一致性。

### 历史版本确实被使用

- `hifi-mid-20k`：运行结束已有 L2、R2，但唯一合格连接是 **L0 + R2**。9 个已有 L×R 组合均有检查记录。
- `hifi-q3-50k`：已有 L3、R2，采用 **L3 + R1**，不是最新的 R2。12 个已有 L×R 组合均被检查；这证明本用例需要保留历史版本，不证明每次新延伸都更差。
- 所有单 gap 用例检查候选排名且只发布一个 winner；预算没有因采用旧版本而返还。

### k-mer 开关对照

在 mid 的 5/20/50 kb 三个位置长度组合上增加 `--kmer_filter`，作业 `43009641`–`43009643`。三例的最终序列 SHA-256 与关闭版本分别一致，轮次和新增长度也一致。作业同时在不同计算节点运行，没有重复性能采样，因此不据耗时宣称 k-mer 加速比。

## 共享限制与冻结侧

| Case / job | 设计 | 实际结果 |
| --- | --- | --- |
| `hifi-limit-round0` / `43009644` | 轮次上限 0 | `round_limit`，0 次延伸；仍执行初始版本连接检查 |
| `hifi-limit-round1` / `43009645` | 左右合计 1 次 | `round_limit`，1 次，新增 14,448 bp；没有给右侧另开 1 次额度 |
| `hifi-limit-length1` / `43009646` | 长度上限 1 bp | 完整执行 1 次，新增 14,448 bp 后 `length_limit`；符合“下一次尝试前检查长度”的语义 |
| `hifi-frozen-left` / `43009791` | 保留只支持右侧的原始 reads | L：1 次 / 0 bp，R：4 次 / 24,942 bp，最终 `exhausted` |
| `hifi-frozen-right` / `43009792` | 保留只支持左侧的原始 reads | R：1 次 / 0 bp，L：3 次 / 24,723 bp，最终 `exhausted` |

冻结侧测试使用 `hifi-mid-50k` 的固定坐标，通过原始 reads 的 primary alignment 筛选左右单侧支持，MAPQ≥20、query aligned fraction≥0.9；不修改 read 序列。它们是缺少证据的受控故障场景，不加入普通填补率分母。两个用例都保存原侧历史并让另一侧继续，最终未闭合也如实保留。

## AutoGapfiller 与 CtgLinker

### HiFi AutoGapfiller（`43009647`）

在完整 2.2 Mb 基线上一次遮蔽 3 个 gap：q1 5 kb、mid 20 kb、q3 50 kb。检测到 3 个任务，各只有一个双向 workspace；全部闭合并整合，最终不含 N。

| Gap | 实际填入长度 | Gap-only 一致性 |
| --- | ---: | ---: |
| gap1 | 5,000 bp | 100.0000% |
| gap2 | 19,984 bp | 99.9050% |
| gap3 | 49,980 bp | 99.9020% |

按每个 gap 的 `replace_start/replace_end` 及最终 family FASTA 重新拼接整个基因组，与输出逐碱基一致；恢复时阶段收据、BAM、GFA、家族 checkpoint 均未改变。

### HiFi CtgLinker（`43009648`）

固定 3 条各 100 kb 的原始 contig，中间间隔各 20 kb，输入方向为 `+ / - / +`。最终一条 scaffold，成员方向恢复为 `+ / - / +`（相对于输入），长度 339,991 bp；AGP 从实际 components 重建成功，成员无遗漏或重复。

端点搜索共享上限，不是整个 contig 项目只允许 12 次：第一次左端搜索用满 12 次未闭合，后续两次右端搜索分别 3 次、2 次完成连接，项目累计 17 次。两处 gap 的单独比对一致性为 **99.5852%**、**99.9750%**，覆盖全部基线 gap。

## 侧翼变化与评估方法

扩大矩阵后发现：**闭合不保证两侧原始 flank 的每个碱基都保持不变**。当前连接使用选中累计版本的左侧 overlap 拷贝，再接右侧 continuation；当累计版本重新覆盖原始 flank 时，可能保留 assembly-consensus 与原始 flank 的差异。这不是本轮新改算法，但需要在结果解释中体现。

例如：

- `hifi-mid-5k` 的 gap 内有 1 个替换；gap 外另有 1 个替换和 5 bp 缺失。
- `hifi-q1-20k` 的 gap 外另有 10 个替换、27 bp 插入、4 bp 缺失。
- AutoGapfiller 的原始非 gap 区域合计另有 6 个替换、19 bp 插入、19 bp 缺失。
- CtgLinker 原始 contig 区域合计另有 10 个替换、1 bp 插入、5 bp 缺失。

因此不能对所有闭合结果直接截取 `final[len(left):-len(right)]` 并称其为 gap 真值检验。本轮新增只读评估 `scripts/assess.py`：

1. 使用已完成的全序列 `minimap2 -cx asm20 --cs=long` 比对，要求唯一 primary、MAPQ≥20 且覆盖 gap 两侧各至少 1 kb；不唯一或未覆盖则标为无法评估，不默认为正确。
2. 用 CIGAR 将基线 gap 边界投影到结果序列，再单独比对该 gap 序列，报告 `matches / block_length`、覆盖度、长度差异和 exact match。边界插入归属采用起点包含、终点不包含。
3. 分别从 `cs` 统计 gap 内和原始非 gap 区域的替换、插入、缺失。完全保留 flank 的用例同时要求坐标投影与直接截取一致。

本地 9 个轻量 CIGAR/cs 检查覆盖正常匹配、插入、缺失、负链坐标及边界处理。比对评估不是独立组装质量认证，低复杂度区仍可能存在等价 indel 放置。

## ONT Salmonella：低覆盖负结果

官方 demo 包下载完成后核对归档结构，并按发布的 MD5 校验成员。只使用带 `runid`、flow-cell 等测序记录的 barcode02；barcode01 header 含 `error-free_length` / `read_identity` 模拟标记，没有拿它作为真实测序证据。

- FASTQ：7,238 条 / 24,192,101 bp，最长 10,715 bp，median 3,304 bp；header 模型 `dna_r10.4.1_e8.2_400bps_hac@v4.2.0`。
- 比较参考：`NC_003197.2`，4,857,450 bp。MAPQ≥20 primary 的 aligned-span 深度 **4.2948×**，覆盖参考 **89.7283%**。这是比对跨度统计，不是无 gap 的精确碱基深度。
- demo 包 SHA-256：`8d2950c4fa63d05b092f8a87a01c78fd81fc40f123a1ea296ef2061543538528`；所用 FASTQ SHA-256：`f8f074b2f3209ded05af66ba6c8aa57e02353820ef23113e5307a85bfc90e70d`。

| 用例组 | Donau jobs | 生物学结果 | 流程检查 |
| --- | --- | --- | --- |
| q1 × 1/5/20 kb | `43009777`–`43009779` | 3 例 `exhausted`，未闭合 | 通过 |
| mid × 1/5/20 kb | `43009780`–`43009782` | 3 例 `exhausted`，未闭合 | 通过 |
| q3 × 1/5/20 kb | `43009783`–`43009785` | 3 例 `exhausted`，未闭合 | 通过 |
| AutoGapfiller 3 gap | `43009786` | 0/3 填补，保留原缺口 | 通过 |

工具无显式错误，但所有 local primary GFA 都无 `S` 段。真实 FASTQ 质量按原始序列/反向互补/唯一片段核对，恢复时质量缓存没有重建。该负结果说明此输入未提供足够的当前算法可用延伸证据，不能单凭它区分覆盖不足、短 reads、参考差异或算法限制。

ONT 准备作业 `43009766` 因验证脚本错误地假设参考文件只有一条记录而失败，在开始填补前终止。修正为显式选择 `NC_003197.2` 后，准备作业 `43009771` 完成；没有改 DEGAP 或放宽连接阈值。失败日志保留。

## 已通过的检查与证据

29 个用例均验证输入/源码 hash、结果 manifest、家族历史、共享计数、AGP 重建和完成后的恢复不重复计算。CtgLinker 另验证原始成员；AutoGapfiller 另验证检测数量、状态传递和全基因组替换。正常未闭合和预算停止按其实际状态接受，不要求每个 case 都 `closed`。

完整运行输出、reads、GFA/BAM 不加入 Git。逐用例结果在远端 `evidence/<case>.validation.json`，坐标投影结果在 `evidence/<case>.assessment.json`。本地只读证据：

| 文件（位于 `/tmp/degap-expanded-20261122/`） | SHA-256 |
| --- | --- |
| `hifi-evidence-complete.tar.gz` | `5f4bdd6c9d53114f420d1029dfd210a47288407a503e14bff7a2eb4081057ed0` |
| `hifi-assessments.tar.gz` | `a79d02db45e5bd16d7dececdcd61c053e3fc7de74befaf5bfdf47855a221a084` |
| `supplementary-evidence.tar.gz` | `5d09789155f0c18d7924c6e1c6135fa6eae47a3dfe1804136d8e3fbb69ac9452` |

测试源码从前轮已验证目录只读复用；逐文件 hash 位于 `evidence/source-files.sha256.json`，对应提交 `c85007b`。没有新修改算法、工具版本、HTML 模板或客户数据。本轮不重复此前 181 项全回归，先做必要的结果及评估器轻量检查；最终检查随新增 ONT 阶段统一汇总。

## 下一阶段：先记录计划，再推进

- [x] 先提交本文件，固定以上已完成设计、负结果及评估边界。
- [ ] 下载并通过 ENA 发布 MD5/字节数核验 `ERR12342517`：Klebsiella pneumoniae MGH 78578，R10.4.1 / Dorado rerio，68,001 reads / 417,332,669 bp，压缩文件 398,270,970 bytes。下载尚在进行，不把 `.part` 当作输入。
- [ ] 使用同名菌株独立参考 `CP000647.1`；先核对实际 read 长度、参考覆盖和 gap 邻域深度，记录 provenance。与参考存在的样本差异不自动称为 assembly error。
- [ ] 在参考固定 1/4、1/2、3/4 位置交叉 5/20/50 kb，预先锁定 60 kb flank、共享 12 次 / 300 kb、4 threads；不移动 gap 或为成功调低阈值。
- [ ] 增加该数据的 AutoGapfiller 三 gap 场景，核对 gap-only 序列、原始 flank 变化、质量恢复、共享历史和重跑复用；成功与未闭合均报告。
- [ ] 所有新增用例终止并核验后，更新本文件及最终证据清单，再提交当前分支。

若下载未完成或输入完整性/来源无法确认，停在该数据准备阶段，不用模拟 reads 替代实测、不把已有 HiFi 成功包装成 ONT 成功。
