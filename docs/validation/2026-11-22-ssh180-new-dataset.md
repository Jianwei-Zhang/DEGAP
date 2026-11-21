# SSH 180：更换公开数据后的填补验证

日期：2026-11-22。分支：`feat/bidirectional-gap-extension`。

## 结论

换用 hifiasm 官方 HiFi 示例 reads 后，两个预先指定的 **20 kb 人工缺口均已闭合**，都由 `L1` 与 `R1` 连接。填入序列相对遮蔽前基线的一致性分别为 **99.970%** 和 **99.670%**；两者都不是逐碱基完全一致。

本次检验的是**真实 reads 上的人工缺口填补和流程一致性**。基线也由同一批 reads 组装，不是独立完成图，不能据此宣称天然缺口填补率、真实错误率或相对旧版本的提升。另行检查的 ONT 小型测试数据未产生可用基线，未进入 DEGAP 填补阶段，不能算作 ONT 填补成功。

## 数据与预先固定的测试设计

此前使用的 `SRR25241091` 没有复用。本次 HiFi 输入来自 hifiasm 官方示例：

```text
https://github.com/chhylp123/hifiasm/releases/download/v0.7/chr11-2M.fa.gz
```

- 5,626 条 reads，共 84,042,542 bp；长度范围 2,416–26,785 bp。
- 原始压缩文件 SHA-256：`c1ff1393a9abd2f245119e081e5c8f20d8141a0f38b70f070f0eacf470f6e0f2`。
- 用 `hifiasm -o baseline -t 4 -f0 chr11-2M.fa.gz` 建立基线，得到一条 2,217,142 bp contig：`ptg000001l`。
- 基线 FASTA SHA-256：`b7e3e0904f7566e0c790702c754c5c45498ccc67bcb34a5be4ce0ece77268a63`。
- 在最长基线 contig 的约 1/3、2/3 位置各遮蔽 20 kb，左右各保留 60 kb。坐标、flank、预算在运行 DEGAP 前写入 `evidence/case-plan.json` 和各 case 的 `*.inputs.json`；没有根据结果移动缺口或调低连接阈值。
- 两例均使用全部 HiFi reads。共享上限为左右合计 12 次尝试、300,000 bp 延伸，4 线程；不传 `--flag`，不启用 k-mer 或 depth 过滤。

坐标为 **0-based、右端不包含**：

| Case | 基线中的 gap 区间 | Gap 大小 | 左 / 右 flank |
| --- | --- | ---: | ---: |
| `hifi-1` | `[729047, 749047)` | 20,000 bp | 各 60,000 bp |
| `hifi-2` | `[1468094, 1488094)` | 20,000 bp | 各 60,000 bp |

## 作业与结果

运行在 SSH 别名 `180` 的 ARM 计算节点，经 Donau 提交，不在登录节点运行组装。任务目录：

```text
/share/home/xbzhang26/cgat-reads-debug/degap-newdata-20261122/
```

使用 site modules：hifiasm `0.25.0-r726`、minimap2 `2.28-r1209`、SAMtools `1.21`、MUMmer `4.0.0rc1`、SeqKit `2.8.2`。每个 job 先 `module purge` 再加载指定模块，避免继承环境使 Conda 工具覆盖实际测试版本。

| 作业 | 用途 | 结果 |
| --- | --- | --- |
| `43009310` | HiFi 基线组装 | `SUCCEEDED` |
| `43009314` | `hifi-1` 填补、比对、恢复检查 | `SUCCEEDED` |
| `43009315` | `hifi-2` 填补、比对、恢复检查 | `SUCCEEDED` |
| `43009325` | 输入、源码、收据汇总核验 | `SUCCEEDED` |

| Case | 状态 / 胜出版本 | L + R 尝试 | L + R 新增碱基 | 填入长度 | Gap-only 一致性 |
| --- | --- | ---: | ---: | ---: | ---: |
| `hifi-1` | `closed` / `L1 + R1` | 1 + 1 | 13,856 + 16,254 | 19,999 bp | 99.970% |
| `hifi-2` | `closed` / `L1 + R1` | 1 + 1 | 19,471 + 13,062 | 19,982 bp | 99.670% |

- **一致性只计算 gap 内填入序列**，没有把两侧共 120 kb 的原始 flank 算进去。用 `minimap2 -cx asm20 --cs=long` 比对，指标为 PAF `matches / alignment block length`。
- 两例均为正向、MAPQ 60，覆盖全部填入序列及基线的全部 20,000 bp gap。分别为 `19996 / 20002` 和 `19947 / 20013`；两侧 flank 均逐碱基保留。
- 家族闭合 overlap 的 identity 是另一个指标：分别为 99.91% / 10,111 bp、99.52% / 12,551 bp，不能代替上面的 gap-only 检验。
- 每例实际调用 hifiasm 两次，均有非空 primary-contig GFA、无显式工具错误；不是仅依靠 exit 0 判定成功。
- 初次 DEGAP 命令耗时分别为 33.071、33.351 秒，不含基线准备、排队、结果比对和恢复检查；这不是新旧版本性能对照。

## 完整性与恢复检查

两例均通过：

- `family-result.json` 的 manifest/hash 校验；历史版本内容与长度不变。
- 左右轮次与新增长度之和等于共享计数；全部 4 个已有 L×R 版本配对均被检查，只保留一个最佳 winner。
- `process.agp` 从 `family-components.fa` 逐段重建最终 FASTA。
- `--resume_auto` 重跑不改变阶段收据、BAM、GFA、家族 checkpoint 和结果文件的内容及修改时间；没有重复 mapping/assembly。
- 原始输入 SHA-256 保持不变。HiFi 用例不创建 ONT quality index，不能视作 ONT 质量恢复测试。
- 实际 `family-report.html` 与结果中的家族计数、胜出版本、identity/overlap、历史长度及 SHA-256 对齐。

在该任务目录中复查命令与证据：

```text
scripts/environment.sh
scripts/prepare_hifi.py
scripts/prepare_cases.py
scripts/run_case.py
scripts/hifi-1.dsub
scripts/hifi-2.dsub
evidence/hifi-1.validation.json
evidence/hifi-2.validation.json
evidence/acceptance.json
runs/hifi-1/family-report.html
runs/hifi-2/family-report.html
```

## ONT 补充数据：未进入填补测试

另下载 ONT 官方 `epi2me-labs/wf-bacterial-genomes` 的 FASTQ 测试数据，固定仓库提交 `db10ac064bf4600a31821c0cab596d768874e23a`。原始路径、下载地址和 hash 保存在 `evidence/input-sources.json`、`evidence/downloads.sha256`、`evidence/acceptance.json`。

- barcode01：1,750 条 / 6,422,808 bp，最长 13,682 bp；barcode02：1,250 条 / 4,623,556 bp，最长 12,679 bp。
- FASTQ header 标注 `dna_r10.4.1_e8.2_400bps_hac@v3.5.2`。配套 `reference.subseq.fa.gz` 的 header 明示参考被随机修改，因此没有将它当作准确率真值。
- barcode01 基线作业 `43009316` 中，`hifiasm --ont -f0` 正常完成、无显式错误，但三个 primary GFA 均无 `S` 段。准备脚本因没有基线 contig 主动失败，预定 ONT 人工缺口没有生成，也没有运行 DEGAP。barcode02 只作输入统计，未组装。
- 之前的准备作业 `43009313` 是测试脚本诊断命令把 `hifiasm` 写成 `ontasm`，在组装前失败；修正该拼写后才提交 `43009316`，并未修改 DEGAP。
- 尝试取得较大的官方 ONT demo 包时，远端 300 秒、本地 240 秒下载均超时；仅留下 `.part`，没有把不完整下载当作测试输入。

这部分是数据适用性检查，不计入已执行 GapFiller 用例，也不提供 ONT 填补效果结论。此前 ONT 质量修复及未闭合结果仍以 [2026-11-21 的记录](2026-11-21-ont-quality-repair.md) 为准。

## 源码、文档与统一检查

- 基础提交：`4d18ab73f600d9899db3804b5b0578adfe191519`。本次只改 Markdown 和 HTML 展示，不改填补算法、输入处理、排名或预算。实际测试 `bin/` 快照逐文件 hash 记录在 `evidence/source-files.sha256.json`；源码包 SHA-256：`9be8cfad519e8982eecc83dabcd44b209924aaf88c762f3ca9f8a8ed7f987a6b`。
- 下载的证据归档：`/tmp/degap-newdata-20261122/evidence-review.tar.gz`，SHA-256：`f52babde806fb2ff610ef68e2715f5519f1807931c51dd2ec0ef803290ba635d`。原始 reads、完整运行目录保留在远端，不加入 Git。
- AGY `gemini-4-neon-high` 精简 Markdown，并校正 AutoGapfiller 表头、TelSeeker 无限预算显示和动态文本转义。随后人工核对参数、输出路径、本地文档链接及实际报告。
- AST 比较确认只改了展示函数；新增回归覆盖 `Filled`/`Direction` 表头、无限/零预算、动态标识符及切换语言后的转义。
- 统一自动检查：**181 项通过，0 跳过，275.821 秒**；日志 `/tmp/degap-docs-final-regression-20261122.log`。本次未启动浏览器，没有把此前浏览器截图作为本次渲染验收。

```bash
# 在 CGAT 工作区，使用已配置的本机验证环境
PATH=/tmp/degap-family-validation-tools:$PATH \
DEGAP_VALIDATION_ARTIFACTS=/tmp/degap-docs-final-fixtures-20261122 \
PYTHONWARNINGS=ignore \
micromamba run -n GRT python -m unittest discover -s DEGAP2.0/tests -p 'test_*.py'
```
