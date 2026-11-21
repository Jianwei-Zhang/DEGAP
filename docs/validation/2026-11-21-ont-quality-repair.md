# ONT 输入质量与零退出码错误修复

日期：2026-11-21。分支：`feat/bidirectional-gap-extension`。
修复提交：`c0fd8e3cdc60dd62a9be0986477bd071e5e992d1`。

**此前的 FASTA/FASTQ 输入错误和零退出码错误漏检已修复，179 项回归通过，SSH 180 同输入、同工具版本复测通过质量恢复与工具执行检查。两个真实缺口仍未闭合；不将缺陷修复等同于填补能力改善。**

修复前记录见 `docs/validation/2026-11-21-ssh180-real-ont.md`，原运行目录和失败证据保持不变。

## 修改范围

- 新增 `bin/ONTAssemblyReads.py`：从一个或多个原始 FASTQ/FASTQ.gz 取得真实 Phred+33 质量值；按筛选 FASTA 的顺序恢复记录，校验 ID 和序列。反向互补时同步反转质量；裁剪仅允许唯一可定位的片段，歧义、缺失 ID、重复 ID、序列不符和损坏输入均报错。
- 原始 reads 按输入内容身份建立压缩 SQLite 索引，左右家族/后续轮次共享；加锁、原子发布、哈希验证，不为每轮反复扫描原始压缩 FASTQ。索引懒构建，未进入 ONT 单独组装的任务无需构建。
- `FamilyGap.prepare_context()` 把原始 ONT 输入（包括质量值）纳入内容契约。仅质量变化也会拒绝复用旧检查点；质量缓存损坏不能静默重建成一次“成功恢复”。
- `FindExtensionContigs` 的纯 ONT 及 mixed→ONT-only 路径使用恢复后的 FASTQ 和 `--ont`。mixed→HiFi-only 不再错误添加 ONT；仍有 HiFi 的 mixed `--ul` 路径保留既有 FASTA 兼容行为。
- hifiasm 的 stdout/stderr 出现明确错误时，即使返回 0 也抛出工具错误；缺少 primary-contig GFA 也不能伪造空结果。错误不落成成功 assembler/assembly 收据，不发布普通 `exhausted` 结果；已完成 reads 阶段可保留用于恢复。
- 删除两个旧转换入口里的假质量 fallback，兼容入口也执行相同的真实质量与序列检查。HiFi FASTA 支持不变。
- 未修改延伸阈值、闭合判据、左右调度顺序、预算规则或生物学 fallback 策略；未发布版本、未替换 CGAT 捆绑运行时。

### 输入与升级约束

ONT **单独组装**必须能取得原始 FASTQ 质量值；只有 FASTA 时明确报错，不以统一 Q 值伪装实测质量。HiFi + ONT 的 `--ul` 补充输入仍允许 ONT FASTA，但某轮若转入 ONT-only，同样需要原始 ONT FASTQ。未启动组装的第零代检查不要求建立质量索引。

修复后的输入/代码契约与旧检查点不同。保留旧失败目录，使用新输出目录复测；不要搬入先前被误判成功的组装收据。

## 自动化验证

统一检查 **179 项通过，0 跳过，293.157 秒**：

```bash
PATH=/tmp/degap-family-validation-tools:$PATH \
DEGAP_VALIDATION_ARTIFACTS=/tmp/degap-ont-fix-final-fixtures \
PYTHONWARNINGS=ignore \
micromamba run -n GRT python -m unittest discover -s DEGAP2.0/tests -p 'test_*.py'
```

日志：`/tmp/degap-ont-fix-full.log`。这是现有本地 GRT 工具环境，没有安装或更新共享依赖。

新增/强化的检查包括：

- 正反向质量、唯一裁剪、多文件 gzip、原始/筛选 ID 重复、缺失/不匹配序列、畸形 FASTQ、FASTA 拒绝、缓存复用/损坏、质量单独变化。
- hifiasm stdout/stderr 报错但退出 0、缺少 GFA、失败阶段不能生成完成收据。
- ONT 合成 CLI 用例必须有真实 hifiasm GFA 序列段，而非仅靠 read-common fallback 闭合；选用明确生成的合成 FASTQ 质量，仅用于该合成 fixture，不是生产质量补造。
- GapFiller ONT 完成后恢复、质量改变后拒绝恢复；AutoGapfiller ONT 全流程/恢复；CtgLinker ONT 完整 scaffold 重建；既有 HiFi/mixed、TelSeeker、报告及家族回归保留。

## SSH 180 真实数据复测

正式同版本复测目录：

```text
~/cgat-reads-debug/degap-ont-repair-20261121-v2/
```

- 同一 `focused.ont.fastq`：2,325 条 reads / 44,928,074 bp。SHA-256 为 `b8d5410f38c481fc25478782e4e812dd65b982a4afd64a75ddf137c38be2a1a7`，通过软链接只读复用，没有再次扫描/比对 11 GB 原始 FASTQ。
- 同样的 gap、60 kb 初始两侧和左右合计 12 次尝试 / 300 kb 新增上限。
- ARM `arm` 队列，每任务 4 CPU；实际 hifiasm 约 16 GiB 峰值，4 CPU 对应的约 20 GiB 内存限额有实际必要。
- 保持集群 hifiasm `0.25.0-r726`、minimap2 `2.28-r1209`、samtools `1.21`、MUMmer `4.0.0rc1`、seqkit `2.8.2`。任务脚本清除继承的 module 状态后加载指定模块，并记录实际可执行路径、断言关键版本。
- 测试源包 SHA-256：`883054b08a264217f0f0cd1d40a58d071387b841fead7e98033ad408804221ed`。传输校验通过；源清单中全部 `bin/`、`tests/` 文件逐个与修复提交核对一致，`VERIFIED_FIX_COMMIT` 记录提交。

| 样例 | 作业 | 总尝试 | 新增 | hifiasm 显式错误 | 生物学结果 |
|---|---:|---:|---:|---:|---|
| 原有 N-gap | 43009160 | 3 | L 4,847 bp；R 0 | 0 | 未闭合，`exhausted` |
| 20 kb 人工缺口 + 真实 reads | 43009161 | 2 | 两侧均为 0 | 0 | 未闭合，`exhausted` |

五次组装均使用 `assembler/ont.extension.fastq`，逐条核对原始序列与质量（包含反向互补），hifiasm 正常执行到完成并生成 GFA；不再出现 FASTA 输入错误或缺失组装输出。全部 GFA 的 primary contig 序列段数为 0，不能声称组装出了可用 contig。原有 gap 三次筛选分别为 1、4、14 条 reads；人工缺口两次均为 6 条。当前剩余问题是这些轮次没有形成可用局部组装，不是本次输入格式错误；是否需要调整筛选/局部组装策略，须另行验证，未在本次修复中偷偷放宽阈值。

两个任务的清单/历史哈希、所有已生成跨家族版本比较、预算合计、AGP 重建、离线报告、完成后恢复检查通过。恢复时阶段产物及质量索引内容和 mtime 不变。仍未覆盖真实数据成功闭合、深层历史择优、填补率改善或真实数据下其他模式的成效。

### 保留的先导复测

`degap-ont-repair-20261121-v1` 的 43009157/43009158 也消除了输入错误并通过质量检查，但调度器继承的 module 状态使重复加载未恢复 PATH 优先级，实际命中了 Conda minimap2 2.31。没有将该轮当作同工具版本对照。保留其证据，v2 在任务自身清理 module 状态后复核；没有改写用户 shell 配置或共享环境。

## 最终证据与结论边界

远端 `evidence/fix-acceptance.json` 分别记录：

- `tool_input_and_error_handling_fix: passed`
- `biological_gap_closure: not_achieved`

`evidence/{natural,masked}.validation.json`、实际工具路径/版本、调度器终态、源码清单、逐轮 hifiasm stderr/GFA、原始质量核对与恢复证据均保留。工具正常完成但没有 contig，和此前工具明确报错却被误当成普通停止，是两个不同结果。
