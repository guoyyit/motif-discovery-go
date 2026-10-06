# motif-discovery-go

本地核酸序列模体发现与统计显著性评估服务（Go，无第三方依赖）。
给定一组 FASTA 序列，发现显著富集的短序列模式，以位置权重矩阵（PWM）
表示，并给出基于序列重排零分布的显著性结论。

## 功能总览

- FASTA 解析：大小写混用、注释行（`;`/`#`）、空行；非法字符返回带
  `行:列` 位置的 `ParseError`；序列统一归一化为大写。空名称头、
  只有头没有序列、数据出现在任何头之前、非法碱基、超长单行、超大记录
  都会被精确定位并整体拒绝（不返回部分序列集合）。
- PWM：由等长比对片段以“计数 + 背景加权伪计数”构建对数几率矩阵；
  支持原始得分、归一化得分（`[0,1]`）、可配置阈值与一致性序列。
- 模体发现：多个模体逐个发现（EM + 占用区间互斥），每个模体进行
  多轮初始化（一轮确定性数据驱动种子 + 多轮随机种子），取似然最高者。
- 显著性：对每条序列做保组成（N 位置固定）的 Fisher–Yates 重排，
  用与观测相同的流程拟合零分布，输出经验 p 值并做 Benjamini–Hochberg FDR 校正。
- 低复杂度（默认开启）：均聚物、周期 2–4 的简单串联重复、低香农熵区域
  会被掩码为 `N`，扫描时跳过，因此默认配置下 AT 串联重复等简单重复
  绝不会被误报为显著模体，无需调用方开启任何开关。
- 并发：多任务独立派生随机流，互不共享随机状态；单个任务失败或取消
  不影响其他任务。
- 中断/恢复：工作单元边界原子落盘检查点；取消后可用同进程或新进程
  `Resume`，恢复结果与不中断运行逐位一致，且无部分结果、无重复计算。
- 有界拟合：EM 有可配置的收敛容差与迭代上限，退化/近似重复输入也在
  有限步内结束；达到上限以可区分的 `ErrEMNotConverged` 失败，不返回半成品。
- 严格校验：缺省字段填默认值，显式越界值直接报错（注明字段与范围），
  绝不静默替换为默认值。
- 可复现：所有随机过程由 `Config.Seed` 显式控制；相同种子与输入
  （含相同配置）输出逐位一致。

## 快速开始

```go
seqs, err := motif.ParseFASTA(f)
if err != nil { /* 解析错误带行/列位置 */ }

result, err := motif.Discover(seqs, motif.Config{
    Seed:              42,
    MotifLength:       10,
    ExpectMotifs:      3,
    InitRounds:        10,
    ConvergenceTol:    1e-4,
    MaxEMIterations:   1000,
    ScoreThreshold:    motif.Float64(0.7),
    SignificanceAlpha: motif.Float64(0.05),
    Shuffles:          200,
    Pseudocount:       0.5,
    // MaskLowComplexity 缺省即启用保护；需要关闭时显式 motif.Bool(false)。
})
for _, m := range result.Motifs {
    fmt.Println(m.Matrix.Consensus(), m.PValue, m.AdjustedPValue, m.Significant)
    // 可核对统计量：m.Statistic（观测统计量）与 m.NullCount（零分布样本数）
}
```

需要中断/恢复或并发多个任务时使用 `Service`：

```go
svc, _ := motif.NewService("./state")
_ = svc.Start(ctx, "job-1", seqs, cfg)
_ = svc.Cancel("job-1")          // 在单元边界安全中断
_ = svc.Resume(ctx, "job-1")     // 同进程或重启后的新 Service 均可
snap, _ := svc.Wait(ctx, "job-1")
```

## 配置项

| 字段 | 默认值 | 含义与合法范围 |
| --- | --- | --- |
| `Seed` | `0` | 全部随机过程的主种子；相同种子与输入输出逐位一致 |
| `MotifLength` | `10` | 模体长度 L，最小 4 |
| `ExpectMotifs` | `3` | 期望模体数量 K（实际数量可能更少）；`>= 1` |
| `InitRounds` | `10` | 每个模体 EM 的随机初始化轮数（另含 1 轮确定性种子）；`>= 1` |
| `ConvergenceTol` | `1e-4` | EM 对数似然收敛的绝对容差；有限且 `>= 0` |
| `MaxEMIterations` | `1000` | 单次 EM 拟合的迭代上限；`>= 1`。达到上限以 `ErrEMNotConverged` 失败 |
| `ScoreThreshold` | `0.7` | 位点归一化得分阈值，指针类型 `*float64`，取值 `(0,1]`，用 `motif.Float64(...)` 设置 |
| `SignificanceAlpha` | `0.05` | BH FDR 校正后的显著性阈值，指针 `*float64`，取值 `(0,1)` |
| `Shuffles` | `200` | 零分布重排次数 B；`>= 1`，最小可达 p 值为 `1/(B+1)` |
| `Pseudocount` | `0.5` | PWM 伪计数强度，按背景频率加权；有限且 `>= 0` |
| `Background` | 自动估计 | A/C/G/T 频率；零值时从输入序列估计，显式值须为和为 1 的分布 |
| `MaskLowComplexity` | `true` | 三态指针 `*bool`：`nil` 缺省=启用，`motif.Bool(true)` 显式启用（与缺省逐位一致），`motif.Bool(false)` 显式关闭 |

### 缺省值与非法值的区分

整型/浮点标量字段以零值表示“未配置”，会被填入上表默认值；
`ScoreThreshold`、`SignificanceAlpha` 与 `MaskLowComplexity` 使用指针，
以区分“未设置”（`nil` → 填默认值）与“显式给出的值”。

显式给出的越界值**直接报错**，错误信息注明字段名与合法取值范围，
绝不静默替换为默认值。例如 `ScoreThreshold=0`、`SignificanceAlpha=1.5`、
`Shuffles=-1`、`MotifLength=3` 都会被拒绝。

## 统计语义

- **默认防误报**：低复杂度掩码在默认路径即启用。AT 串联重复等简单重复
  在进入 EM 前被掩码为 `N`，因此默认配置下不会把重复片段报告为显著模体，
  也不需要调用方额外开启开关。缺省路径与显式 `MaskLowComplexity=true`
  的结论逐位一致；仅当显式 `motif.Bool(false)` 时才关闭该保护。
- **零假设**：序列无模体结构，仅具有与输入相同的碱基组成。
  重排仅打乱 A/C/G/T（N 保留原位，因此掩码结构不变）。
- **观测统计量**：对每个发现的模体，取其在每条序列中的最佳归一化
  PWM 得分之和（结果中的 `Motif.Statistic`）。
- **零分布（拟合强度一致）**：每个重排副本都用与观测**完全相同**的
  EM 流程拟合——相同的初始化轮数（`InitRounds`，外加 1 轮确定性种子）、
  相同的收敛容差（`ConvergenceTol`）、相同的迭代上限
  （`MaxEMIterations`），并按相同的多模体占用顺序逐个拟合。
  零分布不会使用更弱的拟合，因此显著性结论不会随拟合强度漂移。
- **经验 p 值**（单侧，Davison–Hinkley +1 修正）：
  `p = (1 + #{null >= observed}) / (1 + B)`，
  因而有限 B 下 `p >= 1/(B+1)`，不会报告虚假的零 p 值。
  `Motif.NullCount` 给出该模体实际使用的零样本数（不超过 B；个别重排
  可用窗口不足提前停止时偏少，此时结论偏保守），结合 `Statistic` 与
  `PValue` 即可反解超出次数进行核对。
- **多重检验**：同一任务的全部候选模体统一做 Benjamini–Hochberg
  FDR 校正，输出 `AdjustedPValue`；当
  `AdjustedPValue <= SignificanceAlpha` 时标记 `Significant`。
- **位点报告**：最终输出的位点独立于 EM 的每序列一位点分配，
  而是扫描全部窗口，报告所有归一化得分达到 `ScoreThreshold` 的出现。

## 迭代上限与收敛语义

- 单次 EM 拟合最多执行 `MaxEMIterations` 个 E/M 轮次（默认 1000）。
- 当相邻两次对数似然的绝对变化量 `<= ConvergenceTol` 时判定收敛。
- 达到上限仍未收敛，或似然退化为 NaN/Inf 时，该轮初始化被丢弃；
  若一次模体拟合的全部初始化轮次（含确定性种子轮）都未收敛，
  `Discover`/`Service` 以可 `errors.Is(err, motif.ErrEMNotConverged)`
  区分的错误失败：任务状态为 `failed`，不挂起、不返回半成品结果。
- 未显式配置 `MaxEMIterations`/`ConvergenceTol` 时沿用默认值，
  常规合法输入行为与此前一致。

## 可复现性设计

- 随机数基于 `math/rand/v2` 的 PCG128；随机流编号采用两级层级
  （65536 个高级流 × 2^48 子流），用途包括：模体 k 的初始化、
  重排 b 的序列置换、重排 b 中模体 k 的 EM 初始化。
- 每个工作单元只依赖 `(Seed, 用途, 索引)` 派生出的独立随机流，
  与调度顺序、是否中断无关。
- 跳过已完成单元不会改变后续单元看到的随机序列，这是“恢复后
  逐位一致”的根本保证。

## 检查点与恢复

- 状态目录下每个任务一个 `<id>.json`，采用“写临时文件 + 原子 rename”。
- 工作单元 = 发现一个模体，或完成一次重排；单元完成后立即落盘。
- `Cancel` 只在单元边界生效：不会中断半个 EM 或留下半个重排结果。
- `Resume` 直接复用已完成的模体矩阵与累积零分布，从未完成单元开始重演。

## 低复杂度判定

固定 8bp 滑动窗口，命中任一条件即掩码：

- 含长度 ≥ 6 的均聚物（如 `AAAAAA`）；
- 含周期 2–4、至少完整重复 3 次的串联重复（如 `ATATAT`、`ATGATGATG`）；
- 忽略 N 后的香农熵低于 1 bit。

掩码区域在扫描中被视为未知碱基 `N` 跳过；掩码会沿重复周期向两侧
延伸，覆盖超出窗口边界的部分重复单元。

## 容量边界

- FASTA 单行上限 `MaxFASTALineBytes = 1 MiB`，单条记录序列数据上限
  `MaxFASTARecordBytes = 64 MiB`；越界字节所在的 `行:列` 即报错位置，
  解析器在定长缓冲内流式读取，不会因超长行无限占用内存。
- 计算量随输入规模线性/可控增长：
  - 发现阶段约为 `K × InitRounds × MaxEMIterations × Σ(序列窗口数) × L`，
    但 EM 通常远在 `MaxEMIterations` 前收敛；
  - 显著性阶段为 `B × K` 次与观测同强度的拟合，外加 `B × Σ序列长度`
    的 Fisher–Yates 重排。增大 `Shuffles` 或 `InitRounds` 只线性增加
    单元数量，每个单元独立落盘，内存占用不随 B 累积（仅保留每个模体
    的一维零分布分数，`O(K·B)` 个 float64）。
- 中断/恢复以工作单元为粒度，因此容量增大后内存与运行时间均可由
  上述参数解释；单个任务失败或取消不影响并发的其他任务。

## 本地验证

无需外部服务、无网络：

```
scripts/verify.sh    # = go build ./... && go vet ./... && go test -count=1 ./...
go test ./...        # 功能与统计测试
go test -race ./...  # 并发数据竞争检查
```

`scripts/verify.sh` 的退出码语义：

- 任一构建/vet/测试用例失败时立即以**非零退出码**结束，只输出失败信息，
  **绝不**打印任何“全部通过”回显；
- 仅当全部检查真正通过时才输出唯一成功标记 `ALL_CHECKS_PASSED` 并以 0 退出。

覆盖范围：默认配置下简单重复不被误报（且与显式开启防护逐位一致）、
观测与零分布拟合强度一致（初始化轮数/收敛容差/迭代上限/占用顺序）、
非法配置被拒绝而非被默认值替换、退化与近似重复输入下 EM 有限步终止、
达到迭代上限返回可区分的 `ErrEMNotConverged`、畸形 FASTA 的位置化报错、
合成植入模体的位点召回与一致性序列、显著性校准（随机背景假阳性率）、
强信号检出率、固定种子逐位可复现、多任务并发隔离与随机流互不干扰，
以及在发现/重排阶段多个中断点（含增大 Shuffles/InitRounds 的高负载）
的取消—恢复逐位一致（含跨 Service 实例恢复）、单任务失败不影响他任务，
并验证验收脚本本身在失败时非零退出且不打印成功回显。
