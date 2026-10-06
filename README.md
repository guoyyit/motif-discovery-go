# motif-discovery-go

本地核酸序列模体发现与统计显著性评估服务（Go，无第三方依赖）。
给定一组 FASTA 序列，发现显著富集的短序列模式，以位置权重矩阵（PWM）
表示，并给出基于序列重排零分布的显著性结论。

## 功能总览

- FASTA 解析：大小写混用、注释行（`;`/`#`）、空行；非法字符、空名称头、
  有头无序列、数据出现在头之前、超长单行与超大记录都会返回带
  `行:列` 位置（从 1 开始）的 `ParseError`；出错时不返回部分序列集合。
- PWM：由等长比对片段以“计数 + 背景加权伪计数”构建对数几率矩阵；
  支持原始得分、归一化得分（`[0,1]`）、可配置阈值与一致性序列。
- 模体发现：多个模体逐个发现（EM + 占用区间互斥），每个模体进行
  多轮初始化（一轮确定性数据驱动种子 + 多轮随机种子），取似然最高者。
- 显著性：对每条序列做保组成（N 位置固定）的 Fisher–Yates 重排，
  用与观测完全相同的拟合流程与拟合强度拟合零分布，输出经验 p 值并做
  Benjamini–Hochberg FDR 校正。
- 低复杂度防护（默认开启）：均聚物、周期 2–4 的简单串联重复、低香农熵
  区域会被掩码为 `N`，扫描时跳过，因此 AT 串联重复等简单重复在默认
  配置下不会被误报为显著模体，无需调用方开启任何开关。
- 有限步拟合：EM 在可配置的迭代上限内必然结束；达到上限未收敛时以
  可区分的 `ErrNotConverged` 返回，不挂起、不返回半成品结果。
- 并发：多任务独立派生随机流，互不共享随机状态；单个任务失败或取消
  不影响其他任务。
- 中断/恢复：工作单元边界原子落盘检查点；取消后可用同进程或新进程
  `Resume`，恢复结果与不中断运行逐位一致，且无部分结果、无重复计算。
- 可复现：所有随机过程由 `Config.Seed` 显式控制；相同种子与输入
  （含相同配置）输出逐位一致。

## 快速开始

```go
seqs, err := motif.ParseFASTA(f)
if err != nil { /* 解析错误带行/列位置 */ }

// 零值配置即可使用：所有字段取默认值，低复杂度防护默认开启。
result, err := motif.Discover(seqs, motif.Config{Seed: 42})
for _, m := range result.Motifs {
    fmt.Println(m.Matrix.Consensus(), m.PValue, m.AdjustedPValue, m.Significant)
}
```

显式配置示例：

```go
result, err := motif.Discover(seqs, motif.Config{
    Seed:              42,
    MotifLength:       10,
    ExpectMotifs:      3,
    InitRounds:        10,
    ConvergenceTol:    1e-4,
    MaxIterations:     500,
    ScoreThreshold:    0.7,
    SignificanceAlpha: 0.05,
    Shuffles:          200,
    Pseudocount:       0.5,
})
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

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `Seed` | `0` | 全部随机过程的主种子；相同种子与输入输出逐位一致 |
| `MotifLength` | `10` | 模体长度 L，最小 4 |
| `ExpectMotifs` | `3` | 期望模体数量 K（实际数量可能更少） |
| `InitRounds` | `10` | 每个模体 EM 的随机初始化轮数（另含 1 轮确定性种子）；观测与零分布共用 |
| `ConvergenceTol` | `1e-4` | EM 对数似然收敛的绝对容差 |
| `MaxIterations` | `500` | 单次 EM 的最大迭代步数；超限未收敛按 `ErrNotConverged` 处理 |
| `ScoreThreshold` | `0.7` | 位点归一化得分阈值，取值 `(0,1]` |
| `SignificanceAlpha` | `0.05` | BH FDR 校正后的显著性阈值，取值 `(0,1]` |
| `Shuffles` | `200` | 零分布重排次数 B；最小可达 p 值为 `1/(B+1)` |
| `Pseudocount` | `0.5` | PWM 伪计数强度，按背景频率加权 |
| `Background` | 自动估计 | A/C/G/T 频率；零值时从输入序列估计 |
| `MaskLowComplexity` | 默认掩码 | 历史兼容字段：显式 `true` 与默认行为一致 |
| `NoMaskLowComplexity` | `false` | 显式置 `true` 才关闭低复杂度掩码（恢复旧行为） |

### 配置校验规则

- **零值 = 未配置**：自动填入上表默认值，`Config{}` 即可运行。
- **显式越界值 = 错误**：负数计数参数（`MotifLength`/`ExpectMotifs`/
  `InitRounds`/`MaxIterations`/`Shuffles`）、超出 `(0,1]` 的
  `ScoreThreshold`/`SignificanceAlpha`、负的 `ConvergenceTol`/
  `Pseudocount`、任何 NaN，以及 `MaskLowComplexity` 与
  `NoMaskLowComplexity` 同时设置，都会立即返回错误并指明字段名与
  合法取值范围——绝不静默替换为默认值后继续运行。
- 已有的显式合法配置按原语义生效，校验不会改动任何显式字段。

## 统计语义

- **零假设**：序列无模体结构，仅具有与输入相同的碱基组成。
  重排仅打乱 A/C/G/T（N 保留原位，因此掩码结构不变）。
- **观测统计量**：对每个发现的模体，取其在每条序列中的最佳归一化
  PWM 得分之和。
- **零分布**：每个重排副本都用与观测**完全相同**的拟合流程拟合——
  相同的初始化轮数（`InitRounds`）、相同的收敛容差
  （`ConvergenceTol`）、相同的迭代上限（`MaxIterations`）与相同的
  多模体占用顺序。零分布不使用更弱的拟合强度，因此显著性结论不会
  随拟合强度设置而漂移；相同输入与相同种子下重复运行逐位一致。
- **经验 p 值**（单侧，Davison–Hinkley +1 修正）：
  `p = (1 + #{null >= observed}) / (1 + B)`，
  因而有限 B 下 `p >= 1/(B+1)`，不会报告虚假的零 p 值。
- **多重检验**：同一任务的全部候选模体统一做 Benjamini–Hochberg
  FDR 校正，输出 `AdjustedPValue`；当且仅当
  `AdjustedPValue <= SignificanceAlpha` 时标记 `Significant`。
  未达显著的模体同样返回完整的 `PValue`/`AdjustedPValue`/
  `Information` 等统计量，结论可核对。
- **位点报告**：最终输出的位点独立于 EM 的每序列一位点分配，
  而是扫描全部窗口，报告所有归一化得分达到 `ScoreThreshold` 的出现。

## 收敛与迭代上限

- 每次 EM 拟合最多执行 `MaxIterations` 次 E 步（默认 500）；
  对数似然相邻两次变化不超过 `ConvergenceTol` 即判定收敛。
- 退化输入（如全部窗口相同）会在常数步内收敛；极端参数（如远超
  浮点精度的容差）下拟合在迭代上限处停止，不会挂起。
- 达到上限未收敛的初始化轮次被丢弃；所有轮次均不收敛时，拟合以
  `ErrNotConverged`（可用 `errors.Is` 区分）失败，任务进入
  `failed` 状态，不返回任何半成品模型或部分结果。

## 低复杂度判定（默认开启）

固定 8bp 滑动窗口，命中任一条件即掩码：

- 含长度 ≥ 6 的均聚物（如 `AAAAAA`）；
- 含周期 2–4、至少完整重复 3 次的串联重复（如 `ATATAT`、`ATGATGATG`）；
- 忽略 N 后的香农熵低于 1 bit。

掩码区域在扫描中被视为未知碱基 `N` 跳过；掩码会沿重复周期向两侧
延伸，覆盖超出窗口边界的部分重复单元。默认路径与显式
`MaskLowComplexity: true` 的结论逐位一致；确需旧行为时显式设置
`NoMaskLowComplexity: true`。

## FASTA 容量边界

- 单行最长 1 MiB（`maxFASTALineLength`），超出即报错，位置为越界列。
- 单条记录最多 64 Mi 碱基（`maxFASTARecordBases`），超出即报错。
- 任何解析错误都带从 1 开始的行号与列号，且不返回部分序列集合。

## 可复现性设计

- 随机数基于 `math/rand/v2` 的 PCG128；随机流编号采用两级层级
  （65536 个高级流 × 2^48 子流），用途包括：模体 k 的初始化、
  重排 b 的序列置换、重排 b 中模体 k 的 EM 初始化。
- 每个工作单元只依赖 `(Seed, 用途, 索引)` 派生出的独立随机流，
  与调度顺序、是否中断无关；并发任务之间不共享随机状态。
- 跳过已完成单元不会改变后续单元看到的随机序列，这是“恢复后
  逐位一致”的根本保证。

## 检查点与恢复

- 状态目录下每个任务一个 `<id>.json`，采用“写临时文件 + 原子 rename”。
- 工作单元 = 发现一个模体，或完成一次重排；单元完成后立即落盘。
- `Cancel` 只在单元边界生效：不会中断半个 EM 或留下半个重排结果。
- `Resume` 直接复用已完成的模体矩阵与累积零分布，从未完成单元开始重演。
- 增大重排次数、初始化轮数或并发任务数时，内存随 `Shuffles × 模体数`
  线性增长（仅保存零分布统计量），运行时间随工作单元数线性增长；
  单个任务的失败或取消不影响其他任务。

## 本地验证

```
sh verify.sh        # 构建 + vet + 全部测试；任一步失败即非零退出
```

验收流程只有在构建、静态检查与测试全部通过时才打印
`VERIFY OK` 并以 0 退出；任一用例失败都会以非零退出码结束，
且不会输出任何表示全部通过的回显。

也可以直接运行：

```
go test ./...        # 功能与统计测试
go test -race ./...  # 并发数据竞争检查
```

测试覆盖：合成植入模体的位点召回与一致性序列、显著性校准
（随机背景下的假阳性率）、强信号检出率、默认配置下简单重复不被
误报、观测与零分布拟合强度一致、非法配置被拒绝而非替换、退化输入
下拟合有限步终止、畸形 FASTA 的位置化报错、固定种子逐位可复现、
多任务并发隔离与高负载下的取消—恢复逐位一致（含跨 Service 实例
恢复），以及验收流程在失败用例下不输出成功结论。全部测试离线运行，
不依赖外部服务与网络。
