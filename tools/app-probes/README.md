# 其它工程的校验脚本（暂存）

这两个文件是 `FastAcquisition` 里 Savitzky–Golay 重写时的校验脚本，`using FastAcquisition.Analysis.Common`
与 `MathNet.Numerics`，**不属于本库**，也不参与本仓库的构建与测试。

原先误放在 `tests\` 下。保留在这里只是避免丢失；正式归宿是 `FastAcquisition` 自己的测试工程。

| 文件 | 用途 |
|---|---|
| `SgVerify.cs` | 新旧 SG 实现一致性、与 SciPy `savgol_filter(mode='interp')` 对照、参数校验 |
| `SgBench.cs` | SG 的耗时与内存分配测量 |
