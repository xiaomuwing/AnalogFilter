# FAMOS 对照工作台（人工执行，不属于自动化测试）

这里的探针用来回答「换成本库之后，波形会不会和 imc FAMOS 不一样」。它们**不参与**解决方案构建，
也不进 `dotnet test`：需要 32 位进程、本机已注册的 `ImcFamosLib` COM 组件，以及 30 万点的
`data\11.csv` 实测数据（在 `FastAcquisition\data\` 下）。

自动化部分覆盖的是 SciPy 正确性；FAMOS 对齐属于行为验证，结论已写进
[`../../README.md`](../../README.md) 第五节。

## 编译

无 csproj，用 `csc` 直接编译（`AnalogFilter.dll` 先 Release 构建好）：

```powershell
$csc = "C:\Windows\Microsoft.NET\Framework\v4.0.30319\csc.exe"
& $csc /nologo /platform:x86 /optimize+ /out:AllFamosProbe.exe `
      AllFamosProbe.cs /r:..\..\src\AnalogFilter\bin\Release\AnalogFilter.dll `
      /r:ImcFamosLib.dll
```

## 文件

| 文件 | 用途 |
|---|---|
| `AcfProbe.cs` | 用 FAMOS 的 `Spec.ACF` 导出参考值 |
| `AcfVerify.cs` | 本库 ACF 与 FAMOS 导出的二进制结果逐点比较 |
| `AllFamosProbe.cs` | 四类滤波器 × 四频带 × 多阶数，FAMOS 与 SciPy 三方对比 |
| `ChebFamosProbe.cs` | 切比雪夫 `characteristic=2` 的参数语义探测 |
| `CritProbe.cs` / `CritAllProbe.cs` / `CritSweep.cs` | 临界阻尼（FAMOS 第 3 类）的 −3 dB 约定与频带行为 |
| `FamOpsProbe.cs` | FAMOS 的阶数/频带约定探测 |
| `FamosMathVerify.cs` | FAMOS 输出的数值复核 |
| `ThreeWay.cs` + `threeway_compare.py` / `threeway_attribution.py` | 读取 `11.csv`，跑 FAMOS 与本库，再用 SciPy 算第三路 |
| `BenchOurs.cs` / `BenchTheirs.cs` | 同一数据上的耗时对比（本库 vs FAMOS） |
