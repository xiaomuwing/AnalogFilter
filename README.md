# AnalogFilter

.NET Framework 4.8 类库，实现**模拟原型 IIR 滤波器**的设计、变换、离散化与执行。
算法逐行移植自 **SciPy 1.18.1** 的 `scipy.signal`（`_filter_design` / `_signaltools`，BSD-3-Clause）。

**零外部依赖**：不引用任何 NuGet 包、不 P/Invoke 原生库，可直接被 32 位 WinForms 程序引用。

---

## 一、覆盖的 SciPy API

| SciPy | 本库 | 说明 |
|---|---|---|
| `buttap(n)` | `AnalogPrototypes.Butterworth(n)` | 模拟低通原型，1 rad/s |
| `cheb1ap(n, rp)` | `AnalogPrototypes.ChebyshevINormalized(n, rp)` | 切比雪夫 I 型 |
| `cheb2ap(n, rs)` | `AnalogPrototypes.ChebyshevIINormalized(n, rs)` | 切比雪夫 II 型 |
| `besselap(n, norm)` | `AnalogPrototypes.Bessel(n, norm)` | 支持 `phase` / `delay` / `mag` |
| `lp2lp_zpk` / `lp2hp_zpk` | `AnalogTransforms.LowPassToLowPass` / `LowPassToHighPass` | 模拟域频率变换 |
| `lp2bp_zpk` / `lp2bs_zpk` | `AnalogTransforms.LowPassToBandPass` / `LowPassToBandStop` | 同上 |
| `bilinear_zpk` | `BilinearTransform.Zpk(zpk, fs)` | 双线性变换 |
| `bilinear` | `BilinearTransform.Polynomial(b, a, fs)` | 多项式形式 |
| `zpk2tf` / `tf2zpk` | `TransferFunctionConversion.Zpk2Tf` / `Tf2Zpk` | 形式互转 |
| `zpk2sos` | `SecondOrderSections.Zpk2Sos` | 数字用 nearest，模拟用 minimal |
| `normalize` | `TransferFunctionConversion.Normalize` | 含 BadCoefficients 处理 |
| `lfilter` | `SosFiltering.LFilter`（含 `zi`/`zf`） | DF-II 转置 |
| `sosfilt` | `SosFiltering.SosFilt`（含 `zi`） | SOS 级联 |
| `lfilter_zi` | `SosFiltering.LFilterZi` | 累积和闭式 |
| `sosfilt_zi` | `SosFiltering.SosFiltZi` | |
| `filtfilt` | `SosFiltering.FiltFilt` | 零相位（method='pad'） |
| `sosfiltfilt` | `SosFiltering.SosFiltFilt` | 零相位，支持 odd/even/constant 延拓与自定义 padlen |
| `freqs` / `freqs_zpk` | `FrequencyResponse.Freqs` | 模拟频响 |
| `freqz` / `freqz_zpk` / `freqz_sos` | `FrequencyResponse.Freqz` | 数字频响 |
| `iirfilter` | `FilterDesign.DesignZpk` / `DesignSos` / `DesignTf` | 统一入口 |
| `butter` / `cheby1` / `cheby2` / `bessel` | `FilterDesign.Butterworth` / `ChebyshevI` / `ChebyshevII` / `Bessel` | 便捷封装 |
| —（SciPy `correlate` 的循环自相关） | `Spectrum.Autocorrelation` | **按 imc FAMOS 的 `FFL_ACF` 算法**实现：归一化循环自相关 |
| —（SciPy 无此原型） | `FilterDesign.CriticallyDamped` | **扩展功能**：N 个重实极点，−3 dB 落在截止频率；四种频带均与 FAMOS 一致 |
| —（SciPy 无此归一化） | `FilterDesign.ChebyshevI3dB` | **FAMOS 约定**的切比雪夫 I 型：−3 dB 在截止频率、直流增益为 1 |

## 二、尚未覆盖

| 缺失 | 说明 |
|---|---|
| `ellipap` / `ellip` | 椭圆滤波器（FAMOS 无此类型） |
| `buttord` / `cheb1ord` / `cheb2ord` / `ellipord` | 指标驱动自动定阶 |
| `lp2lp` / `lp2hp` / `lp2bp` / `lp2bs`（b,a 形式） | 仅有 ZPK 形式（设计流程用不到） |
| `sos2zpk` / `tf2sos` / `sos2tf` / `zpk2ss` / `ss2zpk` | 形式互转补完 |
| `findfreqs` / `impulse` / `step` / `bode` | 响应分析辅助 |
| `iirnotch` / `iirpeak` / `iircomb` | 数字域特殊滤波器 |
| `filtfilt(method='gust')` | 仅实现 `pad` |

---

## 三、精度验证（基准：SciPy 1.18.1）

> **验证基准是 SciPy，不是 FAMOS。**
> 本库的算法来自 SciPy，因此"是否正确"以 SciPy 的输出为准；FAMOS 只作为"行为参考"，
> 用于回答"能否直接替换现有软件"，不作为正确性判据。

全部用例的参考值由本机 SciPy 生成，逐值比对。

| 验证集 | 覆盖 | 结果 |
|---|---|---|
| 数字域全矩阵 | 4 类型 × 4 频带 × 阶数 3/4，比对 SOS + TF 系数 + 频响 + 滤波输出 | **32/32 通过，≤ 7.4e-14** |
| 模拟域（`analog=true`） | 巴特沃斯 2/3/4/6 阶、切比雪夫 I/II 2/4 阶、贝塞尔 2/4 阶，比对极点/零点/增益 | **10/10 通过，≤ 5.9e-16** |
| 零相位滤波 | 5 个滤波器 × {odd, even, constant} × 默认/自定义 padlen + `filtfilt` | **25/25 通过，≤ 1.8e-13** |
| 巴特沃斯高阶 | 10 阶 / 16 阶 | ≤ 4.0e-16 |

### 贝塞尔原型的精度随阶数变化

贝塞尔极点依赖 `kve`（第二类修正贝塞尔函数）迭代求根。SciPy 调用 Cephes 的复数实现，
本库用半整数阶闭式求值 —— 该闭式在复参数、高阶时存在相消，精度随阶数下降：

| 阶数 | 极点最大相对偏差 | 可用性 |
|---|---|---|
| 2 | 1.1e-16 | ✅ |
| 4 | 5.9e-16 | ✅ |
| 6 | 1.8e-15 | ✅ |
| 8 | 1.1e-13 | ✅ |
| 10 | 9.6e-13 | ✅ |
| 12 | 6.3e-12 | ✅ |
| 16 | 1.3e-09 | ⚠️ 精度下降 |
| 18 | 1.0e-08 | ⚠️ 精度下降 |
| 20 | 1.0e-07 | ⚠️ 精度下降（界面允许的最高阶） |

Aberth 迭代上限取 40 次、Newton 抛光 20 次：闭式 `kve` 在复参数下有相消，|Δ| 在 1e-14
处停滞，再迭代只是空转。相比 200/60 的配置，6 阶设计耗时从 637 µs 降到 140 µs，精度不变。

若需要与 SciPy 同级的 1e-17 精度（尤其 16 阶以上），唯一可靠途径是引入 Cephes/XSF 的原生
复数 `kve`（BSD-3，可用 MSVC 编译），这会破坏当前"零依赖"特性，属可选的增强项。

---

## 四、用法示例

```csharp
using AnalogFilter;

// 1) 数字低通：巴特沃斯 4 阶，截止 = 0.2 × Nyquist
double[][] sos = FilterDesign.DesignSos(4, new[] { 0.2 });

// 2) 零相位滤波（避免相位滞后，适合峰值/峰值时刻分析）
double[] filtered = SosFiltering.SosFiltFilt(sos, data);

// 3) 带通 + 自定义边界延拓长度
double[][] bp = FilterDesign.DesignSos(4, new[] { 0.1, 0.3 },
    FilterType.ChebyshevI, BandType.BandPass, passbandRipple: 1.0);
double[] y = SosFiltering.SosFiltFilt(bp, data, padtype: "odd", padlen: 50);

// 4) 直接给采样率（wn 用 Hz）
double[][] lp = FilterDesign.DesignSos(4, new[] { 50.0 },
    FilterType.Butterworth, BandType.LowPass, sampleRate: 1000.0);

// 5) 模拟原型（角频率 rad/s）
ZpkModel analog = FilterDesign.DesignZpk(3, new[] { 100.0 },
    FilterType.Bessel, BandType.LowPass, analog: true);

// 6) 频响
double[] w = { 0.05, 0.1, 0.2 };
double[] magDb = FrequencyResponse.MagnitudeDb(FrequencyResponse.Freqz(sos, w));
```

---

## 五、与 FAMOS 的对应关系（行为参考，非正确性判据）

用 `data/11.csv`（Strain1，30 万点，1000 Hz）对 FAMOS COM、SciPy、本库做三方逐值比对。
本节只用于回答"替换 FAMOS 后波形会不会变"，**不改变第三节以 SciPy 为准的验证结论**。

| 项目 | 结论 | 误差 |
|---|---|---|
| 巴特沃斯 低通/高通 | FAMOS 与 SciPy 完全一致 | **0.000%** |
| 巴特沃斯 带通/带阻 | FAMOS 的"阶数"指**总极点数**，等价于 SciPy 用 `order/2` | **0.000%** |
| 贝塞尔 全部频带 | FAMOS 用 **-3 dB 幅值归一化**（SciPy `norm='mag'`），非默认 `phase` | **0.000%** |
| 贝塞尔 带通/带阻 | 同样需 `order/2` | **0.000%** |
| 切比雪夫 | 需用 `ChebyshevI3dB`（FAMOS 约定）：−3 dB 落在截止频率、直流增益为 1 | **0.0002%** |
| 临界阻尼 低通/高通 | 本库扩展功能（SciPy 无此原型），与 FAMOS 一致 | **0.000%** |
| 临界阻尼 带通/带阻 | 同上；−3 dB 归一化做在原型上，带边自然落在 −3 dB | **0.000%** |

因此与 FAMOS 对齐的调用方式：

```csharp
// 巴特沃斯低通：order 直接对应
var sos = FilterDesign.DesignSos(famosOrder, new[] { f / (fs / 2.0) },
                                 FilterType.Butterworth, BandType.LowPass);

// 带通/带阻：阶数取一半（FAMOS 的 order 是总极点数）
var bp = FilterDesign.DesignSos(famosOrder / 2, new[] { f1 / (fs / 2.0), f2 / (fs / 2.0) },
                                FilterType.Butterworth, BandType.BandPass);

// 贝塞尔：必须用 Magnitude 归一化
var bl = FilterDesign.DesignSos(famosOrder, new[] { f / (fs / 2.0) },
                                FilterType.Bessel, BandType.LowPass,
                                besselNorm: BesselNormalization.Magnitude);

// 临界阻尼：N 个重实极点，-3 dB 点落在截止频率上
var cd = FilterDesign.DesignSos(famosOrder, new[] { f / (fs / 2.0) },
                                FilterType.CriticallyDamped, BandType.LowPass);
```

切比雪夫（characteristic = 2）的 `parameter` 语义尚未确定：实测 FAMOS 通带增益峰值约为 1.12（>1），
而 SciPy `cheby1` 的直流增益恰为 1，说明两者归一化方式不同；此外最佳拟合的 `rp` 随频带变化
（低通≈1.5 dB、高通≈2.0 dB、带通≈3.0 dB），不是单一常数的关系。

FAMOS 的"临界阻尼"（characteristic = 3）在 SciPy 与本库中均无对应原型。

## 六、自相关函数（ACF）

`Spectrum.Autocorrelation` 按 **imc FAMOS 3.2**（`Im7Famos.DLL`）的 `FFL_ACF` 算法实现。
算法由静态导出定位 + 受控输入探针逐项确认：

```
n = 2^floor(log2(len(x)))        // 取最长的 2 的幂前缀，输出长度 = n
X = FFT(x[0 .. n-1])
r = IFFT(|X|²)                   // 循环（非偏置）自相关
ACF[k] = r[k] / r[0]             // 归一化到零延迟，[0] 恒为 1
```

要点：

- **不去均值、不加窗**——常量输入的 ACF 恒为 1；带直流偏置时 ACF 被直流主导；
- **循环**而非线性/非偏置——输入被隐式当作周期信号；
- 输出长度是输入长度以下最大的 2 的幂（1000 点输入 → 512 点输出）。

验证（262144 点真实应变数据）：

| 项目 | 结果 |
|---|---|
| 与 FAMOS 的最大绝对差 | **4.16e-7**（lag=47318） |
| RMS 差 | 1.23e-7 |
| 冲激输入 | 1, 0, 0, … |
| 常量输入 | 恒为 1（偏差 0） |
| 单频正弦 ACF[3] vs cos(2π·cycles·3/n) | 差 ≤ 2.2e-16 |
| 耗时（262144 点，x64） | 34 ms |

与 FAMOS 的差异处于其 float32 输出精度量级（本库内部用 double）。

## 七、构建与测试

```powershell
dotnet build AnalogFilter.sln -c Release
dotnet test  AnalogFilter.sln -c Release
```

类库为 SDK 风格工程（net48，零 NuGet 依赖），产物固定在 `src\AnalogFilter\bin\Release\`：
`AnalogFilter.dll` 与同名 `.xml` 文档。测试工程 `tests\AnalogFilter.Tests` 只在测试时需要，
不会被复制到类库输出目录。

### 自动化测试（119 项）

| 测试类 | 覆盖 |
|---|---|
| `SciPyDigitalDesignTests` | 数字设计 39 个用例：SOS 系数、`freqz`、`sosfilt`、`sosfilt_zi` / `lfilter_zi`、三种延拓的零相位滤波 |
| `SciPyAnalogDesignTests` | 模拟设计 32 个用例：极点/零点集合、增益、`freqs` |
| `BesselPrototypeAccuracyTests` | 贝塞尔原型 2–12 阶的极点精度门槛（把上面「贝塞尔原型的精度随阶数变化」表固化下来） |
| `AutocorrelationTests` | ACF 与 numpy 循环自相关对照，以及冲激/常量/单频正弦等定义性质 |
| `DesignInvariantTests` | 不依赖参考值的定义性断言：−3 dB 点、直流增益、稳定性、节数、零相位特性、`freqz` 与逐节多项式一致等 |
| `ArgumentValidationTests` | 非法参数必须抛出明确异常 |

参考值基准为 **SciPy 1.18.1**，由 `tools\gen_reference.py` 生成到
`tests\AnalogFilter.Tests\Data\scipy_reference.txt`（89 个用例，浮点数按 17 位有效数字无损写出）。
需要重新生成时（升级 SciPy 或增加用例）：

```powershell
python -m venv %TEMP%\af_ref_venv
%TEMP%\af_ref_venv\Scripts\pip install numpy scipy
%TEMP%\af_ref_venv\Scripts\python tools\gen_reference.py
```

### 人工探针（不属于自动化测试）

需要 32 位 FAMOS COM 的三方对比与 FAMOS 对齐工作台已移到 `tools\famos\`，
FastAcquisition 的 Savitzky–Golay 校验脚本暂存在 `tools\app-probes\`，
两者都不参与解决方案构建，用法见各自目录的 README。

---

## 八、来源与许可

算法移植自 **SciPy**（BSD-3-Clause，Copyright (c) 2001-2002 Enthought, Inc. and The SciPy Developers），
涉及 `scipy/signal/_filter_design.py` 与 `scipy/signal/_signaltools.py` 的对应函数。

本库为独立 C# 实现，不包含 SciPy 源码；仅在算法与数值约定上与之对齐。
巴特沃斯/切比雪夫/贝塞尔的极点构造、`besselap` 的零点迭代、频率变换的增益公式、
`zpk2sos` 的配对规则、`sosfiltfilt` 的延拓与初始条件，均按 SciPy 1.18.1 的语义实现。
