using System;
using System.Collections.Generic;

namespace AnalogFilter
{
    /// <summary>滤波器类型，对应 <c>scipy.signal.iirfilter(ftype=...)</c>。</summary>
    public enum FilterType
    {
        /// <summary>巴特沃斯（最大平坦幅频）。</summary>
        Butterworth,

        /// <summary>切比雪夫 I 型（通带等波纹，SciPy 约定：波纹边在截止频率）。</summary>
        ChebyshevI,

        /// <summary>切比雪夫 II 型（阻带等波纹）。</summary>
        ChebyshevII,

        /// <summary>贝塞尔（最大平坦群延迟）。</summary>
        Bessel,

        /// <summary>
        /// 临界阻尼（N 个重实极点，<c>H(s) = 1/(1 + s/ω0)^N</c>）。
        /// SciPy 没有这个原型，这里按 FAMOS 的约定实现：ω0 取值使 -3 dB 点正好落在给定截止频率上。
        /// </summary>
        CriticallyDamped,

        /// <summary>
        /// 切比雪夫 I 型，按 <b>FAMOS 约定</b>归一化：−3 dB 点落在截止频率上，且直流增益为 1。
        /// 与 <see cref="ChebyshevI"/> 只差归一化方式。
        /// </summary>
        /// <remarks>新增成员一律追加在末尾，避免改变既有成员的底层数值。</remarks>
        ChebyshevI3dB
    }

    /// <summary>频带类型。</summary>
    public enum BandType
    {
        /// <summary>低通。</summary>
        LowPass,

        /// <summary>高通。</summary>
        HighPass,

        /// <summary>带通。</summary>
        BandPass,

        /// <summary>带阻。</summary>
        BandStop
    }

    /// <summary>
    /// 滤波器设计入口，对应 <c>scipy.signal.iirfilter</c> 及其便捷封装
    /// （butter / cheby1 / cheby2 / bessel）。算法移植自 SciPy 1.18.1，
    /// 与 SciPy 的差别仅在于用 C# 重写了数值过程。
    /// </summary>
    /// <remarks>
    /// 频率约定与 SciPy 完全一致：
    /// <list type="bullet">
    /// <item>数字滤波器：<c>wn</c> 归一化到 Nyquist（0 &lt; wn &lt; 1）；若给出 <c>sampleRate</c>，则 <c>wn</c> 以 Hz 表示。</item>
    /// <item>模拟滤波器：<c>wn</c> 是角频率（rad/s）。</item>
    /// <item>带通/带阻需传入两个频率，且 <c>wn[0] &lt; wn[1]</c>。</item>
    /// </list>
    /// </remarks>
    public static class FilterDesign
    {
        /// <summary>设计滤波器并返回零极点增益形式。</summary>
        public static ZpkModel DesignZpk(
            int order,
            double[] wn,
            FilterType type = FilterType.Butterworth,
            BandType band = BandType.LowPass,
            double passbandRipple = 1.0,
            double stopbandAttenuation = 20.0,
            bool analog = false,
            double? sampleRate = null,
            BesselNormalization besselNorm = BesselNormalization.Phase)
        {
            if (order < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(order), "滤波器阶数不能为负。");
            }
            if (wn == null || wn.Length == 0)
            {
                throw new ArgumentNullException(nameof(wn));
            }
            if (analog && sampleRate.HasValue)
            {
                throw new ArgumentException("模拟滤波器不能用 sampleRate 参数。", nameof(sampleRate));
            }
            if (passbandRipple < 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(passbandRipple), "通带纹波必须非负。");
            }
            if (stopbandAttenuation < 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(stopbandAttenuation), "阻带衰减必须非负。");
            }

            bool needsTwoFrequencies = band == BandType.BandPass || band == BandType.BandStop;
            if (needsTwoFrequencies && wn.Length != 2)
            {
                throw new ArgumentException("带通/带阻必须给出两个临界频率。", nameof(wn));
            }
            if (!needsTwoFrequencies && wn.Length != 1)
            {
                throw new ArgumentException("低通/高通只需给出一个临界频率。", nameof(wn));
            }

            double[] normalized = new double[wn.Length];
            if (sampleRate.HasValue)
            {
                if (sampleRate.Value <= 0.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(sampleRate), "采样率必须大于 0。");
                }
                for (int i = 0; i < wn.Length; i++)
                {
                    normalized[i] = wn[i] / (sampleRate.Value / 2.0);
                }
            }
            else
            {
                Array.Copy(wn, normalized, wn.Length);
            }

            for (int i = 0; i < normalized.Length; i++)
            {
                if (double.IsNaN(normalized[i]) || double.IsInfinity(normalized[i]) || normalized[i] <= 0.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(wn), "临界频率必须大于 0。");
                }
                if (!analog && normalized[i] >= 1.0)
                {
                    throw new ArgumentOutOfRangeException(
                        nameof(wn),
                        "数字滤波器的临界频率必须满足 0 < wn < 1（归一化到 Nyquist），或 0 < f < fs/2。");
                }
            }

            if (needsTwoFrequencies && !(normalized[0] < normalized[1]))
            {
                throw new ArgumentException("wn[0] 必须小于 wn[1]。", nameof(wn));
            }

            // 1) 模拟低通原型（截止 1 rad/s）
            ZpkModel prototype = DesignPrototype(order, type, passbandRipple, stopbandAttenuation, besselNorm);

            // 2) 预畸变（数字滤波器）：s = 4*tan(pi*Wn/2)
            double[] warped = new double[normalized.Length];
            for (int i = 0; i < normalized.Length; i++)
            {
                warped[i] = analog
                    ? normalized[i]
                    : 4.0 * Math.Tan(Math.PI * normalized[i] / 2.0);
            }

            // 临界阻尼：把 -3 dB 归一化做在**原型**上，而不是缩放频率轴。
            // 原型的极点在 -1，|H(jω)| = 1/(1+ω²)^(N/2)，-3 dB 点在 ω = sqrt(2^(1/N)-1)。
            // 先把它归一到 ω = 1 处 -3 dB，之后四种频带变换都自动正确：
            //   LP：极点变 -ωc/s               （与直接缩放频率等价）
            //   HP：极点变 -ωc·s               （同上）
            //   BP/BS：LP→BP/BS 变换本身就把原型的 ω = ±1 映射到两个带边，
            //          因此带边自然落在 -3 dB（之前缩放频率轴的做法在这里是错的）
            if (type == FilterType.CriticallyDamped)
            {
                double s = Math.Sqrt(Math.Pow(2.0, 1.0 / order) - 1.0);
                prototype = AnalogTransforms.LowPassToLowPass(prototype, 1.0 / s);
            }

            // 3) 频率变换
            ZpkModel transformed;
            switch (band)
            {
                case BandType.LowPass:
                    transformed = AnalogTransforms.LowPassToLowPass(prototype, warped[0]);
                    break;

                case BandType.HighPass:
                    transformed = AnalogTransforms.LowPassToHighPass(prototype, warped[0]);
                    break;

                case BandType.BandPass:
                    transformed = AnalogTransforms.LowPassToBandPass(
                        prototype,
                        Math.Sqrt(warped[0] * warped[1]),
                        warped[1] - warped[0]);
                    break;

                case BandType.BandStop:
                    transformed = AnalogTransforms.LowPassToBandStop(
                        prototype,
                        Math.Sqrt(warped[0] * warped[1]),
                        warped[1] - warped[0]);
                    break;

                default:
                    throw new ArgumentOutOfRangeException(nameof(band));
            }

            // 4) 数字域：双线性变换（与 SciPy 一致固定 fs = 2）
            return analog ? transformed : BilinearTransform.Zpk(transformed, 2.0);
        }

        /// <summary>设计滤波器并返回二阶节级联（推荐用于高阶滤波，数值更稳定）。</summary>
        public static double[][] DesignSos(
            int order,
            double[] wn,
            FilterType type = FilterType.Butterworth,
            BandType band = BandType.LowPass,
            double passbandRipple = 1.0,
            double stopbandAttenuation = 20.0,
            bool analog = false,
            double? sampleRate = null,
            BesselNormalization besselNorm = BesselNormalization.Phase)
        {
            ZpkModel zpk = DesignZpk(
                order, wn, type, band, passbandRipple, stopbandAttenuation,
                analog, sampleRate, besselNorm);
            return SecondOrderSections.Zpk2Sos(zpk, analog);
        }

        /// <summary>设计滤波器并返回分子/分母多项式。</summary>
        public static TransferFunction DesignTf(
            int order,
            double[] wn,
            FilterType type = FilterType.Butterworth,
            BandType band = BandType.LowPass,
            double passbandRipple = 1.0,
            double stopbandAttenuation = 20.0,
            bool analog = false,
            double? sampleRate = null,
            BesselNormalization besselNorm = BesselNormalization.Phase)
        {
            return TransferFunctionConversion.Zpk2Tf(DesignZpk(
                order, wn, type, band, passbandRipple, stopbandAttenuation,
                analog, sampleRate, besselNorm));
        }

        /// <summary>巴特沃斯滤波器（等价于 <c>scipy.signal.butter</c>）。</summary>
        public static ZpkModel Butterworth(
            int order, double[] wn, BandType band = BandType.LowPass,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.Butterworth, band, 1.0, 20.0, analog, sampleRate);
        }

        /// <summary>切比雪夫 I 型滤波器（等价于 <c>scipy.signal.cheby1</c>）。</summary>
        public static ZpkModel ChebyshevI(
            int order, double passbandRipple, double[] wn, BandType band = BandType.LowPass,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.ChebyshevI, band, passbandRipple, 20.0, analog, sampleRate);
        }

        /// <summary>
        /// 切比雪夫 I 型滤波器，按 FAMOS 约定归一化（−3 dB 落在截止频率、直流增益为 1）。
        /// </summary>
        public static ZpkModel ChebyshevI3dB(
            int order, double passbandRipple, double[] wn, BandType band = BandType.LowPass,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.ChebyshevI3dB, band, passbandRipple, 20.0, analog, sampleRate);
        }

        /// <summary>切比雪夫 II 型滤波器（等价于 <c>scipy.signal.cheby2</c>）。</summary>
        public static ZpkModel ChebyshevII(
            int order, double stopbandAttenuation, double[] wn, BandType band = BandType.LowPass,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.ChebyshevII, band, 1.0, stopbandAttenuation, analog, sampleRate);
        }

        /// <summary>贝塞尔滤波器（等价于 <c>scipy.signal.bessel</c>）。</summary>
        public static ZpkModel Bessel(
            int order, double[] wn, BandType band = BandType.LowPass,
            BesselNormalization norm = BesselNormalization.Phase,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.Bessel, band, 1.0, 20.0, analog, sampleRate, norm);
        }

        /// <summary>
        /// 临界阻尼滤波器（对应 FAMOS characteristic = 3 的语义）。
        /// SciPy 无此原型，属本库扩展功能：等价于 N 个重实极点，且 -3 dB 点落在 <c>wn</c> 上。
        /// </summary>
        /// <remarks>
        /// 四种频带均与 FAMOS 逐位一致。关键在于 −3 dB 归一化做在<b>原型</b>上而非频率轴上：
        /// LP→BP/BS 变换本身就把原型的 ω = ±1 映射到两个带边，因此带边自然落在 −3 dB。
        /// </remarks>
        public static ZpkModel CriticallyDamped(
            int order, double[] wn, BandType band = BandType.LowPass,
            bool analog = false, double? sampleRate = null)
        {
            return DesignZpk(order, wn, FilterType.CriticallyDamped, band, 1.0, 20.0, analog, sampleRate);
        }

        private static ZpkModel DesignPrototype(
            int order,
            FilterType type,
            double passbandRipple,
            double stopbandAttenuation,
            BesselNormalization besselNorm)
        {
            switch (type)
            {
                case FilterType.Butterworth:
                    return AnalogPrototypes.Butterworth(order);
                case FilterType.ChebyshevI:
                    return AnalogPrototypes.ChebyshevINormalized(order, passbandRipple);
                case FilterType.ChebyshevI3dB:
                    return AnalogPrototypes.ChebyshevI3dB(order, passbandRipple);
                case FilterType.ChebyshevII:
                    return AnalogPrototypes.ChebyshevIINormalized(order, stopbandAttenuation);
                case FilterType.Bessel:
                    return AnalogPrototypes.Bessel(order, besselNorm);
                case FilterType.CriticallyDamped:
                    return AnalogPrototypes.CriticallyDamped(order);
                default:
                    throw new ArgumentOutOfRangeException(nameof(type));
            }
        }
    }
}
