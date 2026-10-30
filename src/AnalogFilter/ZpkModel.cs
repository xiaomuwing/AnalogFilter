using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 零极点增益（ZPK）形式的传递函数：
    /// <c>H(s) = k * Π(s - z_i) / Π(s - p_i)</c>。
    /// 对应 <c>scipy.signal</c> 里的 <c>(z, p, k)</c> 三元组。
    /// </summary>
    public sealed class ZpkModel
    {
        /// <summary>用零点、极点和增益构造 ZPK 模型。</summary>
        /// <param name="zeros">零点集合，允许为空。</param>
        /// <param name="poles">极点集合，允许为空。</param>
        /// <param name="gain">增益系数 k。</param>
        public ZpkModel(IReadOnlyList<Complex> zeros, IReadOnlyList<Complex> poles, double gain)
        {
            Zeros = ToArray(zeros);
            Poles = ToArray(poles);
            Gain = gain;
        }

        /// <summary>零点。</summary>
        public Complex[] Zeros { get; }

        /// <summary>极点。</summary>
        public Complex[] Poles { get; }

        /// <summary>增益。</summary>
        public double Gain { get; }

        /// <summary>极点阶数减零点阶数（相对阶），必须 ≥ 0。</summary>
        public int RelativeDegree => Poles.Length - Zeros.Length;

        private static Complex[] ToArray(IReadOnlyList<Complex> source)
        {
            if (source == null || source.Count == 0)
            {
                return Array.Empty<Complex>();
            }

            Complex[] result = new Complex[source.Count];
            for (int i = 0; i < source.Count; i++)
            {
                result[i] = source[i];
            }
            return result;
        }
    }

    /// <summary>
    /// 分子/分母多项式形式的传递函数：<c>H(s) = B(s) / A(s)</c>，系数<strong>降幂</strong>。
    /// 对应 <c>scipy.signal</c> 里的 <c>(b, a)</c>。
    /// </summary>
    public sealed class TransferFunction
    {
        /// <summary>
        /// 用分子、分母系数（降幂）构造传递函数。构造时按 <c>scipy.signal.normalize</c> 归一化，
        /// 并裁掉分子的前导零。
        /// </summary>
        /// <param name="numerator">分子系数（降幂），不能为空。</param>
        /// <param name="denominator">分母系数（降幂），首项不能为 0。</param>
        public TransferFunction(double[] numerator, double[] denominator)
        {
            double[] b;
            double[] a;
            TransferFunctionConversion.Normalize(numerator, denominator, out b, out a);
            Numerator = b;
            Denominator = a;
        }

        /// <summary>分子系数（降幂）。</summary>
        public double[] Numerator { get; }

        /// <summary>分母系数（降幂）。</summary>
        public double[] Denominator { get; }
    }

    /// <summary>ZPK 与 TF 之间的相互转换，以及系数归一化。</summary>
    public static class TransferFunctionConversion
    {
        /// <summary>
        /// 零极点增益 → 分子/分母多项式。等价于 <c>scipy.signal.zpk2tf</c>。
        /// 要求零极点共轭对称（物理可实现滤波器的必然结果），否则会抛出异常。
        /// </summary>
        public static void Zpk2Tf(
            IReadOnlyList<Complex> zeros,
            IReadOnlyList<Complex> poles,
            double gain,
            out double[] b,
            out double[] a)
        {
            Complex[] numerator = ComplexPolynomial.FromRoots(zeros ?? Array.Empty<Complex>());
            Complex[] denominator = ComplexPolynomial.FromRoots(poles ?? Array.Empty<Complex>());

            b = ToReal(numerator, gain, "分子");
            a = ToReal(denominator, 1.0, "分母");
        }

        /// <summary>零极点增益 → 分子/分母多项式。等价于 <c>scipy.signal.zpk2tf</c>。</summary>
        public static TransferFunction Zpk2Tf(ZpkModel zpk)
        {
            double[] b;
            double[] a;
            Zpk2Tf(zpk.Zeros, zpk.Poles, zpk.Gain, out b, out a);
            return new TransferFunction(b, a);
        }

        /// <summary>
        /// 分子/分母多项式 → 零极点增益。等价于 <c>scipy.signal.tf2zpk</c>。
        /// </summary>
        public static ZpkModel Tf2Zpk(double[] numerator, double[] denominator)
        {
            double[] b;
            double[] a;
            Normalize(numerator, denominator, out b, out a);

            double gain = b[0];
            Complex[] scaled = new Complex[b.Length];
            for (int i = 0; i < b.Length; i++)
            {
                scaled[i] = new Complex(b[i] / gain, 0.0);
            }

            Complex[] complexA = new Complex[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                complexA[i] = new Complex(a[i], 0.0);
            }

            Complex[] zeros = ComplexPolynomial.Roots(scaled);
            Complex[] poles = ComplexPolynomial.Roots(complexA);
            return new ZpkModel(zeros, poles, gain);
        }

        /// <summary>
        /// 归一化：去掉分母前导零，两者同除分母首项。
        /// 等价于 <c>scipy.signal.normalize</c>（含分子前导零的 BadCoefficients 处理）。
        /// </summary>
        public static void Normalize(double[] numerator, double[] denominator, out double[] b, out double[] a)
        {
            if (denominator == null || denominator.Length == 0)
            {
                throw new ArgumentException("分母不能为空。", nameof(denominator));
            }
            if (numerator == null || numerator.Length == 0)
            {
                throw new ArgumentException("分子不能为空。", nameof(numerator));
            }

            double[] den = TrimLeadingZeros(denominator);
            double den0 = den[0];
            if (den0 == 0.0)
            {
                throw new ArgumentException("分母的首项系数不能为 0。", nameof(denominator));
            }

            a = new double[den.Length];
            for (int i = 0; i < den.Length; i++)
            {
                a[i] = den[i] / den0;
            }

            b = new double[numerator.Length];
            for (int i = 0; i < numerator.Length; i++)
            {
                b[i] = numerator[i] / den0;
            }

            // 分子前导零（相对量级 ≤ 1e-14）判定与裁剪，与 scipy 的 BadCoefficients 分支一致。
            double maxAbs = 0.0;
            for (int i = 0; i < b.Length; i++)
            {
                maxAbs = Math.Max(maxAbs, Math.Abs(b[i]));
            }

            int leadingZeros = 0;
            for (int i = 0; i < b.Length; i++)
            {
                if (Math.Abs(b[i]) <= 1e-14 * Math.Max(1.0, maxAbs))
                {
                    leadingZeros++;
                }
                else
                {
                    break;
                }
            }

            if (leadingZeros > 0)
            {
                if (leadingZeros >= b.Length)
                {
                    leadingZeros = b.Length - 1;
                }
                LastNormalizeWarning = "分子系数条件数很差，结果可能没有意义（BadCoefficients）。";
                double[] trimmed = new double[b.Length - leadingZeros];
                Array.Copy(b, leadingZeros, trimmed, 0, trimmed.Length);
                b = trimmed;
            }
            else
            {
                LastNormalizeWarning = null;
            }
        }

        /// <summary>最近一次归一化告警（对应 scipy 的 BadCoefficients 警告）。</summary>
        public static string LastNormalizeWarning { get; private set; }

        private static double[] TrimLeadingZeros(double[] coeffs)
        {
            int start = 0;
            while (start < coeffs.Length - 1 && coeffs[start] == 0.0)
            {
                start++;
            }

            if (start == 0)
            {
                return coeffs;
            }

            double[] result = new double[coeffs.Length - start];
            Array.Copy(coeffs, start, result, 0, result.Length);
            return result;
        }

        private static double[] ToReal(Complex[] coeffs, double scale, string what)
        {
            double[] result = new double[coeffs.Length];
            for (int i = 0; i < coeffs.Length; i++)
            {
                Complex value = coeffs[i] * scale;
                double magnitude = Complex.Abs(value);
                if (Math.Abs(value.Imaginary) > 1e-9 * Math.Max(1.0, magnitude))
                {
                    throw new InvalidOperationException(
                        what + "系数出现不可忽略的虚部，零极点不是共轭对称的（无法表示成实系数滤波器）。");
                }
                result[i] = value.Real;
            }
            return result;
        }
    }
}
