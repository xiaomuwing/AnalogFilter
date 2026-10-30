using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 双线性变换（模拟 → 数字）。
    /// 移植自 SciPy 1.18.1 的 <c>bilinear</c> / <c>bilinear_zpk</c>，BSD-3-Clause。
    /// </summary>
    /// <remarks>
    /// SciPy 在数字滤波器设计流程中固定使用 <c>fs = 2.0</c>（归一化到 Nyquist），
    /// 因此 <c>s = 4·tan(π·Wn/2)·(z−1)/(z+1)</c>；本库的 <see cref="FilterDesign"/> 亦沿用该约定。
    /// </remarks>
    public static class BilinearTransform
    {
        /// <summary>零极点增益形式的双线性变换。等价于 <c>scipy.signal.bilinear_zpk</c>。</summary>
        public static ZpkModel Zpk(ZpkModel analog, double fs)
        {
            if (analog == null)
            {
                throw new ArgumentNullException(nameof(analog));
            }
            ValidateRate(fs);

            int degree = AnalogTransforms.RelativeDegree(analog);
            double fs2 = 2.0 * fs;

            Complex[] zeros = new Complex[analog.Zeros.Length + degree];
            for (int i = 0; i < analog.Zeros.Length; i++)
            {
                zeros[i] = (fs2 + analog.Zeros[i]) / (fs2 - analog.Zeros[i]);
            }
            // 无穷远处的零点映射到 Nyquist（z = -1）。
            for (int i = 0; i < degree; i++)
            {
                zeros[analog.Zeros.Length + i] = -Complex.One;
            }

            Complex[] poles = new Complex[analog.Poles.Length];
            for (int i = 0; i < analog.Poles.Length; i++)
            {
                poles[i] = (fs2 + analog.Poles[i]) / (fs2 - analog.Poles[i]);
            }

            Complex numerator = Complex.One;
            for (int i = 0; i < analog.Zeros.Length; i++)
            {
                numerator *= fs2 - analog.Zeros[i];
            }
            Complex denominator = Complex.One;
            for (int i = 0; i < analog.Poles.Length; i++)
            {
                denominator *= fs2 - analog.Poles[i];
            }

            return new ZpkModel(zeros, poles, analog.Gain * (numerator / denominator).Real);
        }

        /// <summary>
        /// 多项式形式的双线性变换。等价于 <c>scipy.signal.bilinear</c>。
        /// </summary>
        public static TransferFunction Polynomial(double[] numerator, double[] denominator, double fs)
        {
            ValidateRate(fs);
            if (numerator == null || numerator.Length == 0)
            {
                throw new ArgumentException("分子不能为空。", nameof(numerator));
            }
            if (denominator == null || denominator.Length == 0)
            {
                throw new ArgumentException("分母不能为空。", nameof(denominator));
            }

            double[] b = TrimLeadingZeros(numerator);
            double[] a = TrimLeadingZeros(denominator);

            // 把 2*fs 的因子拆到分子分母各一份，降低大 fs/高阶时的溢出风险（与 SciPy 一致）。
            double fac = Math.Sqrt(fs * 2.0);

            // 升幂表示：zp1 = (1 + z)/fac，zm1 = (z - 1)*fac
            double[] zp1 = { 1.0 / fac, 1.0 / fac };
            double[] zm1 = { -fac, fac };

            int order = Math.Max(a.Length, b.Length) - 1;

            // 系数按升幂索引 q 参与运算（Python 里的 b[::-1] 即此顺序）。
            double[] numAsc = new double[order + 1];
            for (int q = 0; q < b.Length; q++)
            {
                int power = order - q;
                double[] term = Multiply(Power(zp1, power), Power(zm1, q));
                AddScaled(numAsc, term, b[b.Length - 1 - q]);
            }

            double[] denAsc = new double[order + 1];
            for (int p = 0; p < a.Length; p++)
            {
                int power = order - p;
                double[] term = Multiply(Power(zp1, power), Power(zm1, p));
                AddScaled(denAsc, term, a[a.Length - 1 - p]);
            }

            return new TransferFunction(Reverse(numAsc), Reverse(denAsc));
        }

        private static void ValidateRate(double fs)
        {
            if (double.IsNaN(fs) || double.IsInfinity(fs) || fs <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(fs), "采样率必须是大于 0 的有限数值。");
            }
        }

        private static double[] TrimLeadingZeros(double[] coeffs)
        {
            int start = 0;
            while (start < coeffs.Length - 1 && coeffs[start] == 0.0)
            {
                start++;
            }
            double[] result = new double[coeffs.Length - start];
            Array.Copy(coeffs, start, result, 0, result.Length);
            return result;
        }

        private static double[] Reverse(double[] values)
        {
            double[] result = new double[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                result[i] = values[values.Length - 1 - i];
            }
            return result;
        }

        private static double[] Power(double[] polynomial, int exponent)
        {
            double[] result = { 1.0 };
            for (int i = 0; i < exponent; i++)
            {
                result = Multiply(result, polynomial);
            }
            return result;
        }

        private static double[] Multiply(double[] left, double[] right)
        {
            double[] result = new double[left.Length + right.Length - 1];
            for (int i = 0; i < left.Length; i++)
            {
                if (left[i] == 0.0)
                {
                    continue;
                }
                for (int j = 0; j < right.Length; j++)
                {
                    result[i + j] += left[i] * right[j];
                }
            }
            return result;
        }

        private static void AddScaled(double[] target, double[] term, double scale)
        {
            int offset = target.Length - term.Length;
            for (int i = 0; i < term.Length; i++)
            {
                target[offset + i] += term[i] * scale;
            }
        }
    }
}
