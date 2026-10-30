using System;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 频谱与相关分析。目前提供自相关函数（ACF）。
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="Autocorrelation(double[])"/> 按 imc FAMOS 的算法实现
    /// （<c>Im7Famos.DLL</c> v3.2 的 <c>FFL_ACF</c>），经黑盒探针逐项确认：
    /// </para>
    /// <list type="bullet">
    /// <item>取输入最长的 2 的幂前缀：n = 2^floor(log2(len))，输出长度等于 n；</item>
    /// <item>用 <b>循环</b>（非线性和非偏置）自相关：r[k] = IFFT(|FFT(x)|²)[k]；</item>
    /// <item><b>不去均值、不加窗</b>：常量输入的 ACF 恒为 1，带直流偏置的信号 ACF 被直流主导；</item>
    /// <item>归一化到零延迟：输出 [k] = r[k] / r[0]，因此 [0] 恒为 1。</item>
    /// </list>
    /// <para>
    /// 与 FAMOS 的实测差异约 4e-7（FAMOS 内部用 float32，本库用 double）。
    /// </para>
    /// </remarks>
    public static class Spectrum
    {
        /// <summary>
        /// 归一化自相关函数，算法与 imc FAMOS 的 <c>Spec.ACF</c> 一致。
        /// </summary>
        /// <param name="samples">输入样本。</param>
        /// <returns>
        /// 长度为「输入长度以下最大的 2 的幂」的数组，<c>[0]</c> 为 1。
        /// </returns>
        public static double[] Autocorrelation(double[] samples)
        {
            if (samples == null)
            {
                throw new ArgumentNullException(nameof(samples));
            }
            return Autocorrelation(samples, LargestPowerOfTwo(samples.Length));
        }

        /// <summary>
        /// 归一化自相关函数，指定参与计算的样本数（必须是 2 的幂）。
        /// </summary>
        /// <param name="samples">输入样本。</param>
        /// <param name="length">参与计算的样本数，必须是 2 的幂且不超过样本长度。</param>
        public static double[] Autocorrelation(double[] samples, int length)
        {
            if (samples == null)
            {
                throw new ArgumentNullException(nameof(samples));
            }
            if (length < 2)
            {
                throw new ArgumentOutOfRangeException(nameof(length), "长度必须不小于 2。");
            }
            if (length > samples.Length)
            {
                throw new ArgumentException(
                    "参与计算的样本数（" + length + "）不能超过输入长度（" + samples.Length + "）。",
                    nameof(length));
            }
            if (!IsPowerOfTwo(length))
            {
                throw new ArgumentException("长度必须是 2 的幂。", nameof(length));
            }

            int n = length;

            // 1) 前向 FFT
            Complex[] spectrum = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                spectrum[i] = new Complex(samples[i], 0.0);
            }
            Fft(spectrum);

            // 2) 功率谱 |X|²
            double[] power = new double[n];
            for (int i = 0; i < n; i++)
            {
                power[i] = spectrum[i].Real * spectrum[i].Real
                         + spectrum[i].Imaginary * spectrum[i].Imaginary;
            }

            // 3) 循环自相关。r[k] = (1/N)·Σ P[j]·e^{+2πijk/N} 为实数，
            //    而 FFT(P)[k] = Σ P[j]·e^{-2πijk/N}，故 r[k] = Re(FFT(P)[k]) / N。
            Complex[] again = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                again[i] = new Complex(power[i], 0.0);
            }
            Fft(again);

            double[] result = new double[n];
            double zeroLag = again[0].Real / n;
            if (zeroLag == 0.0)
            {
                // 全零输入：FAMOS 的行为亦无定义，这里返回全 0。
                return result;
            }
            for (int i = 0; i < n; i++)
            {
                result[i] = (again[i].Real / n) / zeroLag;
            }
            return result;
        }

        /// <summary>不超过 <paramref name="value"/> 的最大 2 的幂（最小返回 1）。</summary>
        public static int LargestPowerOfTwo(int value)
        {
            if (value < 1)
            {
                return 1;
            }
            int result = 1;
            while (result <= value / 2)
            {
                result <<= 1;
            }
            return result;
        }

        /// <summary>是否为 2 的幂。</summary>
        public static bool IsPowerOfTwo(int value)
        {
            return value > 0 && (value & (value - 1)) == 0;
        }

        /// <summary>
        /// 原地基 2 快速傅里叶变换（Cooley–Tukey，位反转重排）。
        /// 旋转因子预先按 N 计算成表，避免逐级递推带来的误差累积。
        /// </summary>
        private static void Fft(Complex[] buffer)
        {
            int n = buffer.Length;

            // 位反转重排
            for (int i = 1, j = 0; i < n; i++)
            {
                int bit = n >> 1;
                for (; (j & bit) != 0; bit >>= 1)
                {
                    j ^= bit;
                }
                j ^= bit;
                if (i < j)
                {
                    Complex swap = buffer[i];
                    buffer[i] = buffer[j];
                    buffer[j] = swap;
                }
            }

            // 旋转因子表：w[k] = e^{-2πik/N}，k = 0 .. N/2-1
            Complex[] twiddle = new Complex[n >> 1];
            for (int k = 0; k < twiddle.Length; k++)
            {
                double angle = -2.0 * Math.PI * k / n;
                twiddle[k] = new Complex(Math.Cos(angle), Math.Sin(angle));
            }

            for (int size = 2; size <= n; size <<= 1)
            {
                int halfSize = size >> 1;
                int stride = n / size;
                for (int start = 0; start < n; start += size)
                {
                    for (int k = 0; k < halfSize; k++)
                    {
                        Complex w = twiddle[k * stride];
                        Complex even = buffer[start + k];
                        Complex odd = buffer[start + k + halfSize] * w;
                        buffer[start + k] = even + odd;
                        buffer[start + k + halfSize] = even - odd;
                    }
                }
            }
        }
    }
}
