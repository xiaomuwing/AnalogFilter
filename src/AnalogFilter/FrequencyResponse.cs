using System;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 频率响应计算，对应 <c>scipy.signal.freqs</c>（模拟）与 <c>scipy.signal.freqz</c>（数字）。
    /// 常用于设计结果的验证与绘图。
    /// </summary>
    public static class FrequencyResponse
    {
        /// <summary>模拟滤波器的频率响应 <c>H(jw)</c>，w 单位 rad/s。</summary>
        public static Complex[] Freqs(ZpkModel zpk, double[] w)
        {
            if (zpk == null)
            {
                throw new ArgumentNullException(nameof(zpk));
            }

            Complex[] result = new Complex[w.Length];
            for (int i = 0; i < w.Length; i++)
            {
                Complex s = new Complex(0.0, w[i]);
                Complex numerator = new Complex(zpk.Gain, 0.0);
                for (int k = 0; k < zpk.Zeros.Length; k++)
                {
                    numerator *= s - zpk.Zeros[k];
                }
                Complex denominator = Complex.One;
                for (int k = 0; k < zpk.Poles.Length; k++)
                {
                    denominator *= s - zpk.Poles[k];
                }
                result[i] = numerator / denominator;
            }
            return result;
        }

        /// <summary>模拟滤波器的频率响应 <c>H(jw)</c>，由分子/分母多项式给出。</summary>
        public static Complex[] Freqs(double[] b, double[] a, double[] w)
        {
            Complex[] result = new Complex[w.Length];
            for (int i = 0; i < w.Length; i++)
            {
                Complex s = new Complex(0.0, w[i]);
                result[i] = EvaluateDescending(b, s) / EvaluateDescending(a, s);
            }
            return result;
        }

        /// <summary>数字滤波器的频率响应 <c>H(e^jw)</c>，w 单位 rad/sample。</summary>
        public static Complex[] Freqz(ZpkModel zpk, double[] w)
        {
            if (zpk == null)
            {
                throw new ArgumentNullException(nameof(zpk));
            }

            Complex[] result = new Complex[w.Length];
            for (int i = 0; i < w.Length; i++)
            {
                Complex z = Complex.Exp(new Complex(0.0, w[i]));
                Complex numerator = new Complex(zpk.Gain, 0.0);
                for (int k = 0; k < zpk.Zeros.Length; k++)
                {
                    numerator *= z - zpk.Zeros[k];
                }
                Complex denominator = Complex.One;
                for (int k = 0; k < zpk.Poles.Length; k++)
                {
                    denominator *= z - zpk.Poles[k];
                }
                result[i] = numerator / denominator;
            }
            return result;
        }

        /// <summary>数字滤波器的频率响应 <c>H(e^jw)</c>，由 SOS 级联给出。</summary>
        public static Complex[] Freqz(double[][] sos, double[] w)
        {
            SosFiltering.ValidateSos(sos);

            Complex[] result = new Complex[w.Length];
            for (int i = 0; i < w.Length; i++)
            {
                Complex z = Complex.Exp(new Complex(0.0, w[i]));
                Complex total = Complex.One;
                for (int s = 0; s < sos.Length; s++)
                {
                    double[] section = sos[s];
                    Complex numerator = section[0] + section[1] / z + section[2] / (z * z);
                    Complex denominator = section[3] + section[4] / z + section[5] / (z * z);
                    total *= numerator / denominator;
                }
                result[i] = total;
            }
            return result;
        }

        /// <summary>数字滤波器的频率响应 <c>H(e^jw)</c>，由分子/分母多项式给出。</summary>
        public static Complex[] Freqz(double[] b, double[] a, double[] w)
        {
            Complex[] result = new Complex[w.Length];
            for (int i = 0; i < w.Length; i++)
            {
                Complex z = Complex.Exp(new Complex(0.0, w[i]));
                result[i] = EvaluateDescending(b, z) / EvaluateDescending(a, z);
            }
            return result;
        }

        /// <summary>幅频响应（线性幅值）。</summary>
        public static double[] Magnitude(Complex[] response)
        {
            double[] result = new double[response.Length];
            for (int i = 0; i < response.Length; i++)
            {
                result[i] = Complex.Abs(response[i]);
            }
            return result;
        }

        /// <summary>幅频响应（dB，20log10）。</summary>
        public static double[] MagnitudeDb(Complex[] response)
        {
            double[] result = new double[response.Length];
            for (int i = 0; i < response.Length; i++)
            {
                double magnitude = Complex.Abs(response[i]);
                result[i] = 20.0 * Math.Log10(Math.Max(magnitude, double.Epsilon));
            }
            return result;
        }

        private static Complex EvaluateDescending(double[] coefficients, Complex x)
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < coefficients.Length; i++)
            {
                result = result * x + coefficients[i];
            }
            return result;
        }
    }
}
