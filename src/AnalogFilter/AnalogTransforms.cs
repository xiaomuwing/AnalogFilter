using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 模拟域的频率变换（低通原型 → 低通/高通/带通/带阻）。
    /// 逐行移植自 SciPy 1.18.1 的 <c>lp2lp_zpk</c> / <c>lp2hp_zpk</c> /
    /// <c>lp2bp_zpk</c> / <c>lp2bs_zpk</c>，BSD-3-Clause。
    /// 所有频率都是<strong>角频率</strong>（rad/s）；数字滤波器设计时由 <see cref="BilinearTransform"/>
    /// 负责把角频率与采样率联系起来。
    /// </summary>
    public static class AnalogTransforms
    {
        /// <summary>低通 → 低通，截止角频率 wo。等价于 <c>scipy.signal.lp2lp_zpk</c>。</summary>
        public static ZpkModel LowPassToLowPass(ZpkModel prototype, double wo)
        {
            if (prototype == null)
            {
                throw new ArgumentNullException(nameof(prototype));
            }
            ValidatePositive(wo, nameof(wo));

            int degree = RelativeDegree(prototype);
            Complex[] zeros = Scale(prototype.Zeros, wo);
            Complex[] poles = Scale(prototype.Poles, wo);
            double gain = prototype.Gain * Math.Pow(wo, degree);
            return new ZpkModel(zeros, poles, gain);
        }

        /// <summary>低通 → 高通，截止角频率 wo。等价于 <c>scipy.signal.lp2hp_zpk</c>。</summary>
        public static ZpkModel LowPassToHighPass(ZpkModel prototype, double wo)
        {
            if (prototype == null)
            {
                throw new ArgumentNullException(nameof(prototype));
            }
            ValidatePositive(wo, nameof(wo));

            int degree = RelativeDegree(prototype);
            Complex[] zeros = new Complex[prototype.Zeros.Length + degree];
            for (int i = 0; i < prototype.Zeros.Length; i++)
            {
                zeros[i] = wo / prototype.Zeros[i];
            }
            // 低通的原零点在无穷远，反演后落到原点。
            Complex[] poles = new Complex[prototype.Poles.Length];
            for (int i = 0; i < prototype.Poles.Length; i++)
            {
                poles[i] = wo / prototype.Poles[i];
            }

            double gain = prototype.Gain * (ProdNegative(prototype.Zeros) / ProdNegative(prototype.Poles)).Real;
            return new ZpkModel(zeros, poles, gain);
        }

        /// <summary>
        /// 低通 → 带通，中心角频率 wo、带宽 bw（均为 rad/s）。
        /// 等价于 <c>scipy.signal.lp2bp_zpk</c>。
        /// </summary>
        public static ZpkModel LowPassToBandPass(ZpkModel prototype, double wo, double bw)
        {
            if (prototype == null)
            {
                throw new ArgumentNullException(nameof(prototype));
            }
            ValidatePositive(wo, nameof(wo));
            ValidatePositive(bw, nameof(bw));

            int degree = RelativeDegree(prototype);
            double woSquared = wo * wo;

            Complex[] bpZeros = SplitBySqrt(Scale(prototype.Zeros, bw / 2.0), woSquared);
            Complex[] bpPoles = SplitBySqrt(Scale(prototype.Poles, bw / 2.0), woSquared);

            // 低通的 degree 个无穷远零点，搬移后落在原点。
            Complex[] zeros = new Complex[bpZeros.Length + degree];
            Array.Copy(bpZeros, zeros, bpZeros.Length);

            double gain = prototype.Gain * Math.Pow(bw, degree);
            return new ZpkModel(zeros, bpPoles, gain);
        }

        /// <summary>
        /// 低通 → 带阻，中心角频率 wo、带宽 bw（均为 rad/s）。
        /// 等价于 <c>scipy.signal.lp2bs_zpk</c>。
        /// </summary>
        public static ZpkModel LowPassToBandStop(ZpkModel prototype, double wo, double bw)
        {
            if (prototype == null)
            {
                throw new ArgumentNullException(nameof(prototype));
            }
            ValidatePositive(wo, nameof(wo));
            ValidatePositive(bw, nameof(bw));

            int degree = RelativeDegree(prototype);
            double woSquared = wo * wo;

            Complex[] inverseZeros = new Complex[prototype.Zeros.Length];
            for (int i = 0; i < prototype.Zeros.Length; i++)
            {
                inverseZeros[i] = (bw / 2.0) / prototype.Zeros[i];
            }
            Complex[] inversePoles = new Complex[prototype.Poles.Length];
            for (int i = 0; i < prototype.Poles.Length; i++)
            {
                inversePoles[i] = (bw / 2.0) / prototype.Poles[i];
            }

            Complex[] bsZeros = SplitBySqrt(inverseZeros, woSquared);
            Complex[] bsPoles = SplitBySqrt(inversePoles, woSquared);

            // 无穷远处的零点搬到阻带中心 ±j*wo。
            Complex[] zeros = new Complex[bsZeros.Length + 2 * degree];
            Array.Copy(bsZeros, zeros, bsZeros.Length);
            for (int i = 0; i < degree; i++)
            {
                zeros[bsZeros.Length + i] = new Complex(0.0, wo);
                zeros[bsZeros.Length + degree + i] = new Complex(0.0, -wo);
            }

            double gain = prototype.Gain * (ProdNegative(prototype.Zeros) / ProdNegative(prototype.Poles)).Real;
            return new ZpkModel(zeros, bsPoles, gain);
        }

        /// <summary>相对阶（极点数 − 零点数），必须 ≥ 0。</summary>
        public static int RelativeDegree(ZpkModel zpk)
        {
            int degree = zpk.Poles.Length - zpk.Zeros.Length;
            if (degree < 0)
            {
                throw new ArgumentException("传递函数非真（零点多于极点），无法进行频率变换。");
            }
            return degree;
        }

        /// <summary>把 z、p 整体缩放到 bw/2 后，按 ±sqrt(...) 生成两路共轭分支。</summary>
        private static Complex[] SplitBySqrt(Complex[] source, double woSquared)
        {
            Complex[] result = new Complex[source.Length * 2];
            for (int i = 0; i < source.Length; i++)
            {
                Complex root = Complex.Sqrt(source[i] * source[i] - woSquared);
                result[i] = source[i] + root;
                result[source.Length + i] = source[i] - root;
            }
            return result;
        }

        private static Complex[] Scale(Complex[] source, double factor)
        {
            Complex[] result = new Complex[source.Length];
            for (int i = 0; i < source.Length; i++)
            {
                result[i] = source[i] * factor;
            }
            return result;
        }

        private static Complex ProdNegative(Complex[] values)
        {
            Complex product = Complex.One;
            for (int i = 0; i < values.Length; i++)
            {
                product *= -values[i];
            }
            return product;
        }

        private static void ValidatePositive(double value, string name)
        {
            if (double.IsNaN(value) || double.IsInfinity(value) || value <= 0.0)
            {
                throw new ArgumentOutOfRangeException(name, "频率必须是大于 0 的有限数值。");
            }
        }
    }
}
