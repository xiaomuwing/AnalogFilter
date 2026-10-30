using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 复数多项式工具：系数按<strong>降幂</strong>排列（与 numpy.poly / scipy 的 b、a 约定一致），
    /// 即 <c>c[0]*x^n + c[1]*x^(n-1) + ... + c[n]</c>。
    /// </summary>
    internal static class ComplexPolynomial
    {
        /// <summary>由根构造首一多项式（降幂）。等价于 <c>numpy.poly(roots)</c>。</summary>
        public static Complex[] FromRoots(IReadOnlyList<Complex> roots)
        {
            Complex[] coeffs = { Complex.One };
            if (roots == null || roots.Count == 0)
            {
                return coeffs;
            }

            foreach (Complex root in roots)
            {
                Complex[] next = new Complex[coeffs.Length + 1];
                next[0] = coeffs[0];
                for (int k = 1; k < coeffs.Length; k++)
                {
                    next[k] = coeffs[k] - root * coeffs[k - 1];
                }
                next[coeffs.Length] = -root * coeffs[coeffs.Length - 1];
                coeffs = next;
            }

            return coeffs;
        }

        /// <summary>Horner 法求值（降幂系数）。</summary>
        public static Complex Evaluate(IReadOnlyList<Complex> coeffs, Complex x)
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < coeffs.Count; i++)
            {
                result = result * x + coeffs[i];
            }
            return result;
        }

        /// <summary>去掉多项式<strong>前导</strong>零系数，至少保留一项。</summary>
        public static Complex[] TrimLeadingZeros(IReadOnlyList<Complex> coeffs)
        {
            int start = 0;
            while (start < coeffs.Count - 1 && coeffs[start] == Complex.Zero)
            {
                start++;
            }

            Complex[] result = new Complex[coeffs.Count - start];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = coeffs[start + i];
            }
            return result;
        }

        /// <summary>
        /// 求多项式全部复根（Aberth–Ehrlich 迭代）。
        /// 对首一化的降幂系数工作，适用于滤波器设计中阶数不太高（N ≤ 24）的多项式。
        /// </summary>
        public static Complex[] Roots(IReadOnlyList<Complex> coeffs)
        {
            Complex[] c = TrimLeadingZeros(coeffs);
            int degree = c.Length - 1;
            if (degree <= 0)
            {
                return Array.Empty<Complex>();
            }

            // 化为首一多项式，避免量级差异影响迭代。
            Complex leading = c[0];
            if (leading == Complex.Zero)
            {
                return Array.Empty<Complex>();
            }

            Complex[] monic = new Complex[c.Length];
            for (int i = 0; i < c.Length; i++)
            {
                monic[i] = c[i] / leading;
            }

            if (degree == 1)
            {
                return new[] { -monic[1] };
            }

            // 初值：以 Cauchy 界为半径、均匀分布在圆上，并用黄金角错开对称性。
            double radius = 1.0;
            for (int i = 1; i < monic.Length; i++)
            {
                radius = Math.Max(radius, 1.0 + Complex.Abs(monic[i]));
            }

            Complex[] roots = new Complex[degree];
            const double goldenAngle = 2.39996322972865332;
            for (int i = 0; i < degree; i++)
            {
                double angle = i * goldenAngle;
                roots[i] = new Complex(radius * Math.Cos(angle), radius * Math.Sin(angle));
            }

            Complex[] derivative = Derivative(monic);
            double tolerance = 1e-14;
            int maxIterations = 2000;

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                double maxStep = 0.0;
                for (int i = 0; i < degree; i++)
                {
                    Complex p = Evaluate(monic, roots[i]);
                    Complex dp = Evaluate(derivative, roots[i]);

                    if (Complex.Abs(dp) < 1e-300)
                    {
                        // 导数为零：轻推一下再继续。
                        roots[i] += new Complex(1e-6 * radius, 1e-6 * radius);
                        continue;
                    }

                    Complex ratio = p / dp;
                    Complex sum = Complex.Zero;
                    for (int j = 0; j < degree; j++)
                    {
                        if (j == i)
                        {
                            continue;
                        }

                        Complex diff = roots[i] - roots[j];
                        if (Complex.Abs(diff) < 1e-300)
                        {
                            diff = new Complex(1e-12, -1e-12);
                        }
                        sum += Complex.One / diff;
                    }

                    Complex denominator = Complex.One - ratio * sum;
                    Complex step = Complex.Abs(denominator) < 1e-300 ? ratio : ratio / denominator;
                    roots[i] -= step;

                    double stepMagnitude = Complex.Abs(step) / (1.0 + Complex.Abs(roots[i]));
                    if (stepMagnitude > maxStep)
                    {
                        maxStep = stepMagnitude;
                    }
                }

                if (maxStep < tolerance)
                {
                    return PolishAndSort(monic, derivative, roots, 2);
                }
            }

            // 未完全收敛时做最后一轮牛顿打磨，尽量给出可用结果。
            return PolishAndSort(monic, derivative, roots, 3);
        }

        private static Complex[] PolishAndSort(Complex[] monic, Complex[] derivative, Complex[] roots, int rounds)
        {
            for (int round = 0; round < rounds; round++)
            {
                for (int i = 0; i < roots.Length; i++)
                {
                    Complex dp = Evaluate(derivative, roots[i]);
                    if (Complex.Abs(dp) < 1e-300)
                    {
                        continue;
                    }
                    roots[i] -= Evaluate(monic, roots[i]) / dp;
                }
            }

            // 归一到实轴附近，消除极小的虚部残留（真实滤波器极点成共轭对）。
            for (int i = 0; i < roots.Length; i++)
            {
                if (Math.Abs(roots[i].Imaginary) < 1e-12 * (1.0 + Complex.Abs(roots[i])))
                {
                    roots[i] = new Complex(roots[i].Real, 0.0);
                }
            }

            Array.Sort(roots, (a, b) =>
            {
                int cmp = a.Real.CompareTo(b.Real);
                return cmp != 0 ? cmp : a.Imaginary.CompareTo(b.Imaginary);
            });
            return roots;
        }

        private static Complex[] Derivative(IReadOnlyList<Complex> coeffs)
        {
            int degree = coeffs.Count - 1;
            if (degree <= 0)
            {
                return new[] { Complex.Zero };
            }

            Complex[] result = new Complex[degree];
            for (int i = 0; i < degree; i++)
            {
                result[i] = coeffs[i] * (degree - i);
            }
            return result;
        }
    }
}
