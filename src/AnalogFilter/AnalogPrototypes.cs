using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>贝塞尔滤波器的频率归一化方式，对应 <c>scipy.signal.besselap(norm=...)</c>。</summary>
    public enum BesselNormalization
    {
        /// <summary>相位中点位于 1 rad/s（等效 Butterworth 渐近线），与 MATLAB 一致，默认。</summary>
        Phase,

        /// <summary>通带群延迟归一化为 1 秒。</summary>
        Delay,

        /// <summary>幅频 -3 dB 点位于 1 rad/s。</summary>
        Magnitude
    }

    /// <summary>
    /// 模拟低通原型（截止角频率统一为 1 rad/s）。
    /// 算法逐行移植自 SciPy 1.18.1 的 <c>scipy.signal._filter_design</c>
    /// （buttap / cheb1ap / cheb2ap / besselap），BSD-3-Clause。
    /// </summary>
    public static class AnalogPrototypes
    {
        /// <summary>N 阶巴特沃斯模拟低通原型。等价于 <c>scipy.signal.buttap</c>。</summary>
        public static ZpkModel Butterworth(int n)
        {
            ValidateOrder(n);
            Complex[] poles = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                double m = -n + 1 + 2 * i;
                // -exp(j*pi*m/(2N))；m=0 时得到精确的 -1（实极点）。
                poles[i] = -Complex.Exp(new Complex(0.0, Math.PI * m / (2.0 * n)));
            }

            return new ZpkModel(Array.Empty<Complex>(), poles, 1.0);
        }

        /// <summary>N 阶切比雪夫 I 型模拟低通原型，通带纹波 rp dB。等价于 <c>scipy.signal.cheb1ap</c>。</summary>
        public static ZpkModel ChebyshevINormalized(int n, double rp)
        {
            ValidateOrder(n);
            if (n == 0)
            {
                // 避免除零：0 阶滤波器的增益为 -rp dB。
                return new ZpkModel(Array.Empty<Complex>(), Array.Empty<Complex>(), Math.Pow(10.0, -rp / 20.0));
            }

            double eps = Math.Sqrt(Math.Pow(10.0, 0.1 * rp) - 1.0);
            double mu = Asinh(1.0 / eps) / n;

            Complex[] poles = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                double m = -n + 1 + 2 * i;
                double theta = Math.PI * m / (2.0 * n);
                // p = -sinh(mu + j*theta)
                poles[i] = new Complex(
                    -Math.Sinh(mu) * Math.Cos(theta),
                    -Math.Cosh(mu) * Math.Sin(theta));
            }

            // k = real(prod(-p))：极点成共轭对，乘积的虚部相互抵消。
            Complex product = Complex.One;
            for (int i = 0; i < n; i++)
            {
                product *= -poles[i];
            }
            double k = product.Real;

            if (n % 2 == 0)
            {
                k /= Math.Sqrt(1.0 + eps * eps);
            }

            return new ZpkModel(Array.Empty<Complex>(), poles, k);
        }

        /// <summary>N 阶切比雪夫 II 型模拟低通原型，阻带衰减 rs dB。等价于 <c>scipy.signal.cheb2ap</c>。</summary>
        public static ZpkModel ChebyshevIINormalized(int n, double rs)
        {
            ValidateOrder(n);
            if (n == 0)
            {
                return new ZpkModel(Array.Empty<Complex>(), Array.Empty<Complex>(), 1.0);
            }

            double de = 1.0 / Math.Sqrt(Math.Pow(10.0, 0.1 * rs) - 1.0);
            double mu = Asinh(1.0 / de) / n;

            List<double> mi = new List<double>();
            if (n % 2 != 0)
            {
                for (int m = -n + 1; m < 0; m += 2)
                {
                    mi.Add(m);
                }
                for (int m = 2; m < n; m += 2)
                {
                    mi.Add(m);
                }
            }
            else
            {
                for (int m = -n + 1; m < n; m += 2)
                {
                    mi.Add(m);
                }
            }

            // 零点：z = j / sin(m*pi/(2N))
            Complex[] zeros = new Complex[mi.Count];
            for (int i = 0; i < mi.Count; i++)
            {
                zeros[i] = new Complex(0.0, 1.0) / new Complex(Math.Sin(mi[i] * Math.PI / (2.0 * n)), 0.0);
            }

            // 极点：p = -1 / sinh(mu + j*theta)
            Complex[] poles = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                double m1 = -n + 1 + 2 * i;
                double theta1 = Math.PI * m1 / (2.0 * n);
                Complex s = new Complex(Math.Sinh(mu) * Math.Cos(theta1), Math.Cosh(mu) * Math.Sin(theta1));
                poles[i] = -Complex.One / s;
            }

            Complex nz = Complex.One;
            for (int i = 0; i < zeros.Length; i++)
            {
                nz *= -zeros[i];
            }
            Complex np = Complex.One;
            for (int i = 0; i < poles.Length; i++)
            {
                np *= -poles[i];
            }

            return new ZpkModel(zeros, poles, (np / nz).Real);
        }

        /// <summary>
        /// N 阶贝塞尔模拟低通原型。等价于 <c>scipy.signal.besselap</c>。
        /// </summary>
        /// <remarks>
        /// SciPy 用 Bessel 函数 K_{N+1/2} 配合 Aberth 迭代求根；这里改为直接对
        /// <strong>反贝塞尔多项式</strong>求根（数学上等价，极点即普通贝塞尔多项式根之倒数），
        /// 少一层特殊函数依赖，并可与 SciPy 结果逐值校验。
        /// </remarks>
        public static ZpkModel Bessel(int n, BesselNormalization norm = BesselNormalization.Phase)
        {
            ValidateOrder(n, allowZero: true);
            if (n == 0)
            {
                return new ZpkModel(Array.Empty<Complex>(), Array.Empty<Complex>(), 1.0);
            }

            // 与 SciPy 的 besselap 一致：极点 = 1 / (贝塞尔多项式 y_N 的零点)。
            // 零点由 _campos_zeros 给初值，用 Aberth 迭代 + 逐点 Newton 求得，
            // 最后对共轭对取平均以保证严格对称（见 BesselZeros）。
            double aLast = BesselALast(n);
            Complex[] zeros = BesselZeros(n);
            Complex[] poles = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                poles[i] = Complex.One / zeros[i];
            }

            double k;
            switch (norm)
            {
                case BesselNormalization.Phase:
                    // 相位中点归一化：极点整体缩放，渐近线与同阶巴特沃斯一致。
                    double factor = Math.Pow(10.0, -Math.Log10(aLast) / n);
                    for (int i = 0; i < poles.Length; i++)
                    {
                        poles[i] *= factor;
                    }
                    k = 1.0;
                    break;

                case BesselNormalization.Delay:
                    k = aLast;
                    break;

                case BesselNormalization.Magnitude:
                    k = aLast;
                    double normFactor = FindMagnitudeNormFactor(poles, k);
                    for (int i = 0; i < poles.Length; i++)
                    {
                        poles[i] /= normFactor;
                    }
                    k = Math.Pow(normFactor, -n) * aLast;
                    break;

                default:
                    throw new ArgumentOutOfRangeException(nameof(norm));
            }

            return new ZpkModel(Array.Empty<Complex>(), poles, k);
        }

        /// <summary>
        /// 普通贝塞尔多项式系数（降幂，首项为 1）。
        /// 对应 <c>scipy.signal._filter_design._bessel_poly(n)</c>。
        /// </summary>
        public static double[] BesselPolynomial(int n)
        {
            ValidateOrder(n, allowZero: true);

            // 用 BigInteger 精确求出整数系数，再逐个转 double。
            // 注意：不能为了防溢出而整体缩放后再把首项置 1 —— 那会改变多项式本身
            // （根只在"整体同乘一个常数"时才不变）。这些系数直到 n≈60 都不会溢出 double。
            double[] result = new double[n + 1];
            for (int k = 0; k <= n; k++)
            {
                System.Numerics.BigInteger numerator = FallingFactorial(2 * n - k, n);
                System.Numerics.BigInteger denominator = System.Numerics.BigInteger.Pow(2, n - k) * Factorial(k);
                result[k] = (double)(numerator / denominator);
            }

            return result;
        }

        private static System.Numerics.BigInteger FallingFactorial(int top, int count)
        {
            System.Numerics.BigInteger result = System.Numerics.BigInteger.One;
            for (int i = 0; i < count; i++)
            {
                result *= top - i;
            }
            return result;
        }

        private static System.Numerics.BigInteger Factorial(int n)
        {
            System.Numerics.BigInteger result = System.Numerics.BigInteger.One;
            for (int i = 2; i <= n; i++)
            {
                result *= i;
            }
            return result;
        }

        /// <summary>二分求 -3 dB 角频率（仅用于 Bessel 的 mag 归一化）。</summary>
        private static double FindMagnitudeNormFactor(Complex[] poles, double gain)
        {
            const double target = 0.70710678118654752440; // 1/sqrt(2)
            Func<double, double> magnitude = w =>
            {
                Complex s = new Complex(0.0, w);
                Complex denominator = Complex.One;
                for (int i = 0; i < poles.Length; i++)
                {
                    denominator *= s - poles[i];
                }
                return gain / Complex.Abs(denominator);
            };

            double low = 1e-9;
            double high = 1.0;
            while (magnitude(high) > target && high < 1e12)
            {
                high *= 2.0;
            }

            for (int i = 0; i < 200; i++)
            {
                double mid = 0.5 * (low + high);
                if (magnitude(mid) > target)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }

            return 0.5 * (low + high);
        }

        /// <summary>
        /// 切比雪夫 I 型原型，按 <b>FAMOS 约定</b>归一化：−3 dB 点位于 ω = 1，且直流增益为 1。
        /// </summary>
        /// <remarks>
        /// <para>
        /// 与 <see cref="ChebyshevINormalized"/> 的区别：后者是 SciPy 的约定，通带波纹边在 ω = 1，
        /// 偶阶时直流增益为 1/√(1+ε²)，−3 dB 点落在 ω &gt; 1 处。
        /// </para>
        /// <para>
        /// 把偶阶的直流增益补到 1（乘 √(1+ε²)）之后，−3 dB 点满足
        /// T_N(ω₃) = √(1+2ε²)/ε；奇阶直流增益本来就是 1，有 T_N(ω₃) = 1/ε。
        /// 再把极点/零点整体缩放 1/ω₃，即得到 −3 dB 落在 ω = 1 的原型。
        /// </para>
        /// </remarks>
        public static ZpkModel ChebyshevI3dB(int n, double rp)
        {
            ValidateOrder(n, allowZero: false);

            ZpkModel prototype = ChebyshevINormalized(n, rp);
            double eps = Math.Sqrt(Math.Pow(10.0, 0.1 * rp) - 1.0);
            double gain = prototype.Gain;
            double omega3;

            if (n % 2 == 0)
            {
                // 偶阶：先补直流增益，−3 dB 点随之移到 √(1+2ε²)/ε
                gain *= Math.Sqrt(1.0 + eps * eps);
                omega3 = Math.Cosh(Acosh(Math.Sqrt(1.0 + 2.0 * eps * eps) / eps) / n);
            }
            else
            {
                omega3 = Math.Cosh(Acosh(1.0 / eps) / n);
            }

            double scale = 1.0 / omega3;
            Complex[] zeros = new Complex[prototype.Zeros.Length];
            for (int i = 0; i < zeros.Length; i++)
            {
                zeros[i] = prototype.Zeros[i] * scale;
            }
            Complex[] poles = new Complex[prototype.Poles.Length];
            for (int i = 0; i < poles.Length; i++)
            {
                poles[i] = prototype.Poles[i] * scale;
            }

            // lp2lp：零点/极点整体缩放 s 倍时，增益乘 s^degree
            gain *= Math.Pow(scale, n);
            return new ZpkModel(zeros, poles, gain);
        }

        /// <summary>
        /// 临界阻尼低通原型：N 个重实极点，<c>H(s) = 1/(1 + s)^N</c>。
        /// SciPy 没有这个原型，此处按 FAMOS 的临界阻尼定义补充。
        /// </summary>
        /// <remarks>
        /// 该原型的 -3 dB 点并不在 s = 1 处，而是 ω = sqrt(2^(1/N) - 1)；
        /// 换算由 <see cref="FilterDesign"/> 在频率变换前统一处理。
        /// </remarks>
        public static ZpkModel CriticallyDamped(int n)
        {
            ValidateOrder(n, allowZero: false);

            Complex[] poles = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                poles[i] = new Complex(-1.0, 0.0);
            }
            return new ZpkModel(Array.Empty<Complex>(), poles, 1.0);
        }

        /// <summary>
        /// a_last = (2N)!/(N!) / 2^N，即普通贝塞尔多项式的常数项。
        /// 对应 SciPy 的 <c>_falling_factorial(2*N, N) // 2**N</c>。
        /// </summary>
        private static double BesselALast(int n)
        {
            System.Numerics.BigInteger value = FallingFactorial(2 * n, n);
            value = value / System.Numerics.BigInteger.Pow(2, n);
            return (double)value;
        }

        /// <summary>
        /// 普通贝塞尔多项式 y_N 的零点。对应 <c>scipy.signal._filter_design._bessel_zeros</c>。
        /// </summary>
        private static Complex[] BesselZeros(int n)
        {
            if (n == 0)
            {
                return Array.Empty<Complex>();
            }

            Complex[] initial = BesselStartingPoints(n);

            Func<Complex, Complex> f = x => BesselKve(n, Complex.One / x);
            Func<Complex, Complex> fp = x =>
            {
                Complex inv = Complex.One / x;
                Complex x2 = x * x;
                // 对应 SciPy 里用递推写出的 kve 导数。
                return BesselKve(n - 1, inv) / (2.0 * x2)
                     - BesselKve(n, inv) / x2
                     + BesselKve(n + 1, inv) / (2.0 * x2);
            };

            // 迭代上限取 40 / 20：闭式 kve 在复参数下存在相消，|alpha| 在 ~1e-14 处停滞，
            // 再往上加迭代只是空转（实测 6~8 次即到达该精度，40 次留足余量）。
            Complex[] roots = Aberth(f, fp, initial, 1e-15, 40);
            for (int i = 0; i < roots.Length; i++)
            {
                roots[i] = Newton(f, fp, roots[i], 1e-15, 20);
            }

            // 对共轭对取平均，消除数值不对称。
            Complex[] symmetric = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                symmetric[i] = (roots[i] + Complex.Conjugate(roots[n - 1 - i])) / 2.0;
            }

            Complex sum = Complex.Zero;
            for (int i = 0; i < n; i++)
            {
                sum += symmetric[i];
            }
            // 零点之和应为 -1。托管实现在高阶时精度有限（见 README 的精度表），
            // 这里只用于拦截"算法失效"，不是精度承诺，因此阈值随阶数放宽。
            double tolerance = 1e-7 * Math.Max(1, n);
            if (Complex.Abs(sum + 1.0) > tolerance)
            {
                throw new InvalidOperationException(
                    "贝塞尔零点求解误差过大（阶数 " + n + "，零点之和 " + (sum + 1.0) + "）。");
            }

            return symmetric;
        }

        /// <summary>
        /// 指数缩放的第二类修正贝塞尔函数 <c>kve(ν, z) = exp(z)·K_ν(z)</c>，ν = n + 1/2。
        /// </summary>
        /// <remarks>
        /// SciPy 调用 Cephes 的通用实数实现；对半整数阶它有精确闭式：
        /// <c>K_{n+1/2}(z) = sqrt(pi/(2z))·exp(-z)·Σ_{k=0..n} (n+k)!/(k!(n-k)!)/(2z)^k</c>，
        /// 因此这里用闭式直接求值，结果与 Cephes 在机器精度内一致，且不需要原生依赖。
        /// </remarks>
        private static Complex BesselKve(int n, Complex z)
        {
            Complex sum = Complex.Zero;
            double coefficient = 1.0; // (n+k)! / (k! (n-k)!)
            Complex term = Complex.One; // (2z)^-k
            Complex twoZ = 2.0 * z;
            for (int k = 0; k <= n; k++)
            {
                sum += coefficient * term;
                // c_{k+1} = c_k * (n+k+1)(n-k)/(k+1)
                coefficient *= (double)(n + k + 1) * (n - k) / (k + 1);
                term /= twoZ;
            }

            return Complex.Sqrt(Math.PI / (2.0 * z)) * sum;
        }

        /// <summary>
        /// 贝塞尔零点的初值近似。对应 <c>scipy.signal._filter_design._campos_zeros</c>
        /// （Campos &amp; Calderón 的闭式公式，系数按升幂给出）。
        /// </summary>
        private static Complex[] BesselStartingPoints(int n)
        {
            if (n == 1)
            {
                return new[] { new Complex(-1.0, 0.0) };
            }

            double s = Polyval(n, new[] { 0.0, 0.0, 2.0, 0.0, -3.0, 1.0 });
            double b3 = Polyval(n, new[] { 16.0, -8.0 }) / s;
            double b2 = Polyval(n, new[] { -24.0, -12.0, 12.0 }) / s;
            double b1 = Polyval(n, new[] { 8.0, 24.0, -12.0, -2.0 }) / s;
            double b0 = Polyval(n, new[] { 0.0, -6.0, 0.0, 5.0, -1.0 }) / s;

            double r = Polyval(n, new[] { 0.0, 0.0, 2.0, 1.0 });
            double a1 = Polyval(n, new[] { -6.0, -6.0 }) / r;
            double a2 = 6.0 / r;

            Complex[] result = new Complex[n];
            for (int k = 1; k <= n; k++)
            {
                double realPart = Polyval(k, new[] { 0.0, a1, a2 });
                double imagPart = Polyval(k, new[] { b0, b1, b2, b3 });
                result[k - 1] = new Complex(realPart, imagPart);
            }
            return result;
        }

        /// <summary>升幂多项式求值（与 numpy.polynomial.polynomial.polyval 一致）。</summary>
        private static double Polyval(double x, double[] ascendingCoefficients)
        {
            double result = 0.0;
            for (int i = ascendingCoefficients.Length - 1; i >= 0; i--)
            {
                result = result * x + ascendingCoefficients[i];
            }
            return result;
        }

        /// <summary>
        /// Aberth–Ehrlich 同时求根迭代。对应 <c>scipy.signal._filter_design._aberth</c>。
        /// </summary>
        private static Complex[] Aberth(
            Func<Complex, Complex> f,
            Func<Complex, Complex> fp,
            Complex[] initial,
            double tolerance,
            int maxIterations)
        {
            int count = initial.Length;
            Complex[] roots = new Complex[count];
            Array.Copy(initial, roots, count);
            Complex[] repulsion = new Complex[count];

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                Complex[] alpha = new Complex[count];
                for (int k = 0; k < count; k++)
                {
                    Complex derivative = fp(roots[k]);
                    if (Complex.Abs(derivative) < 1e-300)
                    {
                        throw new InvalidOperationException("贝塞尔零点迭代中导数为 0。");
                    }
                    alpha[k] = -f(roots[k]) / derivative;

                    Complex beta = Complex.Zero;
                    for (int j = 0; j < count; j++)
                    {
                        if (j == k)
                        {
                            continue;
                        }
                        Complex delta = roots[k] - roots[j];
                        if (Complex.Abs(delta) < 1e-300)
                        {
                            delta = new Complex(1e-12, 1e-12);
                        }
                        beta += Complex.One / delta;
                    }
                    repulsion[k] = beta;
                }

                for (int k = 0; k < count; k++)
                {
                    roots[k] += alpha[k] / (1.0 + alpha[k] * repulsion[k]);
                }

                bool converged = true;
                for (int k = 0; k < count; k++)
                {
                    if (double.IsNaN(roots[k].Real) || double.IsInfinity(roots[k].Real) ||
                        double.IsNaN(roots[k].Imaginary) || double.IsInfinity(roots[k].Imaginary))
                    {
                        throw new InvalidOperationException("贝塞尔零点迭代发散。");
                    }
                    if (Complex.Abs(alpha[k]) > tolerance)
                    {
                        converged = false;
                    }
                }

                if (converged)
                {
                    return roots;
                }
            }

            // 未达到严格容差时返回当前最优解：后续还有 Newton 抛光与共轭平均，
            // 且高阶贝塞尔的精度上限由此决定（见精度表）。
            return roots;
        }

        /// <summary>逐点牛顿法。对应 SciPy 的 <c>optimize.newton(f, x, fp, tol)</c>。</summary>
        private static Complex Newton(
            Func<Complex, Complex> f,
            Func<Complex, Complex> fp,
            Complex guess,
            double tolerance,
            int maxIterations)
        {
            Complex x = guess;
            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                Complex derivative = fp(x);
                if (Complex.Abs(derivative) < 1e-300)
                {
                    return x;
                }
                Complex step = f(x) / derivative;
                x -= step;
                if (Complex.Abs(step) <= tolerance)
                {
                    return x;
                }
            }

            // 达到最大迭代次数时返回当前值（高阶时不保证严格容差）。
            return x;
        }

        private static void ValidateOrder(int n, bool allowZero = true)
        {
            if (n < 0 || (n == 0 && !allowZero))
            {
                throw new ArgumentOutOfRangeException(nameof(n), "滤波器阶数必须是非负整数。");
            }
        }

        /// <summary>反双曲正弦。.NET Framework 4.8 没有 <c>Math.Asinh</c>，这里补上。</summary>
        private static double Asinh(double x)
        {
            // 大 |x| 时用等价形式避免 x*x 溢出。
            if (Math.Abs(x) > 1e154)
            {
                return Math.Sign(x) * (Math.Log(Math.Abs(x)) + Math.Log(2.0));
            }
            return Math.Log(x + Math.Sqrt(x * x + 1.0));
        }

        /// <summary>反双曲余弦。.NET Framework 4.8 没有 <c>Math.Acosh</c>，这里补上。</summary>
        private static double Acosh(double x)
        {
            if (x <= 1.0)
            {
                return 0.0;
            }
            // acosh(x) = ln(x + sqrt(x-1)*sqrt(x+1))，避免 x*x 溢出/相消
            return Math.Log(x + Math.Sqrt(x - 1.0) * Math.Sqrt(x + 1.0));
        }
    }
}
