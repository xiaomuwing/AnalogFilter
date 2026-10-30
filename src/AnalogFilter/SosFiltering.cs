using System;

namespace AnalogFilter
{
    /// <summary>
    /// 二阶节（SOS）级联滤波与直接 II 型转置（DF2T）差分方程求解。
    /// 算法与 <c>scipy.signal.sosfilt</c> / <c>scipy.signal.lfilter</c> 一致。
    /// </summary>
    public static class SosFiltering
    {
        /// <summary>
        /// 用 SOS 级联对信号滤波（零初始状态）。等价于 <c>scipy.signal.sosfilt(sos, x)</c>。
        /// </summary>
        public static double[] SosFilt(double[][] sos, double[] x)
        {
            double[][] finalState;
            return SosFilt(sos, x, null, out finalState);
        }

        /// <summary>
        /// 用 SOS 级联对信号滤波，可指定/获取每节的状态（长度 2）。
        /// 等价于 <c>scipy.signal.sosfilt(sos, x, zi=...)</c>。
        /// </summary>
        public static double[] SosFilt(double[][] sos, double[] x, double[][] initialState, out double[][] finalState)
        {
            ValidateSos(sos);
            if (x == null)
            {
                throw new ArgumentNullException(nameof(x));
            }

            int sectionCount = sos.Length;
            finalState = new double[sectionCount][];
            if (initialState == null)
            {
                for (int s = 0; s < sectionCount; s++)
                {
                    finalState[s] = new double[2];
                }
            }
            else
            {
                if (initialState.Length != sectionCount)
                {
                    throw new ArgumentException("初始状态的长度必须等于节数。", nameof(initialState));
                }
                for (int s = 0; s < sectionCount; s++)
                {
                    finalState[s] = new[] { initialState[s][0], initialState[s][1] };
                }
            }

            double[] output = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                double value = x[i];
                for (int s = 0; s < sectionCount; s++)
                {
                    double[] section = sos[s];
                    double[] state = finalState[s];

                    // 直接 II 型转置；a0 已在 ValidateSos 中确认为 1。
                    double y = section[0] * value + state[0];
                    state[0] = section[1] * value - section[4] * y + state[1];
                    state[1] = section[2] * value - section[5] * y;
                    value = y;
                }
                output[i] = value;
            }

            return output;
        }

        /// <summary>
        /// 直接 II 型转置的差分方程求解。等价于 <c>scipy.signal.lfilter(b, a, x)</c>。
        /// </summary>
        public static double[] LFilter(double[] b, double[] a, double[] x)
        {
            double[] finalState;
            return LFilter(b, a, x, null, out finalState);
        }

        /// <summary>
        /// 带初始状态的直接 II 型转置求解，并可取出终态。
        /// 等价于 <c>scipy.signal.lfilter(b, a, x, zi=...)</c>。
        /// </summary>
        public static double[] LFilter(double[] b, double[] a, double[] x, double[] zi, out double[] zf)
        {
            if (b == null || b.Length == 0)
            {
                throw new ArgumentException("分子不能为空。", nameof(b));
            }
            if (a == null || a.Length == 0)
            {
                throw new ArgumentException("分母不能为空。", nameof(a));
            }
            if (x == null)
            {
                throw new ArgumentNullException(nameof(x));
            }
            if (a[0] == 0.0)
            {
                throw new ArgumentException("分母首项不能为 0。", nameof(a));
            }

            int order = Math.Max(b.Length, a.Length) - 1;
            double[] bb = new double[order + 1];
            double[] aa = new double[order + 1];
            Array.Copy(b, bb, b.Length);
            Array.Copy(a, aa, a.Length);

            double a0 = aa[0];
            for (int i = 0; i <= order; i++)
            {
                bb[i] /= a0;
                aa[i] /= a0;
            }

            double[] state = new double[order];
            if (zi != null)
            {
                if (zi.Length != order)
                {
                    throw new ArgumentException("初始状态长度必须为 max(len(a), len(b)) - 1。", nameof(zi));
                }
                Array.Copy(zi, state, order);
            }

            double[] output = new double[x.Length];
            for (int n = 0; n < x.Length; n++)
            {
                double input = x[n];
                double y = bb[0] * input + (order > 0 ? state[0] : 0.0);

                for (int k = 0; k < order - 1; k++)
                {
                    state[k] = bb[k + 1] * input - aa[k + 1] * y + state[k + 1];
                }
                if (order > 0)
                {
                    state[order - 1] = bb[order] * input - aa[order] * y;
                }

                output[n] = y;
            }

            zf = state;
            return output;
        }

        /// <summary>
        /// 单位阶跃响应的稳态初始条件。等价于 <c>scipy.signal.lfilter_zi</c>。
        /// </summary>
        /// <remarks>
        /// SciPy 1.18 用累积和闭式求解（早期版本用线性方程组）；这里沿用闭式，因此无需任何原生线性代数库。
        /// </remarks>
        public static double[] LFilterZi(double[] b, double[] a)
        {
            if (b == null || b.Length == 0)
            {
                throw new ArgumentException("分子不能为空。", nameof(b));
            }
            if (a == null || a.Length == 0)
            {
                throw new ArgumentException("分母不能为空。", nameof(a));
            }
            if (a[0] == 0.0)
            {
                throw new ArgumentException("分母首项不能为 0。", nameof(a));
            }

            int n = Math.Max(a.Length, b.Length);
            double[] bb = new double[n];
            double[] aa = new double[n];
            double a0 = a[0];
            for (int i = 0; i < b.Length; i++)
            {
                bb[i] = b[i] / a0;
            }
            for (int i = 0; i < a.Length; i++)
            {
                aa[i] = a[i] / a0;
            }

            double sumA = 0.0;
            for (int i = 0; i < n; i++)
            {
                sumA += aa[i];
            }
            if (sumA == 0.0)
            {
                throw new ArgumentException("滤波器在 z = 1 处有极点（sum(a) = 0），不稳定。", nameof(a));
            }

            double sumB = 0.0;
            for (int i = 0; i < n; i++)
            {
                sumB += bb[i];
            }
            double yInf = sumB / sumA;

            // flip(cumsum(flip(b - y_inf*a)))[1:]
            double[] cumulative = new double[n];
            double running = 0.0;
            for (int i = n - 1; i >= 0; i--)
            {
                running += bb[i] - yInf * aa[i];
                cumulative[i] = running;
            }

            double[] result = new double[n - 1];
            for (int i = 1; i < n; i++)
            {
                result[i - 1] = cumulative[i];
            }
            return result;
        }

        /// <summary>
        /// SOS 级联的稳态初始条件。等价于 <c>scipy.signal.sosfilt_zi</c>。
        /// </summary>
        public static double[][] SosFiltZi(double[][] sos)
        {
            ValidateSos(sos);

            double[][] result = new double[sos.Length][];
            double scale = 1.0;
            for (int s = 0; s < sos.Length; s++)
            {
                double[] b = { sos[s][0], sos[s][1], sos[s][2] };
                double[] a = { sos[s][3], sos[s][4], sos[s][5] };
                double[] section = LFilterZi(b, a);
                for (int i = 0; i < section.Length; i++)
                {
                    section[i] *= scale;
                }
                result[s] = section;

                // 本节直流增益 H(1) = sum(b)/sum(a)，逐级级联。
                scale *= (b[0] + b[1] + b[2]) / (a[0] + a[1] + a[2]);
            }
            return result;
        }

        /// <summary>
        /// 零相位滤波（前向 + 反向各滤一次）。等价于 <c>scipy.signal.sosfiltfilt</c>。
        /// </summary>
        /// <param name="sos">二阶节级联，每行 <c>[b0, b1, b2, a0, a1, a2]</c>，<c>a0 = 1</c>。</param>
        /// <param name="x">输入信号。</param>
        /// <param name="padtype">边界延拓方式：<c>odd</c>（默认）/ <c>even</c> / <c>constant</c>，传 null 表示不延拓。</param>
        /// <param name="padlen">延拓长度，null 表示 SciPy 的默认值（3 × ntaps）。</param>
        /// <returns>零相位滤波结果，长度与 <paramref name="x"/> 相同。</returns>
        public static double[] SosFiltFilt(double[][] sos, double[] x, string padtype = "odd", int? padlen = null)
        {
            ValidateSos(sos);
            if (x == null)
            {
                throw new ArgumentNullException(nameof(x));
            }

            int nSections = sos.Length;
            int ntaps = 2 * nSections + 1;
            int zeroB2 = 0;
            int zeroA2 = 0;
            for (int s = 0; s < nSections; s++)
            {
                if (sos[s][2] == 0.0)
                {
                    zeroB2++;
                }
                if (sos[s][5] == 0.0)
                {
                    zeroA2++;
                }
            }
            ntaps -= Math.Min(zeroB2, zeroA2);

            int edge;
            double[] extended = ValidatePad(padtype, padlen, x, ntaps, out edge);

            double[][] zi = SosFiltZi(sos);

            double[][] forwardState;
            double[] y = SosFilt(sos, extended, ScaleState(zi, extended[0]), out forwardState);

            double[][] backwardState;
            double[] reversed = Reverse(y);
            double[] y2 = SosFilt(sos, reversed, ScaleState(zi, y[y.Length - 1]), out backwardState);

            double[] result = Reverse(y2);
            if (edge > 0)
            {
                double[] trimmed = new double[result.Length - 2 * edge];
                Array.Copy(result, edge, trimmed, 0, trimmed.Length);
                result = trimmed;
            }
            return result;
        }

        /// <summary>
        /// 零相位滤波（分子/分母多项式版本）。等价于 <c>scipy.signal.filtfilt</c>（method='pad'）。
        /// </summary>
        public static double[] FiltFilt(double[] b, double[] a, double[] x, string padtype = "odd", int? padlen = null)
        {
            if (b == null || b.Length == 0)
            {
                throw new ArgumentException("分子不能为空。", nameof(b));
            }
            if (a == null || a.Length == 0)
            {
                throw new ArgumentException("分母不能为空。", nameof(a));
            }
            if (x == null)
            {
                throw new ArgumentNullException(nameof(x));
            }

            int ntaps = Math.Max(a.Length, b.Length);
            int edge;
            double[] extended = ValidatePad(padtype, padlen, x, ntaps, out edge);

            double[] zi = LFilterZi(b, a);
            double[] scaled = new double[zi.Length];
            for (int i = 0; i < zi.Length; i++)
            {
                scaled[i] = zi[i] * extended[0];
            }

            double[] zf;
            double[] y = LFilter(b, a, extended, scaled, out zf);

            double[] reversed = Reverse(y);
            double[] scaled2 = new double[zi.Length];
            for (int i = 0; i < zi.Length; i++)
            {
                scaled2[i] = zi[i] * y[y.Length - 1];
            }
            double[] y2 = LFilter(b, a, reversed, scaled2, out zf);

            double[] result = Reverse(y2);
            if (edge > 0)
            {
                double[] trimmed = new double[result.Length - 2 * edge];
                Array.Copy(result, edge, trimmed, 0, trimmed.Length);
                result = trimmed;
            }
            return result;
        }

        private static double[][] ScaleState(double[][] zi, double scale)
        {
            double[][] result = new double[zi.Length][];
            for (int s = 0; s < zi.Length; s++)
            {
                result[s] = new[] { zi[s][0] * scale, zi[s][1] * scale };
            }
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

        /// <summary>等价于 <c>scipy.signal._signaltools._validate_pad</c>。</summary>
        private static double[] ValidatePad(string padtype, int? padlen, double[] x, int ntaps, out int edge)
        {
            if (padtype != null && padtype != "even" && padtype != "odd" && padtype != "constant")
            {
                throw new ArgumentException("padtype 只能是 odd / even / constant，或 null。", nameof(padtype));
            }

            if (padtype == null)
            {
                padlen = 0;
            }

            edge = padlen.HasValue ? padlen.Value : ntaps * 3;

            if (x.Length <= edge)
            {
                throw new ArgumentException(
                    "输入长度（" + x.Length + "）必须大于 padlen（" + edge + "）。", nameof(padlen));
            }

            if (padtype == null || edge <= 0)
            {
                return x;
            }

            switch (padtype)
            {
                case "even":
                    return EvenExt(x, edge);
                case "odd":
                    return OddExt(x, edge);
                default:
                    return ConstExt(x, edge);
            }
        }

        /// <summary>奇对称延拓。等价于 <c>scipy.signal._signaltools.odd_ext</c>。</summary>
        private static double[] OddExt(double[] x, int n)
        {
            if (n < 1)
            {
                return x;
            }
            if (n > x.Length - 1)
            {
                throw new ArgumentException("延拓长度不能超过 length(x) - 1。");
            }

            double[] ext = new double[x.Length + 2 * n];
            double leftEnd = x[0];
            for (int j = 0; j < n; j++)
            {
                ext[j] = 2.0 * leftEnd - x[n - j];
            }
            Array.Copy(x, 0, ext, n, x.Length);
            double rightEnd = x[x.Length - 1];
            for (int j = 0; j < n; j++)
            {
                ext[n + x.Length + j] = 2.0 * rightEnd - x[x.Length - 2 - j];
            }
            return ext;
        }

        /// <summary>偶对称延拓。等价于 <c>scipy.signal._signaltools.even_ext</c>。</summary>
        private static double[] EvenExt(double[] x, int n)
        {
            if (n < 1)
            {
                return x;
            }
            if (n > x.Length - 1)
            {
                throw new ArgumentException("延拓长度不能超过 length(x) - 1。");
            }

            double[] ext = new double[x.Length + 2 * n];
            for (int j = 0; j < n; j++)
            {
                ext[j] = x[n - j];
            }
            Array.Copy(x, 0, ext, n, x.Length);
            for (int j = 0; j < n; j++)
            {
                ext[n + x.Length + j] = x[x.Length - 2 - j];
            }
            return ext;
        }

        /// <summary>常数延拓。等价于 <c>scipy.signal._signaltools.const_ext</c>。</summary>
        private static double[] ConstExt(double[] x, int n)
        {
            if (n < 1)
            {
                return x;
            }

            double[] ext = new double[x.Length + 2 * n];
            for (int j = 0; j < n; j++)
            {
                ext[j] = x[0];
            }
            Array.Copy(x, 0, ext, n, x.Length);
            for (int j = 0; j < n; j++)
            {
                ext[n + x.Length + j] = x[x.Length - 1];
            }
            return ext;
        }

        /// <summary>校验 SOS 结构：每行 6 个系数且 a0 = 1。</summary>
        public static void ValidateSos(double[][] sos)
        {
            if (sos == null || sos.Length == 0)
            {
                throw new ArgumentException("SOS 不能为空。", nameof(sos));
            }

            for (int i = 0; i < sos.Length; i++)
            {
                if (sos[i] == null || sos[i].Length != 6)
                {
                    throw new ArgumentException("SOS 每行必须是 6 个系数 [b0,b1,b2,a0,a1,a2]。", nameof(sos));
                }
                if (sos[i][3] != 1.0)
                {
                    throw new ArgumentException("SOS 的 a0 必须为 1。", nameof(sos));
                }
            }
        }
    }
}
