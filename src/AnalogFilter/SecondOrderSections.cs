using System;
using System.Collections.Generic;
using System.Numerics;

namespace AnalogFilter
{
    /// <summary>
    /// 零极点增益 → 二阶节（SOS）级联，以及配套的共轭配对工具。
    /// 移植自 SciPy 1.18.1 的 <c>zpk2sos</c> / <c>_cplxreal</c> /
    /// <c>_nearest_real_complex_idx</c> / <c>_single_zpksos</c>，BSD-3-Clause。
    /// </summary>
    /// <remarks>
    /// SOS 每行 6 个系数：<c>[b0, b1, b2, a0, a1, a2]</c>，且 <c>a0 = 1</c>，
    /// 与 <c>scipy.signal.sosfilt</c> 的约定一致。
    /// 数字滤波器使用 "nearest" 配对，模拟滤波器使用 "minimal" 配对。
    /// </remarks>
    public static class SecondOrderSections
    {
        /// <summary>把 ZPK 转换成 SOS 级联。等价于 <c>scipy.signal.zpk2sos</c>（pairing 默认值）。</summary>
        /// <param name="zpk">要转换的零极点增益模型。</param>
        /// <param name="analog">true 表示模拟滤波器（minimal 配对）；false 表示数字滤波器（nearest 配对）。</param>
        /// <returns>SOS 级联，每行 <c>[b0, b1, b2, a0, a1, a2]</c>，且 <c>a0 = 1</c>。</returns>
        public static double[][] Zpk2Sos(ZpkModel zpk, bool analog = false)
        {
            if (zpk == null)
            {
                throw new ArgumentNullException(nameof(zpk));
            }
            return Zpk2Sos(zpk.Zeros, zpk.Poles, zpk.Gain, analog, analog ? "minimal" : "nearest");
        }

        /// <summary>把 ZPK 转换成 SOS 级联。等价于 <c>scipy.signal.zpk2sos</c>。</summary>
        /// <param name="zeros">零点集合，允许为空。</param>
        /// <param name="poles">极点集合，允许为空。</param>
        /// <param name="gain">增益系数 k。</param>
        /// <param name="analog">true 表示模拟滤波器；false 表示数字滤波器。</param>
        /// <param name="pairing">配对方式：<c>nearest</c> / <c>keep_odd</c> / <c>minimal</c>，null 表示按 <paramref name="analog"/> 取默认值。</param>
        /// <returns>SOS 级联，每行 <c>[b0, b1, b2, a0, a1, a2]</c>。</returns>
        public static double[][] Zpk2Sos(
            IReadOnlyList<Complex> zeros,
            IReadOnlyList<Complex> poles,
            double gain,
            bool analog = false,
            string pairing = null)
        {
            if (pairing == null)
            {
                pairing = analog ? "minimal" : "nearest";
            }
            if (pairing != "nearest" && pairing != "keep_odd" && pairing != "minimal")
            {
                throw new ArgumentException("pairing 只能是 nearest / keep_odd / minimal。", nameof(pairing));
            }
            if (analog && pairing != "minimal")
            {
                throw new ArgumentException("模拟滤波器的 zpk2sos 只支持 minimal 配对。", nameof(pairing));
            }

            List<Complex> zList = new List<Complex>(zeros ?? Array.Empty<Complex>());
            List<Complex> pList = new List<Complex>(poles ?? Array.Empty<Complex>());

            if (zList.Count == 0 && pList.Count == 0)
            {
                // 与 scipy 的退化分支一致。
                return analog
                    ? new[] { new[] { 0.0, 0.0, gain, 0.0, 0.0, 1.0 } }
                    : new[] { new[] { gain, 0.0, 0.0, 1.0, 0.0, 0.0 } };
            }

            int sectionCount;
            if (pairing != "minimal")
            {
                // 先补齐极点（用原长度判断），再按补齐后的极点数补齐零点。
                int padPoles = Math.Max(zList.Count - pList.Count, 0);
                for (int i = 0; i < padPoles; i++)
                {
                    pList.Add(Complex.Zero);
                }

                int padZeros = Math.Max(pList.Count - zList.Count, 0);
                for (int i = 0; i < padZeros; i++)
                {
                    zList.Add(Complex.Zero);
                }

                sectionCount = (Math.Max(pList.Count, zList.Count) + 1) / 2;

                if (pList.Count % 2 == 1 && pairing == "nearest")
                {
                    pList.Add(Complex.Zero);
                    zList.Add(Complex.Zero);
                }
            }
            else
            {
                if (pList.Count < zList.Count)
                {
                    throw new ArgumentException("模拟 zpk2sos 要求极点数不少于零点数。", nameof(poles));
                }
                sectionCount = (pList.Count + 1) / 2;
            }

            // 展开成"每对共轭只保留一个代表 + 实数"的列表。
            zList = ConcatCplxReal(zList);
            pList = ConcatCplxReal(pList);

            double[][] sos = new double[sectionCount][];
            for (int i = 0; i < sectionCount; i++)
            {
                sos[i] = new double[6];
            }

            // 从最后一节往前填，使"最差"的极点排在最后。
            for (int si = sectionCount - 1; si >= 0; si--)
            {
                int p1Index = WorstPoleIndex(pList, analog);
                Complex p1 = pList[p1Index];
                pList.RemoveAt(p1Index);

                if (IsReal(p1) && CountReal(pList) == 0)
                {
                    // 特例 1：只剩一个实极点。
                    if (pairing != "minimal")
                    {
                        int z1Index = NearestRealComplexIndex(zList, p1, "real");
                        Complex z1 = zList[z1Index];
                        zList.RemoveAt(z1Index);
                        sos[si] = SingleZpkSos(new[] { z1, Complex.Zero }, new[] { p1, Complex.Zero }, 1.0);
                    }
                    else if (zList.Count > 0)
                    {
                        int z1Index = NearestRealComplexIndex(zList, p1, "real");
                        Complex z1 = zList[z1Index];
                        zList.RemoveAt(z1Index);
                        sos[si] = SingleZpkSos(new[] { z1 }, new[] { p1 }, 1.0);
                    }
                    else
                    {
                        sos[si] = SingleZpkSos(Array.Empty<Complex>(), new[] { p1 }, 1.0);
                    }
                }
                else if (pList.Count + 1 == zList.Count
                         && !IsReal(p1)
                         && CountReal(pList) == 1
                         && CountReal(zList) == 1)
                {
                    // 特例 2：各剩一个实极点/实零点，必须与复零点配对。
                    int z1Index = NearestRealComplexIndex(zList, p1, "complex");
                    Complex z1 = zList[z1Index];
                    zList.RemoveAt(z1Index);
                    sos[si] = SingleZpkSos(
                        new[] { z1, Complex.Conjugate(z1) },
                        new[] { p1, Complex.Conjugate(p1) },
                        1.0);
                }
                else
                {
                    Complex p2;
                    if (IsReal(p1))
                    {
                        List<int> realIndices = new List<int>();
                        for (int i = 0; i < pList.Count; i++)
                        {
                            if (IsReal(pList[i]))
                            {
                                realIndices.Add(i);
                            }
                        }

                        List<Complex> realPoles = new List<Complex>();
                        foreach (int index in realIndices)
                        {
                            realPoles.Add(pList[index]);
                        }

                        int worstAmongReal = WorstPoleIndex(realPoles, analog);
                        int p2Index = realIndices[worstAmongReal];
                        p2 = pList[p2Index];
                        pList.RemoveAt(p2Index);
                    }
                    else
                    {
                        p2 = Complex.Conjugate(p1);
                    }

                    if (zList.Count > 0)
                    {
                        int z1Index = NearestRealComplexIndex(zList, p1, "any");
                        Complex z1 = zList[z1Index];
                        zList.RemoveAt(z1Index);

                        if (!IsReal(z1))
                        {
                            sos[si] = SingleZpkSos(
                                new[] { z1, Complex.Conjugate(z1) },
                                new[] { p1, p2 },
                                1.0);
                        }
                        else if (zList.Count > 0)
                        {
                            int z2Index = NearestRealComplexIndex(zList, p1, "real");
                            Complex z2 = zList[z2Index];
                            zList.RemoveAt(z2Index);
                            sos[si] = SingleZpkSos(new[] { z1, z2 }, new[] { p1, p2 }, 1.0);
                        }
                        else
                        {
                            sos[si] = SingleZpkSos(new[] { z1 }, new[] { p1, p2 }, 1.0);
                        }
                    }
                    else
                    {
                        sos[si] = SingleZpkSos(Array.Empty<Complex>(), new[] { p1, p2 }, 1.0);
                    }
                }
            }

            // 增益统一放到第一节。
            sos[0][0] *= gain;
            sos[0][1] *= gain;
            sos[0][2] *= gain;

            return sos;
        }

        /// <summary>单个二阶节的 ZPK → SOS 行。</summary>
        private static double[] SingleZpkSos(IReadOnlyList<Complex> zeros, IReadOnlyList<Complex> poles, double gain)
        {
            double[] b;
            double[] a;
            TransferFunctionConversion.Zpk2Tf(zeros, poles, gain, out b, out a);

            double[] sos = new double[6];
            for (int i = 0; i < b.Length && i < 3; i++)
            {
                sos[3 - b.Length + i] = b[i];
            }
            for (int i = 0; i < a.Length && i < 3; i++)
            {
                sos[6 - a.Length + i] = a[i];
            }
            return sos;
        }

        /// <summary>数字域取"最靠近单位圆"的极点，模拟域取"最靠近虚轴"的极点。</summary>
        private static int WorstPoleIndex(IReadOnlyList<Complex> poles, bool analog)
        {
            int best = 0;
            double bestValue = double.MaxValue;
            for (int i = 0; i < poles.Count; i++)
            {
                double value = analog
                    ? Math.Abs(poles[i].Real)
                    : Math.Abs(1.0 - Complex.Abs(poles[i]));
                if (value < bestValue)
                {
                    bestValue = value;
                    best = i;
                }
            }
            return best;
        }

        private static int NearestRealComplexIndex(IReadOnlyList<Complex> candidates, Complex target, string which)
        {
            int[] order = new int[candidates.Count];
            double[] distance = new double[candidates.Count];
            for (int i = 0; i < candidates.Count; i++)
            {
                order[i] = i;
                distance[i] = Complex.Abs(candidates[i] - target);
            }
            Array.Sort(distance, order);

            if (which == "any")
            {
                return order[0];
            }

            bool wantReal = which == "real";
            for (int i = 0; i < order.Length; i++)
            {
                bool isReal = IsReal(candidates[order[i]]);
                if (isReal == wantReal)
                {
                    return order[i];
                }
            }

            throw new InvalidOperationException("找不到所需的实/复零点，零极点配对失败。");
        }

        private static bool IsReal(Complex value) => value.Imaginary == 0.0;

        private static int CountReal(IReadOnlyList<Complex> values)
        {
            int count = 0;
            for (int i = 0; i < values.Count; i++)
            {
                if (IsReal(values[i]))
                {
                    count++;
                }
            }
            return count;
        }

        /// <summary>把 _cplxreal 的两个返回值拼接起来（与 scipy 的 np.concatenate(_cplxreal(x)) 等价）。</summary>
        private static List<Complex> ConcatCplxReal(List<Complex> values)
        {
            List<Complex> complexPart;
            List<double> realPart;
            CplxReal(values, out complexPart, out realPart);

            List<Complex> result = new List<Complex>(complexPart);
            foreach (double r in realPart)
            {
                result.Add(new Complex(r, 0.0));
            }
            return result;
        }

        /// <summary>
        /// 分成"复数共轭对代表（正虚部）"与"实数"两部分。等价于 <c>scipy.signal._cplxreal</c>。
        /// </summary>
        private static void CplxReal(List<Complex> values, out List<Complex> complexPart, out List<double> realPart)
        {
            complexPart = new List<Complex>();
            realPart = new List<double>();
            if (values == null || values.Count == 0)
            {
                return;
            }

            double tolerance = 100.0 * 2.2204460492503131e-16;

            List<Complex> sorted = new List<Complex>(values);
            sorted.Sort((a, b) =>
            {
                int cmp = a.Real.CompareTo(b.Real);
                return cmp != 0 ? cmp : Math.Abs(a.Imaginary).CompareTo(Math.Abs(b.Imaginary));
            });

            List<Complex> nonReal = new List<Complex>();
            foreach (Complex value in sorted)
            {
                if (Math.Abs(value.Imaginary) <= tolerance * Complex.Abs(value))
                {
                    realPart.Add(value.Real);
                }
                else
                {
                    nonReal.Add(value);
                }
            }

            if (nonReal.Count == 0)
            {
                return;
            }

            List<Complex> positive = new List<Complex>();
            List<Complex> negative = new List<Complex>();
            foreach (Complex value in nonReal)
            {
                if (value.Imaginary > 0)
                {
                    positive.Add(value);
                }
                else
                {
                    negative.Add(value);
                }
            }

            if (positive.Count != negative.Count)
            {
                throw new ArgumentException("存在没有共轭配对的复零点/极点。");
            }

            // 按实部把相邻的复根分成若干"同一实部"的段，段内按 |虚部| 排序。
            SortRunsByImaginary(positive, negative, tolerance);

            for (int i = 0; i < positive.Count; i++)
            {
                if (Complex.Abs(positive[i] - Complex.Conjugate(negative[i])) > tolerance * Complex.Abs(negative[i]))
                {
                    throw new ArgumentException("存在没有共轭配对的复零点/极点。");
                }
            }

            for (int i = 0; i < positive.Count; i++)
            {
                complexPart.Add((positive[i] + Complex.Conjugate(negative[i])) / 2.0);
            }
        }

        private static void SortRunsByImaginary(List<Complex> positive, List<Complex> negative, double tolerance)
        {
            int n = positive.Count;
            if (n < 2)
            {
                return;
            }

            int start = 0;
            while (start < n)
            {
                int stop = start + 1;
                while (stop < n &&
                       Math.Abs(positive[stop].Real - positive[stop - 1].Real) <= tolerance * Complex.Abs(positive[stop - 1]))
                {
                    stop++;
                }

                if (stop - start > 1)
                {
                    SortRangeByImaginary(positive, start, stop);
                    SortRangeByImaginary(negative, start, stop);
                }

                start = stop;
            }
        }

        private static void SortRangeByImaginary(List<Complex> values, int start, int stop)
        {
            List<Complex> slice = values.GetRange(start, stop - start);
            slice.Sort((a, b) => Math.Abs(a.Imaginary).CompareTo(Math.Abs(b.Imaginary)));
            for (int i = 0; i < slice.Count; i++)
            {
                values[start + i] = slice[i];
            }
        }
    }
}
