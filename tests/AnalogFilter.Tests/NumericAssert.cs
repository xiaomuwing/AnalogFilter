using System;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 统一的数值比较口径：相对误差，并对接近 0 的元素改用数组峰值的 1e-4 作为下限，
    /// 避免分母趋零时把舍入残差判成失败。
    /// </summary>
    internal static class NumericAssert
    {
        internal static void Close(double[] expected, double[] actual, double relativeTolerance, string context)
        {
            if (expected == null || actual == null)
            {
                Assert.Fail(context + "：缺少数据（参考或实际为 null）");
            }
            if (expected.Length != actual.Length)
            {
                Assert.Fail(string.Format(
                    "{0}：长度 {1} ≠ {2}", context, actual.Length, expected.Length));
            }

            double scale = 0.0;
            for (int i = 0; i < expected.Length; i++)
            {
                scale = Math.Max(scale, Math.Abs(expected[i]));
            }
            double floor = Math.Max(scale * 1e-4, 1e-300);

            double worst = 0.0;
            int worstIndex = -1;
            for (int i = 0; i < expected.Length; i++)
            {
                double a = expected[i];
                double b = actual[i];
                double denominator = Math.Max(Math.Max(Math.Abs(a), Math.Abs(b)), floor);
                double error = Math.Abs(a - b) / denominator;
                if (error > worst)
                {
                    worst = error;
                    worstIndex = i;
                }
            }

            if (worst > relativeTolerance)
            {
                Assert.Fail(string.Format(
                    "{0}：最大相对误差 {1:E3} 超过上限 {2:E1}（索引 {3}：参考 {4:G17}，实际 {5:G17}）",
                    context, worst, relativeTolerance, worstIndex,
                    expected[worstIndex], actual[worstIndex]));
            }
        }

        /// <summary>比较复序列，参考数据按 [re, im, re, im, ...] 交错存放。</summary>
        internal static void CloseComplex(double[] interleaved, Complex[] actual, double relativeTolerance, string context)
        {
            if (interleaved.Length != actual.Length * 2)
            {
                Assert.Fail(string.Format(
                    "{0}：复序列长度 {1} ≠ 参考 {2}", context, actual.Length, interleaved.Length / 2));
            }

            double[] expected = new double[interleaved.Length];
            Array.Copy(interleaved, expected, interleaved.Length);
            double[] values = new double[actual.Length * 2];
            for (int i = 0; i < actual.Length; i++)
            {
                values[2 * i] = actual[i].Real;
                values[2 * i + 1] = actual[i].Imaginary;
            }
            Close(expected, values, relativeTolerance, context);
        }

        internal static void Scalar(double expected, double actual, double relativeTolerance, string context)
        {
            Close(new[] { expected }, new[] { actual }, relativeTolerance, context);
        }
    }
}
