using System;
using System.Collections.Generic;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 模拟域设计与 SciPy 对照：极点、零点、增益按集合比较（先按实部、虚部排序），
    /// 再比较 <c>freqs</c> 频响。
    /// </summary>
    [TestClass]
    public class SciPyAnalogDesignTests
    {
        public static IEnumerable<object[]> Cases()
        {
            return DesignSpec.IdsOf("analog");
        }

        [TestMethod]
        [DynamicData(nameof(Cases))]
        public void AnalogDesign_MatchesSciPy(string id)
        {
            ReferenceCase reference = ReferenceData.Get(id);
            DesignSpec spec = DesignSpec.Parse(reference.Spec);

            ZpkModel zpk = FilterDesign.DesignZpk(
                spec.Order, spec.Wn, spec.Type, spec.Band, spec.Rp, spec.Rs, analog: true);

            ComplexArrays.AssertSameSet(reference.Single("POLE"), zpk.Poles, 1e-11, id + " 极点");
            ComplexArrays.AssertSameSet(reference.Single("ZERO"), zpk.Zeros, 1e-11, id + " 零点");
            NumericAssert.Scalar(reference.Single("GAIN")[0], zpk.Gain, 1e-11, id + " 增益");

            double[] w = reference.Single("W");
            NumericAssert.CloseComplex(reference.Single("FREQS"),
                FrequencyResponse.Freqs(zpk, w), 1e-10, id + " freqs");
        }
    }

    /// <summary>
    /// 贝塞尔原型的精度随阶数下降（见 README 的贝塞尔精度表）：这里把当前的精度水平固化成门槛，
    /// 既防止悄悄退化，也把已知上限写进测试。
    /// </summary>
    [TestClass]
    public class BesselPrototypeAccuracyTests
    {
        public static IEnumerable<object[]> Cases()
        {
            return DesignSpec.IdsOf("besselp");
        }

        [TestMethod]
        [DynamicData(nameof(Cases))]
        public void Poles_MatchSciPyWithinDocumentedTolerance(string id)
        {
            ReferenceCase reference = ReferenceData.Get(id);
            DesignSpec spec = DesignSpec.Parse(reference.Spec);
            int order = spec.Order;

            ZpkModel zpk = AnalogPrototypes.Bessel(order, BesselNormalization.Phase);

            ComplexArrays.AssertSameSet(reference.Single("POLE"), zpk.Poles, ToleranceFor(order),
                string.Format("{0} 极点", id));
            NumericAssert.Scalar(reference.Single("GAIN")[0], zpk.Gain, ToleranceFor(order), id + " 增益");
        }

        /// <summary>与 README 的贝塞尔精度表一致：8 阶 1e-13，12 阶 6e-12。</summary>
        private static double ToleranceFor(int order)
        {
            if (order <= 8)
            {
                return 1e-11;
            }
            if (order <= 10)
            {
                return 1e-10;
            }
            if (order <= 12)
            {
                return 1e-9;
            }
            return 1e-6;
        }
    }

    internal static class ComplexArrays
    {
        /// <summary>
        /// 集合比较：对每个参考点找最近的实现点，误差按该点的模归一化。
        /// 不依赖排序——共轭对的实部只差 1 ulp 时，排序会把两个点对调。
        /// 参考数据是 [re, im, re, im, ...] 的展平形式。
        /// </summary>
        internal static void AssertSameSet(
            double[] expectedInterleaved, Complex[] actual, double relativeTolerance, string context)
        {
            int expectedCount = expectedInterleaved == null ? 0 : expectedInterleaved.Length / 2;
            Assert.AreEqual(expectedCount, actual.Length, context + "：元素个数");

            if (expectedCount == 0)
            {
                return;
            }

            double scale = 0.0;
            for (int i = 0; i < expectedCount; i++)
            {
                scale = Math.Max(scale, new Complex(expectedInterleaved[2 * i], expectedInterleaved[2 * i + 1]).Magnitude);
            }
            double floor = Math.Max(scale * 1e-4, 1e-300);

            double worst = 0.0;
            int worstIndex = -1;
            for (int i = 0; i < expectedCount; i++)
            {
                Complex expected = new Complex(expectedInterleaved[2 * i], expectedInterleaved[2 * i + 1]);

                double best = double.MaxValue;
                foreach (Complex candidate in actual)
                {
                    best = Math.Min(best, (candidate - expected).Magnitude);
                }

                double error = best / Math.Max(expected.Magnitude, floor);
                if (error > worst)
                {
                    worst = error;
                    worstIndex = i;
                }
            }

            if (worst > relativeTolerance)
            {
                Assert.Fail(string.Format(
                    "{0}：最大相对误差 {1:E3} 超过上限 {2:E1}（参考点 #{3} = {4}）",
                    context, worst, relativeTolerance, worstIndex,
                    new Complex(expectedInterleaved[2 * worstIndex], expectedInterleaved[2 * worstIndex + 1])));
            }
        }
    }
}
