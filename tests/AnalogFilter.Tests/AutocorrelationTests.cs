using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 自相关（FAMOS <c>FFL_ACF</c> 语义）：既与 numpy 的循环自相关对照，也固化其定义性质。
    /// </summary>
    [TestClass]
    public class AutocorrelationTests
    {
        public static IEnumerable<object[]> Cases()
        {
            return DesignSpec.IdsOf("acf");
        }

        [TestMethod]
        [DynamicData(nameof(Cases))]
        public void MatchesReference(string id)
        {
            ReferenceCase reference = ReferenceData.Get(id);
            double[] input = reference.Single("X");
            double[] expected = reference.Single("ACF");

            double[] actual = Spectrum.Autocorrelation(input);

            Assert.AreEqual(expected.Length, actual.Length, id + "：输出长度");
            NumericAssert.Close(expected, actual, 1e-12, id + " ACF");
        }

        [TestMethod]
        public void Impulse_ProducesDelta()
        {
            double[] x = new double[64];
            x[0] = 1.0;

            double[] acf = Spectrum.Autocorrelation(x);

            Assert.AreEqual(1.0, acf[0], 1e-15);
            for (int i = 1; i < acf.Length; i++)
            {
                Assert.AreEqual(0.0, acf[i], 1e-15, "冲激输入的 ACF 在非零延迟处应为 0");
            }
        }

        [TestMethod]
        public void ConstantInput_IsAllOnes()
        {
            double[] x = new double[128];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = 3.5;
            }

            double[] acf = Spectrum.Autocorrelation(x);

            for (int i = 0; i < acf.Length; i++)
            {
                Assert.AreEqual(1.0, acf[i], 1e-14, "常量输入的 ACF 恒为 1（不去均值）");
            }
        }

        [TestMethod]
        public void PureSine_MatchesCosineAtIntegerLags()
        {
            const int n = 512;
            const int cycles = 3;
            double[] x = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = System.Math.Sin(2.0 * System.Math.PI * cycles * i / n);
            }

            double[] acf = Spectrum.Autocorrelation(x);

            for (int lag = 0; lag < 16; lag++)
            {
                double expected = System.Math.Cos(2.0 * System.Math.PI * cycles * lag / n);
                Assert.AreEqual(expected, acf[lag], 1e-12, "单频正弦的 ACF 应等于 cos(2π·k·lag/n)");
            }
        }

        [TestMethod]
        public void Length_IsLargestPowerOfTwoBelowInput()
        {
            double[] x = new double[1000];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = System.Math.Sin(i * 0.01);
            }

            double[] acf = Spectrum.Autocorrelation(x);

            Assert.AreEqual(512, acf.Length, "1000 点输入应截取 512 点（最大的 2 的幂）");
            Assert.AreEqual(1.0, acf[0], 1e-15, "归一化后零延迟必须为 1");
        }

        [TestMethod]
        public void LargestPowerOfTwo_Table()
        {
            Assert.AreEqual(1, Spectrum.LargestPowerOfTwo(1));
            Assert.AreEqual(2, Spectrum.LargestPowerOfTwo(3));
            Assert.AreEqual(256, Spectrum.LargestPowerOfTwo(300));
            Assert.AreEqual(1024, Spectrum.LargestPowerOfTwo(1024));
            Assert.IsTrue(Spectrum.IsPowerOfTwo(1024));
            Assert.IsFalse(Spectrum.IsPowerOfTwo(1000));
        }
    }
}
