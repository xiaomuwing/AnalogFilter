using System;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 不依赖 SciPy 参照值的定义性测试：−3 dB 点、直流增益、稳定性、节数、零相位特性等。
    /// 这些断言直接来自滤波器的数学定义，因此能独立于参考数据发现回归。
    /// </summary>
    [TestClass]
    public class DesignInvariantTests
    {
        /// <summary>20·log10(1/√2)。</summary>
        private const double Minus3dB = -3.0102999566398116;

        [TestMethod]
        public void Butterworth_HasMinus3dBAtCutoff_AndUnityDcGain()
        {
            foreach (int order in new[] { 2, 4, 8 })
            {
                double[][] sos = FilterDesign.DesignSos(order, new[] { 0.2 });

                double atCutoff = FrequencyResponse.MagnitudeDb(
                    FrequencyResponse.Freqz(sos, Nyquist(0.2)))[0];
                double atDc = FrequencyResponse.MagnitudeDb(
                    FrequencyResponse.Freqz(sos, new[] { 0.0 }))[0];

                Assert.AreEqual(Minus3dB, atCutoff, 1e-9,
                    "巴特沃斯的 −3 dB 点必须落在截止频率上（N=" + order + "）");
                Assert.AreEqual(0.0, atDc, 1e-12,
                    "巴特沃斯低通的直流增益必须为 1（N=" + order + "）");
            }
        }

        [TestMethod]
        public void ChebyshevI3dB_UsesFamosConvention()
        {
            // FAMOS 约定：−3 dB 在截止频率，且直流增益为 1（区别于 SciPy 的波纹边约定）。
            foreach (int order in new[] { 2, 3, 4, 5 })
            {
                double[][] sos = FilterDesign.DesignSos(
                    order, new[] { 0.2 }, FilterType.ChebyshevI3dB, BandType.LowPass, 1.0);

                double atCutoff = FrequencyResponse.MagnitudeDb(
                    FrequencyResponse.Freqz(sos, Nyquist(0.2)))[0];
                double atDc = FrequencyResponse.MagnitudeDb(
                    FrequencyResponse.Freqz(sos, new[] { 0.0 }))[0];

                Assert.AreEqual(Minus3dB, atCutoff, 1e-9,
                    "ChebyshevI3dB 的 −3 dB 点必须落在截止频率上（N=" + order + "）");
                Assert.AreEqual(0.0, atDc, 1e-12,
                    "ChebyshevI3dB 的直流增益必须为 1（N=" + order + "）");
            }
        }

        [TestMethod]
        public void CriticallyDamped_Minus3dBAtBandEdges_AllBands()
        {
            double[][] lowPass = FilterDesign.DesignSos(
                4, new[] { 0.2 }, FilterType.CriticallyDamped, BandType.LowPass);
            Assert.AreEqual(Minus3dB,
                FrequencyResponse.MagnitudeDb(FrequencyResponse.Freqz(lowPass, Nyquist(0.2)))[0],
                1e-9, "临界阻尼低通的 −3 dB 点");

            double[][] highPass = FilterDesign.DesignSos(
                4, new[] { 0.2 }, FilterType.CriticallyDamped, BandType.HighPass);
            Assert.AreEqual(Minus3dB,
                FrequencyResponse.MagnitudeDb(FrequencyResponse.Freqz(highPass, Nyquist(0.2)))[0],
                1e-9, "临界阻尼高通的 −3 dB 点");

            double[][] bandPass = FilterDesign.DesignSos(
                3, new[] { 0.1, 0.3 }, FilterType.CriticallyDamped, BandType.BandPass);
            double[] bpEdges = FrequencyResponse.MagnitudeDb(
                FrequencyResponse.Freqz(bandPass, Nyquist(0.1, 0.3)));
            Assert.AreEqual(Minus3dB, bpEdges[0], 1e-9, "临界阻尼带通的下带边");
            Assert.AreEqual(Minus3dB, bpEdges[1], 1e-9, "临界阻尼带通的上带边");

            double[][] bandStop = FilterDesign.DesignSos(
                3, new[] { 0.1, 0.3 }, FilterType.CriticallyDamped, BandType.BandStop);
            double[] bsEdges = FrequencyResponse.MagnitudeDb(
                FrequencyResponse.Freqz(bandStop, Nyquist(0.1, 0.3)));
            Assert.AreEqual(Minus3dB, bsEdges[0], 1e-9, "临界阻尼带阻的下带边");
            Assert.AreEqual(Minus3dB, bsEdges[1], 1e-9, "临界阻尼带阻的上带边");
        }

        [TestMethod]
        public void CriticallyDamped_StepResponseIsMonotonic()
        {
            double[][] sos = FilterDesign.DesignSos(
                4, new[] { 0.1 }, FilterType.CriticallyDamped, BandType.LowPass);
            double[] y = SosFiltering.SosFilt(sos, TestSignals.Step(512));

            for (int i = 1; i < y.Length; i++)
            {
                Assert.IsTrue(y[i] >= y[i - 1] - 1e-12,
                    "临界阻尼的阶跃响应不得过冲（索引 " + i + "）");
            }
            Assert.AreEqual(1.0, y[y.Length - 1], 1e-3, "阶跃响应最终应收敛到 1");
        }

        [TestMethod]
        public void Bessel_MagnitudeNormalisation_HasMinus3dBAtCutoff()
        {
            double[][] sos = FilterDesign.DesignSos(
                4, new[] { 0.2 }, FilterType.Bessel, BandType.LowPass,
                1.0, 20.0, false, null, BesselNormalization.Magnitude);

            double atCutoff = FrequencyResponse.MagnitudeDb(
                FrequencyResponse.Freqz(sos, Nyquist(0.2)))[0];

            Assert.AreEqual(Minus3dB, atCutoff, 1e-9,
                "BesselNormalization.Magnitude 的 −3 dB 点必须落在截止频率上");
        }

        [TestMethod]
        public void SosSectionCount_MatchesPoleCount()
        {
            foreach (BandType band in new[] { BandType.LowPass, BandType.HighPass, BandType.BandPass, BandType.BandStop })
            {
                double[] wn = band == BandType.BandPass || band == BandType.BandStop
                    ? new[] { 0.1, 0.3 }
                    : new[] { 0.2 };

                for (int order = 2; order <= 5; order++)
                {
                    ZpkModel zpk = FilterDesign.DesignZpk(order, wn, FilterType.Butterworth, band);
                    double[][] sos = SecondOrderSections.Zpk2Sos(zpk);

                    int expectedSections = (zpk.Poles.Length + 1) / 2;
                    Assert.AreEqual(expectedSections, sos.Length,
                        string.Format("{0} 阶 {1} 的节数", order, band));
                    for (int section = 0; section < sos.Length; section++)
                    {
                        Assert.AreEqual(1.0, sos[section][3], 1e-15, "SOS 的 a0 必须为 1");
                    }
                }
            }
        }

        [TestMethod]
        public void DigitalFilters_AreStable()
        {
            foreach (FilterType type in new[]
                     { FilterType.Butterworth, FilterType.ChebyshevI, FilterType.ChebyshevII, FilterType.Bessel })
            {
                ZpkModel zpk = FilterDesign.DesignZpk(6, new[] { 0.2 }, type);
                foreach (Complex pole in zpk.Poles)
                {
                    Assert.IsTrue(pole.Magnitude < 1.0,
                        string.Format("{0} 的极点必须严格位于单位圆内（|p|={1:G5}）", type, pole.Magnitude));
                }
            }
        }

        [TestMethod]
        public void SosFilt_AgreesWithLFilter()
        {
            double[][] sos = FilterDesign.DesignSos(4, new[] { 0.2 });
            TransferFunction tf = FilterDesign.DesignTf(4, new[] { 0.2 });

            double[] x = TestSignals.DesignInput();
            double[] bySos = SosFiltering.SosFilt(sos, x);
            double[] byBa = SosFiltering.LFilter(tf.Numerator, tf.Denominator, x);

            NumericAssert.Close(byBa, bySos, 1e-10, "SOS 级联与直接 II 型转置应给出同一个滤波器");
        }

        [TestMethod]
        public void SosFiltFilt_IsZeroPhase()
        {
            // TestSignals.Sine 的参数是「每周期的采样数」倒数，即 cycles/sample；
            // 而 Freqz 的 w 是 rad/sample，两者相差 2π。
            const double cyclesPerSample = 0.05;
            const int length = 400;

            double[][] sos = FilterDesign.DesignSos(4, new[] { 0.2 });
            double[] x = TestSignals.Sine(cyclesPerSample, length);
            double[] y = SosFiltering.SosFiltFilt(sos, x);

            // 零相位滤波对纯正弦相当于乘以 |H|²，且不产生相移。
            // 只比较中段，避开延拓边界带来的瞬态。
            Complex h = FrequencyResponse.Freqz(sos, new[] { 2.0 * Math.PI * cyclesPerSample })[0];
            double gain = h.Magnitude * h.Magnitude;

            double worst = 0.0;
            int worstIndex = -1;
            for (int i = length / 4; i < length - length / 4; i++)
            {
                double residual = Math.Abs(gain * x[i] - y[i]);
                if (residual > worst)
                {
                    worst = residual;
                    worstIndex = i;
                }
            }

            Assert.IsTrue(worst <= 1e-6,
                string.Format("零相位滤波后的正弦应与输入同相：最大残留 {0:E2}（索引 {1}）", worst, worstIndex));
        }

        [TestMethod]
        public void SampleRate_ParameterIsEquivalentToNormalisedFrequency()
        {
            double[][] byRate = FilterDesign.DesignSos(
                4, new[] { 50.0 }, FilterType.Butterworth, BandType.LowPass,
                1.0, 20.0, false, 1000.0);
            double[][] byNormalised = FilterDesign.DesignSos(4, new[] { 0.1 });

            Assert.AreEqual(byNormalised.Length, byRate.Length, "两种写法应得到相同的节数");
            for (int section = 0; section < byRate.Length; section++)
            {
                NumericAssert.Close(byNormalised[section], byRate[section], 1e-15,
                    "sampleRate=1000、wn=50 Hz 应等价于 wn=0.1");
            }
        }

        [TestMethod]
        public void Freqz_MatchesDirectPolynomialEvaluation()
        {
            double[] w = { 0.01, 0.1, 0.25, 0.5, 0.9 };
            double[][] sos = FilterDesign.DesignSos(4, new[] { 0.2 });

            Complex[] byImplementation = FrequencyResponse.Freqz(sos, w);

            for (int i = 0; i < w.Length; i++)
            {
                Complex z = Complex.FromPolarCoordinates(1.0, -w[i]);
                Complex response = Complex.One;
                foreach (double[] section in sos)
                {
                    Complex numerator = section[0] + section[1] * z + section[2] * z * z;
                    Complex denominator = section[3] + section[4] * z + section[5] * z * z;
                    response *= numerator / denominator;
                }

                Assert.AreEqual(0.0, (response - byImplementation[i]).Magnitude, 1e-12,
                    "SOS 频响应等于逐节多项式的乘积（w=" + w[i] + "）");
            }
        }

        [TestMethod]
        public void Tf2Zpk_RoundTripsButterworth()
        {
            ZpkModel original = FilterDesign.DesignZpk(4, new[] { 0.2 });
            TransferFunction tf = TransferFunctionConversion.Zpk2Tf(original);
            ZpkModel roundTripped = TransferFunctionConversion.Tf2Zpk(tf.Numerator, tf.Denominator);

            Assert.AreEqual(original.Poles.Length, roundTripped.Poles.Length);
            ComplexArrays.AssertSameSet(
                ToInterleaved(original.Poles), roundTripped.Poles, 1e-8, "ZPK→TF→ZPK 的极点");
            NumericAssert.Scalar(original.Gain, roundTripped.Gain, 1e-10, "ZPK→TF→ZPK 的增益");
        }

        [TestMethod]
        public void AnalogButterworthPrototype_HasExpectedShape()
        {
            ZpkModel zpk = AnalogPrototypes.Butterworth(4);

            Assert.AreEqual(4, zpk.Poles.Length, "4 阶巴特沃斯原型应有 4 个极点");
            Assert.AreEqual(0, zpk.Zeros.Length, "巴特沃斯原型没有零点");
            Assert.AreEqual(1.0, zpk.Gain, 1e-15, "原型的直流增益为 1");
            foreach (Complex pole in zpk.Poles)
            {
                Assert.AreEqual(1.0, pole.Magnitude, 1e-15, "原型的极点都落在单位圆上");
            }
        }

        /// <summary>
        /// 设计 API 的 <c>wn</c> 归一化到 Nyquist（0..1），而 <see cref="FrequencyResponse.Freqz"/>
        /// 的 <c>w</c> 单位是 rad/sample（0..π），两者差一个 π。
        /// </summary>
        private static double[] Nyquist(params double[] normalizedFrequencies)
        {
            double[] result = new double[normalizedFrequencies.Length];
            for (int i = 0; i < normalizedFrequencies.Length; i++)
            {
                result[i] = normalizedFrequencies[i] * Math.PI;
            }
            return result;
        }

        private static double[] ToInterleaved(Complex[] values)
        {
            double[] result = new double[values.Length * 2];
            for (int i = 0; i < values.Length; i++)
            {
                result[2 * i] = values[i].Real;
                result[2 * i + 1] = values[i].Imaginary;
            }
            return result;
        }
    }
}
