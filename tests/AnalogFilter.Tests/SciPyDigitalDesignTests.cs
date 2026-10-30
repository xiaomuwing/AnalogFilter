using System.Collections.Generic;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 数字滤波器设计与 SciPy（<c>tools\gen_reference.py</c> 生成）逐项对照：
    /// SOS 系数、频响、滤波输出、初始状态、零相位滤波。
    /// </summary>
    [TestClass]
    public class SciPyDigitalDesignTests
    {
        private static readonly string[] PadTypes = { "odd", "even", "constant" };

        public static IEnumerable<object[]> Cases()
        {
            return DesignSpec.IdsOf("digital");
        }

        [TestMethod]
        [DynamicData(nameof(Cases))]
        public void Design_MatchesSciPy(string id)
        {
            ReferenceCase reference = ReferenceData.Get(id);
            DesignSpec spec = DesignSpec.Parse(reference.Spec);

            double[][] sos = FilterDesign.DesignSos(
                spec.Order, spec.Wn, spec.Type, spec.Band, spec.Rp, spec.Rs);

            List<double[]> expectedSos = reference.Multi("SOS");
            Assert.AreEqual(expectedSos.Count, sos.Length, id + "：SOS 节数");
            for (int section = 0; section < sos.Length; section++)
            {
                NumericAssert.Close(expectedSos[section], sos[section], 1e-11,
                    string.Format("{0} 第 {1} 节系数", id, section));
            }

            double[] w = reference.Single("W");
            NumericAssert.CloseComplex(reference.Single("FREQZ"),
                FrequencyResponse.Freqz(sos, w), 1e-10, id + " freqz");

            NumericAssert.Close(reference.Single("SOSFILT"),
                SosFiltering.SosFilt(sos, TestSignals.DesignInput()), 1e-11, id + " sosfilt");

            CheckStateInitialisation(reference, spec, sos, id);
            CheckZeroPhase(reference, spec, sos, id);
        }

        private static void CheckStateInitialisation(
            ReferenceCase reference, DesignSpec spec, double[][] sos, string id)
        {
            double[] expectedZi = reference.Single("ZI_SOS");
            if (expectedZi != null)
            {
                NumericAssert.Close(expectedZi, Flatten(SosFiltering.SosFiltZi(sos)), 1e-11,
                    id + " sosfilt_zi");
            }

            double[] expectedLfilterZi = reference.Single("ZI_LFILTER");
            if (expectedLfilterZi != null)
            {
                TransferFunction tf = FilterDesign.DesignTf(
                    spec.Order, spec.Wn, spec.Type, spec.Band, spec.Rp, spec.Rs);
                // 带阻/切比雪夫 II 的 (b,a) 形式条件数差，初始状态是差分出来的，放宽到 1e-8。
                NumericAssert.Close(expectedLfilterZi,
                    SosFiltering.LFilterZi(tf.Numerator, tf.Denominator), 1e-8, id + " lfilter_zi");
            }
        }

        /// <summary>参考数据把各节的初始状态展平成一行，比较前同样展平。</summary>
        private static double[] Flatten(double[][] rows)
        {
            int count = 0;
            foreach (double[] row in rows)
            {
                count += row.Length;
            }

            double[] result = new double[count];
            int cursor = 0;
            foreach (double[] row in rows)
            {
                System.Array.Copy(row, 0, result, cursor, row.Length);
                cursor += row.Length;
            }
            return result;
        }

        private static void CheckZeroPhase(
            ReferenceCase reference, DesignSpec spec, double[][] sos, string id)
        {
            double[] input = TestSignals.ZeroPhaseInput();
            bool hasVariants = false;

            foreach (string padtype in PadTypes)
            {
                double[] expectedDefault = reference.Single("Y_ZP_" + padtype + "_default");
                if (expectedDefault != null)
                {
                    hasVariants = true;
                    NumericAssert.Close(expectedDefault,
                        SosFiltering.SosFiltFilt(sos, input, padtype, null), 1e-9,
                        id + " sosfiltfilt padtype=" + padtype + " padlen=default");
                }

                double[] expected20 = reference.Single("Y_ZP_" + padtype + "_20");
                if (expected20 != null)
                {
                    hasVariants = true;
                    NumericAssert.Close(expected20,
                        SosFiltering.SosFiltFilt(sos, input, padtype, 20), 1e-9,
                        id + " sosfiltfilt padtype=" + padtype + " padlen=20");
                }
            }

            if (!hasVariants)
            {
                return;
            }

            TransferFunction tf = FilterDesign.DesignTf(
                spec.Order, spec.Wn, spec.Type, spec.Band, spec.Rp, spec.Rs);

            double[] expectedOdd20 = reference.Single("Y_FF_odd_20");
            if (expectedOdd20 != null)
            {
                // BA 形式（filtfilt）比 SOS 路径条件数差，且经过两次滤波，放宽到 1e-7。
                NumericAssert.Close(expectedOdd20,
                    SosFiltering.FiltFilt(tf.Numerator, tf.Denominator, input, "odd", 20), 1e-7,
                    id + " filtfilt padtype=odd padlen=20");
            }

            double[] expectedEvenDefault = reference.Single("Y_FF_even_default");
            if (expectedEvenDefault != null)
            {
                NumericAssert.Close(expectedEvenDefault,
                    SosFiltering.FiltFilt(tf.Numerator, tf.Denominator, input, "even", null), 1e-7,
                    id + " filtfilt padtype=even padlen=default");
            }
        }
    }
}
