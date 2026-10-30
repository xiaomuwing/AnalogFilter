using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace AnalogFilter.Tests
{
    /// <summary>参数校验：非法输入必须抛出明确的异常，而不是静默给出错误结果。</summary>
    [TestClass]
    public class ArgumentValidationTests
    {
        [TestMethod]
        public void NegativeOrder_IsRejected()
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => FilterDesign.DesignSos(-1, new[] { 0.2 }));
        }

        [TestMethod]
        public void MissingFrequencies_AreRejected()
        {
            Assert.ThrowsExactly<ArgumentNullException>(
                () => FilterDesign.DesignSos(2, null));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => FilterDesign.DesignSos(2, new[] { 0.0 }));
            Assert.ThrowsExactly<ArgumentException>(
                () => FilterDesign.DesignSos(2, new[] { 0.1, 0.3 }));
        }

        [TestMethod]
        public void FrequencyAtOrAboveNyquist_IsRejectedForDigitalFilters()
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => FilterDesign.DesignSos(2, new[] { 1.0 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => FilterDesign.DesignSos(2, new[] { 600.0 }, sampleRate: 1000.0));
        }

        [TestMethod]
        public void ReversedBandEdges_AreRejected()
        {
            Assert.ThrowsExactly<ArgumentException>(
                () => FilterDesign.DesignSos(2, new[] { 0.3, 0.1 }, FilterType.Butterworth, BandType.BandPass));
        }

        [TestMethod]
        public void AnalogWithSampleRate_IsRejected()
        {
            Assert.ThrowsExactly<ArgumentException>(
                () => FilterDesign.DesignZpk(2, new[] { 100.0 }, FilterType.Butterworth,
                    BandType.LowPass, 1.0, 20.0, true, 1000.0));
        }

        [TestMethod]
        public void InvalidSampleRate_IsRejected()
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => FilterDesign.DesignSos(2, new[] { 50.0 }, FilterType.Butterworth,
                    BandType.LowPass, 1.0, 20.0, false, 0.0));
        }

        [TestMethod]
        public void UnknownPadType_IsRejected()
        {
            double[][] sos = FilterDesign.DesignSos(2, new[] { 0.2 });
            Assert.ThrowsExactly<ArgumentException>(
                () => SosFiltering.SosFiltFilt(sos, TestSignals.DesignInput(), "mirror"));
        }

        [TestMethod]
        public void InputShorterThanPad_IsRejected()
        {
            double[][] sos = FilterDesign.DesignSos(4, new[] { 0.2 });
            Assert.ThrowsExactly<ArgumentException>(
                () => SosFiltering.SosFiltFilt(sos, new[] { 1.0, 2.0 }, "odd", 20));
        }

        [TestMethod]
        public void Autocorrelation_RejectsBadLengths()
        {
            Assert.ThrowsExactly<ArgumentNullException>(
                () => Spectrum.Autocorrelation(null));

            double[] x = new double[100];
            Assert.ThrowsExactly<ArgumentException>(
                () => Spectrum.Autocorrelation(x, 60));
            Assert.ThrowsExactly<ArgumentException>(
                () => Spectrum.Autocorrelation(x, 128));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(
                () => Spectrum.Autocorrelation(x, 0));
        }

        [TestMethod]
        public void Normalize_RejectsZeroLeadingDenominator()
        {
            double[] b;
            double[] a;
            Assert.ThrowsExactly<ArgumentException>(
                () => TransferFunctionConversion.Normalize(new[] { 1.0 }, new[] { 0.0 }, out b, out a));
        }

        [TestMethod]
        public void AnalogSos_RejectsNonMinimalPairing()
        {
            ZpkModel zpk = AnalogPrototypes.Butterworth(2);
            Assert.ThrowsExactly<ArgumentException>(
                () => SecondOrderSections.Zpk2Sos(zpk.Zeros, zpk.Poles, zpk.Gain, true, "nearest"));
        }

        [TestMethod]
        public void UnknownPairing_IsRejected()
        {
            ZpkModel zpk = AnalogPrototypes.Butterworth(2);
            Assert.ThrowsExactly<ArgumentException>(
                () => SecondOrderSections.Zpk2Sos(zpk.Zeros, zpk.Poles, zpk.Gain, false, "bogus"));
        }
    }
}
