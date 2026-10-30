using System;
using System.Collections.Generic;
using System.Globalization;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// <c>SPEC</c> 行解析。格式：
    /// <c>digital|analog &lt;类型&gt; &lt;频带&gt; &lt;阶数&gt; &lt;wn...&gt; &lt;rp&gt; &lt;rs&gt;</c>
    /// 或 <c>besselp &lt;阶数&gt;</c>、<c>acf &lt;长度&gt;</c>。
    /// </summary>
    internal sealed class DesignSpec
    {
        internal string Kind { get; private set; }
        internal FilterType Type { get; private set; }
        internal BandType Band { get; private set; }
        internal int Order { get; private set; }
        internal double[] Wn { get; private set; }
        internal double Rp { get; private set; }
        internal double Rs { get; private set; }

        internal static DesignSpec Parse(string[] spec)
        {
            DesignSpec result = new DesignSpec();
            result.Kind = spec[0];

            if (result.Kind == "digital" || result.Kind == "analog")
            {
                result.Type = ParseType(spec[1]);
                result.Band = ParseBand(spec[2]);
                result.Order = int.Parse(spec[3], CultureInfo.InvariantCulture);

                int count = spec.Length - 6;
                result.Wn = new double[count];
                for (int i = 0; i < count; i++)
                {
                    result.Wn[i] = double.Parse(spec[4 + i], CultureInfo.InvariantCulture);
                }
                result.Rp = double.Parse(spec[spec.Length - 2], CultureInfo.InvariantCulture);
                result.Rs = double.Parse(spec[spec.Length - 1], CultureInfo.InvariantCulture);
                return result;
            }

            if (result.Kind == "besselp" || result.Kind == "acf")
            {
                result.Order = int.Parse(spec[1], CultureInfo.InvariantCulture);
                return result;
            }

            throw new NotSupportedException("未知的 SPEC 类型：" + result.Kind);
        }

        private static FilterType ParseType(string name)
        {
            switch (name)
            {
                case "butter": return FilterType.Butterworth;
                case "cheby1": return FilterType.ChebyshevI;
                case "cheby2": return FilterType.ChebyshevII;
                case "bessel": return FilterType.Bessel;
                default: throw new NotSupportedException("未知的滤波器类型：" + name);
            }
        }

        private static BandType ParseBand(string name)
        {
            switch (name)
            {
                case "lp": return BandType.LowPass;
                case "hp": return BandType.HighPass;
                case "bp": return BandType.BandPass;
                case "bs": return BandType.BandStop;
                default: throw new NotSupportedException("未知的频带：" + name);
            }
        }

        /// <summary>把用例 id 收集成 <c>[DynamicData]</c> 需要的形状。</summary>
        internal static IEnumerable<object[]> IdsOf(string kind)
        {
            foreach (ReferenceCase item in ReferenceData.All)
            {
                if (item.Spec != null && item.Spec[0] == kind)
                {
                    yield return new object[] { item.Id };
                }
            }
        }
    }
}
