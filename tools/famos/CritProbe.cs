using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using AnalogFilter;

/// <summary>
/// 临界阻尼专项验证：FAMOS COM 与 AnalogFilter 在同一份数据上的逐值对比。
/// 带通/带阻按 FAMOS 的"总极点数"约定，原型阶数取一半。
/// </summary>
internal static class CritProbe
{
    private static double[] ReadCsvChannel(string path, int columnIndex, int skipLines)
    {
        List<double> values = new List<double>();
        using (StreamReader reader = new StreamReader(path, Encoding.Default))
        {
            for (int i = 0; i < skipLines; i++)
            {
                reader.ReadLine();
            }
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length == 0) continue;
                string[] parts = line.Split(',');
                if (parts.Length <= columnIndex) continue;
                values.Add(double.Parse(parts[columnIndex].Trim(), CultureInfo.InvariantCulture));
            }
        }
        return values.ToArray();
    }

    private static double[] ToArray(object comData)
    {
        Array array = (Array)comData;
        int lower = array.GetLowerBound(0);
        double[] result = new double[array.Length];
        for (int i = 0; i < array.Length; i++)
        {
            result[i] = Convert.ToDouble(array.GetValue(lower + i));
        }
        return result;
    }

    private static void Compare(string label, double[] reference, double[] candidate, double pp, double rms, int edge)
    {
        int start = edge;
        int stop = reference.Length - edge;
        double maxDiff = 0.0;
        double sumSq = 0.0;
        int count = 0;
        for (int i = start; i < stop; i++)
        {
            double d = reference[i] - candidate[i];
            maxDiff = Math.Max(maxDiff, Math.Abs(d));
            sumSq += d * d;
            count++;
        }
        double rmsDiff = Math.Sqrt(sumSq / count);
        double maxPct = maxDiff / pp * 100.0;
        double rmsPct = rmsDiff / rms * 100.0;
        string flag = maxPct < 1e-6 ? "完全一致" : (maxPct < 0.01 ? "一致" : "有差异");
        Console.WriteLine($"  {label,-46} max={maxPct,8:F4}%  rms={rmsPct,8:F4}%   {flag}");
    }

    private static int Main(string[] args)
    {
        string csvPath = args[0];
        double fs = 1000.0;
        double[] signal = ReadCsvChannel(csvPath, 2, 2);
        Console.WriteLine($"输入 Signal: n={signal.Length}  fs={fs} Hz");

        double min = double.MaxValue, max = double.MinValue, sumSq = 0.0;
        foreach (double v in signal)
        {
            min = Math.Min(min, v);
            max = Math.Max(max, v);
            sumSq += v * v;
        }
        double pp = max - min;
        double rms = Math.Sqrt(sumSq / signal.Length);
        Console.WriteLine($"峰峰值={pp:F6}  RMS={rms:F6}");
        Console.WriteLine();

        var famos = new ImcFamosLib.Famos();

        Console.WriteLine("=== 临界阻尼：FAMOS vs AnalogFilter（跳过两端各 2000 点）===");

        // 低通 / 高通：原型阶数 = FAMOS 阶数
        foreach (bool lowPass in new[] { true, false })
        {
            int order = 4;
            double freq = 50.0;
            var ch = new ImcCoreLib.DChannel();
            ch.Data = (double[])signal.Clone();
            ch.xDelta = 1.0 / fs;

            ImcCoreLib.DChannel famosResult = lowPass
                ? famos.Filter.LowPass(ch, 3, 0.0, order, freq)
                : famos.Filter.HighPass(ch, 3, 0.0, order, freq);

            double[] reference = ToArray(famosResult.Data);
            double[][] sos = FilterDesign.DesignSos(
                order, new[] { freq / (fs / 2.0) },
                FilterType.CriticallyDamped,
                lowPass ? BandType.LowPass : BandType.HighPass);
            double[] candidate = SosFiltering.SosFilt(sos, signal);

            Compare($"临界阻尼 {(lowPass ? "低通" : "高通")} N={order} {freq}Hz",
                    reference, candidate, pp, rms, 2000);
        }

        // 带通 / 带阻：FAMOS 阶数 = 总极点数，原型阶数取一半
        foreach (bool bandPass in new[] { true, false })
        {
            int famosOrder = 4;
            int prototypeOrder = famosOrder / 2;
            double f1 = 20.0, f2 = 80.0;

            var ch = new ImcCoreLib.DChannel();
            ch.Data = (double[])signal.Clone();
            ch.xDelta = 1.0 / fs;

            ImcCoreLib.DChannel famosResult = bandPass
                ? famos.Filter.BandPass(ch, 3, 0.0, famosOrder, f1, f2)
                : famos.Filter.BandStop(ch, 3, 0.0, famosOrder, f1, f2);

            double[] reference = ToArray(famosResult.Data);
            double[][] sos = FilterDesign.DesignSos(
                prototypeOrder, new[] { f1 / (fs / 2.0), f2 / (fs / 2.0) },
                FilterType.CriticallyDamped,
                bandPass ? BandType.BandPass : BandType.BandStop);
            double[] candidate = SosFiltering.SosFilt(sos, signal);

            Compare($"临界阻尼 {(bandPass ? "带通" : "带阻")} FAMOS阶={famosOrder}→原型{prototypeOrder} {f1}/{f2}Hz",
                    reference, candidate, pp, rms, 2000);
        }

        Console.WriteLine();
        Console.WriteLine("=== 对照：更高阶（N=8）低通 ===");
        {
            int order = 8;
            double freq = 50.0;
            var ch = new ImcCoreLib.DChannel();
            ch.Data = (double[])signal.Clone();
            ch.xDelta = 1.0 / fs;
            double[] reference = ToArray(famos.Filter.LowPass(ch, 3, 0.0, order, freq).Data);
            double[][] sos = FilterDesign.DesignSos(order, new[] { freq / (fs / 2.0) },
                FilterType.CriticallyDamped, BandType.LowPass);
            double[] candidate = SosFiltering.SosFilt(sos, signal);
            Compare($"临界阻尼 低通 N={order} {freq}Hz", reference, candidate, pp, rms, 2000);
        }

        return 0;
    }
}
