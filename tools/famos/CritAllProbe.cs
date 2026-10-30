using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Text;
using ImcCoreLib;
using ImcFamosLib;
using AnalogFilter;

/// <summary>
/// 临界阻尼专项验证：四种频带 × 多个阶数，FAMOS 与托管实现逐值对比。
/// 同时验证归一化做在原型上之后，带通/带阻的 -3 dB 带边是否落位正确。
/// </summary>
internal static class CritAllProbe
{
    private static double[] ReadCsv(string path, int column, int skip)
    {
        var list = new List<double>();
        using (var reader = new StreamReader(path, Encoding.Default))
        {
            for (int i = 0; i < skip; i++) reader.ReadLine();
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length == 0) continue;
                string[] p = line.Split(',');
                if (p.Length <= column) continue;
                list.Add(double.Parse(p[column].Trim(), CultureInfo.InvariantCulture));
            }
        }
        return list.ToArray();
    }

    private static double[] ToArray(object data)
    {
        Array a = (Array)data;
        int lo = a.GetLowerBound(0);
        double[] r = new double[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = Convert.ToDouble(a.GetValue(lo + i));
        return r;
    }

    private static int Main(string[] args)
    {
        string bin = args[0];
        string csv = args[1];

        Assembly asm = Assembly.LoadFrom(Path.Combine(bin, "FastAcquisition.FatigueAnalysis.exe"));
        Type t = asm.GetType("FastAcquisition.FatigueAnalysis.ManagedFilter");
        MethodInfo mLow = t.GetMethod("LowPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo mHigh = t.GetMethod("HighPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo mBand = t.GetMethod("BandPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo mStop = t.GetMethod("BandStop", BindingFlags.NonPublic | BindingFlags.Static);


        double[] signal = ReadCsv(csv, 2, 2);
        double fs = 1000.0;
        double min = double.MaxValue, max = double.MinValue, sq = 0;
        foreach (double v in signal) { min = Math.Min(min, v); max = Math.Max(max, v); sq += v * v; }
        double pp = max - min, rms = Math.Sqrt(sq / signal.Length);
        var famos = new Famos();

        Console.WriteLine($"输入 n={signal.Length} 峰峰值={pp:F3} RMS={rms:F3}");
        Console.WriteLine();
        Console.WriteLine("=== 临界阻尼：FAMOS vs ManagedFilter（跳过两端各 2000 点）===");
        Console.WriteLine("  低通/高通 50 Hz；带通/带阻 20–80 Hz");

        foreach (int order in new[] { 2, 3, 4, 5, 6, 8 })
        {
            (string Name, bool Low, bool High, bool Bp, bool Bs) = ("低通", true, false, false, false);
            foreach (var c in new[] { ("低通", 0), ("高通", 1), ("带通", 2), ("带阻", 3) })
            {
                var src1 = new DChannel(); src1.Data = (double[])signal.Clone(); src1.xDelta = 1.0 / fs;
                DChannel reference;
                switch (c.Item2)
                {
                    case 0: reference = famos.Filter.LowPass(src1, 3, 0.0, order, 50.0); break;
                    case 1: reference = famos.Filter.HighPass(src1, 3, 0.0, order, 50.0); break;
                    case 2: reference = famos.Filter.BandPass(src1, 3, 0.0, order, 20.0, 80.0); break;
                    default: reference = famos.Filter.BandStop(src1, 3, 0.0, order, 20.0, 80.0); break;
                }

                var src2 = new DChannel(); src2.Data = (double[])signal.Clone(); src2.xDelta = 1.0 / fs;
                DChannel got;
                switch (c.Item2)
                {
                    case 0: got = (DChannel)mLow.Invoke(null, new object[] { src2, 3, 0.0, order, 50.0 }); break;
                    case 1: got = (DChannel)mHigh.Invoke(null, new object[] { src2, 3, 0.0, order, 50.0 }); break;
                    case 2: got = (DChannel)mBand.Invoke(null, new object[] { src2, 3, 0.0, order, 20.0, 80.0 }); break;
                    default: got = (DChannel)mStop.Invoke(null, new object[] { src2, 3, 0.0, order, 20.0, 80.0 }); break;
                }

                double[] a = ToArray(reference.Data), b = ToArray(got.Data);
                double mx = 0, ss = 0; int n = 0;
                for (int i = 2000; i < a.Length - 2000; i++)
                {
                    double d = a[i] - b[i];
                    mx = Math.Max(mx, Math.Abs(d));
                    ss += d * d;
                    n++;
                }
                Console.WriteLine($"  FAMOS阶={order} {c.Item1}  max={mx / pp * 100,9:F5}%  rms={Math.Sqrt(ss / n) / rms * 100,9:F5}%");
            }
            Console.WriteLine();
        }

        // 带边 -3 dB 自检（不依赖 FAMOS）
        Console.WriteLine("=== 带通/带阻 -3 dB 带边自检（目标 20 / 80 Hz）===");
        foreach (int prototypeOrder in new[] { 1, 2, 3 })
        {
            foreach (bool isBandPass in new[] { true, false })
            {
                double[][] sos = FilterDesign.DesignSos(
                    prototypeOrder, new[] { 20.0 / (fs / 2.0), 80.0 / (fs / 2.0) },
                    FilterType.CriticallyDamped,
                    isBandPass ? BandType.BandPass : BandType.BandStop);

                double[] w = new double[40000];
                for (int i = 0; i < w.Length; i++) w[i] = 1.0 + i * 0.005;
                Complex[] h = FrequencyResponse.Freqz(sos, w.Select(v => 2 * Math.PI * v / fs).ToArray());
                double[] mag = FrequencyResponse.Magnitude(h);

                // 带通找两侧穿越 0.7071 的点；带阻找中心附近的最低点
                var crossings = new List<double>();
                for (int i = 1; i < mag.Length; i++)
                {
                    bool up = mag[i - 1] < 0.7071067811865476 && mag[i] >= 0.7071067811865476;
                    bool down = mag[i - 1] >= 0.7071067811865476 && mag[i] < 0.7071067811865476;
                    if (up || down) crossings.Add(w[i]);
                }
                string text = crossings.Count >= 2
                    ? $"{crossings[0],7:F2} / {crossings[1],7:F2} Hz"
                    : (crossings.Count == 1 ? $"{crossings[0],7:F2} Hz（仅一处）" : "无穿越");
                Console.WriteLine($"  原型N={prototypeOrder} {(isBandPass ? "带通" : "带阻")}: -3dB 穿越 {text}");
            }
        }
        return 0;
    }
}


