using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text;
using ImcCoreLib;
using ImcFamosLib;

/// <summary>
/// 全组合验证：4 种特征 × 4 种频带 × 多个阶数，FAMOS 与托管实现逐值对比。
/// 直接反射调用 ManagedFilter，不经过任何 FAMOS 回退。
/// </summary>
internal static class AllFamosProbe
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

        string[] names = { "巴特沃斯", "贝塞尔", "切比雪夫", "临界阻尼" };
        string[] bands = { "低通", "高通", "带通", "带阻" };

        Console.WriteLine($"输入 n={signal.Length} 峰峰值={pp:F3} RMS={rms:F3}");
        Console.WriteLine("低通/高通 50 Hz；带通/带阻 20–80 Hz；切比雪夫纹波 1 dB");
        Console.WriteLine();
        Console.WriteLine("特征       频带   阶数   max 偏差 %      rms 偏差 %");
        Console.WriteLine(new string('-', 62));

        double worst = 0;
        string worstAt = "-";
        int total = 0, bad = 0;

        foreach (int ch in new[] { 0, 1, 2, 3 })
        {
            foreach (int b in new[] { 0, 1, 2, 3 })
            {
                foreach (int order in new[] { 3, 4, 5, 6 })
                {
                    bool twoSided = b >= 2;
                    double param = ch == 2 ? 1.0 : 0.0;

                    var s1 = new DChannel(); s1.Data = (double[])signal.Clone(); s1.xDelta = 1.0 / fs;
                    DChannel reference;
                    switch (b)
                    {
                        case 0: reference = famos.Filter.LowPass(s1, ch, param, order, 50.0); break;
                        case 1: reference = famos.Filter.HighPass(s1, ch, param, order, 50.0); break;
                        case 2: reference = famos.Filter.BandPass(s1, ch, param, order, 20.0, 80.0); break;
                        default: reference = famos.Filter.BandStop(s1, ch, param, order, 20.0, 80.0); break;
                    }

                    var s2 = new DChannel(); s2.Data = (double[])signal.Clone(); s2.xDelta = 1.0 / fs;
                    DChannel got;
                    object[] call = twoSided
                        ? new object[] { s2, ch, param, order, 20.0, 80.0 }
                        : new object[] { s2, ch, param, order, 50.0 };
                    switch (b)
                    {
                        case 0: got = (DChannel)mLow.Invoke(null, call); break;
                        case 1: got = (DChannel)mHigh.Invoke(null, call); break;
                        case 2: got = (DChannel)mBand.Invoke(null, call); break;
                        default: got = (DChannel)mStop.Invoke(null, call); break;
                    }

                    double[] a = ToArray(reference.Data), c = ToArray(got.Data);
                    double mx = 0, ss = 0; int n = 0;
                    for (int i = 2000; i < a.Length - 2000; i++)
                    {
                        double d = a[i] - c[i];
                        mx = Math.Max(mx, Math.Abs(d));
                        ss += d * d;
                        n++;
                    }
                    double mPct = mx / pp * 100, rPct = Math.Sqrt(ss / n) / rms * 100;
                    total++;
                    if (mPct > worst) { worst = mPct; worstAt = names[ch] + "/" + bands[b] + "/" + order; }
                    if (mPct > 0.01) bad++;
                    Console.WriteLine($"{names[ch],-10} {bands[b],-5} {order,3}   {mPct,10:F6}      {rPct,10:F6}");
                }
            }
        }

        Console.WriteLine(new string('-', 62));
        Console.WriteLine($"共 {total} 组，超过 0.01% 的 {bad} 组，最大 {worst:F6}% ({worstAt})");
        return bad == 0 ? 0 : 1;
    }
}
