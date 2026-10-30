using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text;
using ImcCoreLib;
using ImcFamosLib;

/// <summary>
/// 验证 FastAcquisition 中 ManagedFilter 对"切比雪夫"四种频带的 FAMOS 对齐程度。
/// </summary>
internal static class ChebFamosProbe
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
        MethodInfo lp = t.GetMethod("LowPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo hp = t.GetMethod("HighPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo bp = t.GetMethod("BandPass", BindingFlags.NonPublic | BindingFlags.Static);
        MethodInfo bs = t.GetMethod("BandStop", BindingFlags.NonPublic | BindingFlags.Static);


        double[] signal = ReadCsv(csv, 2, 2);
        double fs = 1000.0;
        double min = double.MaxValue, max = double.MinValue, sq = 0;
        foreach (double v in signal) { min = Math.Min(min, v); max = Math.Max(max, v); sq += v * v; }
        double pp = max - min, rms = Math.Sqrt(sq / signal.Length);

        var famos = new ImcFamosLib.Famos();
        Console.WriteLine($"输入 n={signal.Length} 峰峰值={pp:F3} RMS={rms:F3}");
        Console.WriteLine();
        Console.WriteLine("=== 切比雪夫（characteristic=2, parameter=1 dB）经 ManagedFilter vs FAMOS ===");

        (string Name, bool IsLow, bool IsBp, int Order)[] cases =
        {
            ("低通", true, false, 4), ("高通", false, false, 4),
            ("带通", false, true, 4), ("带阻", false, false, 0),
        };

        foreach (var c in cases)
        {
          try
          {
            string label = c.Name;
            bool twoSided = label == "带通" || label == "带阻";
            int famosOrder = 4;

            var src1 = new DChannel(); src1.Data = (double[])signal.Clone(); src1.xDelta = 1.0 / fs;
            DChannel reference;
            if (label == "低通") reference = famos.Filter.LowPass(src1, 2, 1.0, famosOrder, 50.0);
            else if (label == "高通") reference = famos.Filter.HighPass(src1, 2, 1.0, famosOrder, 50.0);
            else if (label == "带通") reference = famos.Filter.BandPass(src1, 2, 1.0, famosOrder, 20.0, 80.0);
            else reference = famos.Filter.BandStop(src1, 2, 1.0, famosOrder, 20.0, 80.0);

            var src2 = new DChannel(); src2.Data = (double[])signal.Clone(); src2.xDelta = 1.0 / fs;
            DChannel got;
            if (label == "低通") got = (DChannel)lp.Invoke(null, new object[] { src2, 2, 1.0, famosOrder, 50.0 });
            else if (label == "高通") got = (DChannel)hp.Invoke(null, new object[] { src2, 2, 1.0, famosOrder, 50.0 });
            else if (label == "带通") got = (DChannel)bp.Invoke(null, new object[] { src2, 2, 1.0, famosOrder, 20.0, 80.0 });
            else got = (DChannel)bs.Invoke(null, new object[] { src2, 2, 1.0, famosOrder, 20.0, 80.0 });

            double[] a = ToArray(reference.Data), b = ToArray(got.Data);
            double mx = 0, ss = 0; int n = 0;
            for (int i = 2000; i < a.Length - 2000; i++)
            {
                double d = a[i] - b[i];
                mx = Math.Max(mx, Math.Abs(d));
                ss += d * d;
                n++;
            }
            Console.WriteLine($"  {label}  max={mx / pp * 100,8:F5}%  rms={Math.Sqrt(ss / n) / rms * 100,8:F5}%");
          }
          catch (Exception ex)
          {
            Exception e = ex is TargetInvocationException && ex.InnerException != null ? ex.InnerException : ex;
            Console.WriteLine($"  {c.Name}  异常: {e.GetType().Name}: {e.Message}");
          }
        }
        return 0;
    }
}

