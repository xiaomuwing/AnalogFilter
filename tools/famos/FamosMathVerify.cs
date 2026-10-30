using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text;
using ImcCoreLib;
using ImcFamosLib;

/// <summary>
/// 验证 ManagedBasics / ManagedSpectrum 与 FAMOS 的一致性。
/// </summary>
internal static class FamosMathVerify
{
    private static Assembly _asm;
    private static Type _basics;
    private static Type _spectrum;

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

    private static double[] Arr(object data)
    {
        Array a = (Array)data;
        int lo = a.GetLowerBound(0);
        double[] r = new double[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = Convert.ToDouble(a.GetValue(lo + i));
        return r;
    }

    private static DChannel Ch(double[] x, double dt = 0.001)
    {
        var c = new DChannel();
        c.Data = (double[])x.Clone();
        c.xDelta = dt;
        return c;
    }

    private static DChannel Call(string type, string method, DChannel ch)
    {
        Type t = type == "basics" ? _basics : _spectrum;
        MethodInfo m = t.GetMethod(method, BindingFlags.NonPublic | BindingFlags.Static);
        if (m == null) throw new InvalidOperationException("找不到方法 " + method);
        return (DChannel)m.Invoke(null, new object[] { ch });
    }

    private static void Report(string label, double[] reference, double[] actual, double scale)
    {
        if (reference.Length != actual.Length)
        {
            Console.WriteLine($"  {label,-10} 长度不一致: FAMOS={reference.Length} 本实现={actual.Length}");
            return;
        }
        double maxDiff = 0, sumSq = 0;
        bool hasNaN = false;
        for (int i = 0; i < reference.Length; i++)
        {
            if (double.IsNaN(reference[i]) && double.IsNaN(actual[i])) continue;
            if (double.IsNaN(reference[i]) || double.IsNaN(actual[i])) { hasNaN = true; continue; }
            double d = Math.Abs(reference[i] - actual[i]);
            if (d > maxDiff) maxDiff = d;
            sumSq += d * d;
        }
        double rel = scale > 0 ? maxDiff / scale * 100 : maxDiff;
        Console.WriteLine($"  {label,-10} n={reference.Length,-7} max|Δ|={maxDiff:E3}  RMS={Math.Sqrt(sumSq / reference.Length):E3}  " +
                          $"相对量程={rel:E2}%{(hasNaN ? "  [含 NaN 差异]" : "")}");
    }

    private static int Main(string[] args)
    {
        string bin = args[0];
        string csv = args[1];

        _asm = Assembly.LoadFrom(Path.Combine(bin, "FastAcquisition.FatigueAnalysis.exe"));
        _basics = _asm.GetType("FastAcquisition.FatigueAnalysis.ManagedBasics");
        _spectrum = _asm.GetType("FastAcquisition.FatigueAnalysis.ManagedSpectrum");

        double[] all = ReadCsv(csv, 2, 2);
        var famos = new Famos();

        // ---------- 1. Basics ----------
        int n = 4096;
        double[] x = new double[n];
        Array.Copy(all, x, n);
        double range = 0;
        foreach (double v in x) range = Math.Max(range, Math.Abs(v));

        Console.WriteLine($"=== 1. Basics.*（n={n}，跳过两端各 2000 点不做，全量对比）===");
        string[] ops = { "Sin", "Cos", "Tan", "Abs", "Exp", "Ln", "Log", "Sqr", "Sqrt", "Diff", "Int" };
        foreach (string op in ops)
        {
            DChannel reference;
            try
            {
                object basicsObject = famos.Basics;
                MethodInfo fm = basicsObject.GetType().GetMethod(op);
                reference = (DChannel)fm.Invoke(basicsObject, new object[] { Ch(x) });
            }
            catch (Exception ex)
            {
                Console.WriteLine($"  {op,-10} FAMOS 调用失败: {(ex.InnerException ?? ex).Message}");
                continue;
            }
            DChannel ours = Call("basics", op, Ch(x));
            Report(op, Arr(reference.Data), Arr(ours.Data), range);
        }

        // ASin / ACos 需要 [-1,1] 输入
        double[] unit = new double[n];
        for (int i = 0; i < n; i++) unit[i] = Math.Sin(0.001 * i) * 0.9;
        foreach (string op in new[] { "ASin", "ACos" })
        {
            object basicsObject2 = famos.Basics;
            var reference = (DChannel)basicsObject2.GetType().GetMethod(op).Invoke(basicsObject2, new object[] { Ch(unit) });
            var ours = Call("basics", op, Ch(unit));
            Report(op, Arr(reference.Data), Arr(ours.Data), 1.0);
        }

        // ---------- 2. Resample2 ----------
        Console.WriteLine();
        Console.WriteLine("=== 2. Resample2 ===");
        foreach (int len in new[] { 1000, 1024, 1500, 3000 })
        {
            double[] s = new double[len];
            for (int i = 0; i < len; i++) s[i] = Math.Sin(2 * Math.PI * 5 * i / len) * 10;
            var reference = famos.Misc.Resample2(Ch(s));
            var ours = Call("spectrum", "Resample2", Ch(s));
            double[] a = Arr(reference.Data), b = Arr(ours.Data);
            Console.WriteLine($"  n={len,-6} -> FAMOS n={a.Length,-6} dt={reference.xDelta:G9}   本实现 n={b.Length,-6} dt={ours.xDelta:G9}");
            Report("  值", a, b, 10.0);
        }

        // ---------- 3. FFTa ----------
        Console.WriteLine();
        Console.WriteLine("=== 3. FFTa（单边幅值谱）===");
        {
            int m = 1024;
            double[] s = new double[m];
            for (int i = 0; i < m; i++) s[i] = 3.0 + 2.0 * Math.Cos(2 * Math.PI * 64 * i / m) + 0.5 * Math.Sin(2 * Math.PI * 16 * i / m);
            var reference = famos.Spec.FFTa(Ch(s));
            var ours = Call("spectrum", "AmplitudeSpectrum", Ch(s));
            Report("合成信号", Arr(reference.Data), Arr(ours.Data), 3.0);
            double[] a = Arr(reference.Data), b = Arr(ours.Data);
            Console.WriteLine($"    DC: FAMOS={a[0]:G8} 本实现={b[0]:G8}   64 频点: FAMOS={a[64]:G8} 本实现={b[64]:G8}");
        }
        {
            double[] s = new double[n];
            Array.Copy(all, s, n);
            var reference = famos.Spec.FFTa(Ch(s));
            var ours = Call("spectrum", "AmplitudeSpectrum", Ch(s));
            Report("真实数据", Arr(reference.Data), Arr(ours.Data), range);
        }

        // ---------- 4. ACF（Resample2 + 自相关）----------
        Console.WriteLine();
        Console.WriteLine("=== 4. ACF（FAMOS: ACF∘Resample2）===");
        {
            int m = 65536;
            double[] s = new double[m];
            Array.Copy(all, s, m);
            var reference = famos.Spec.ACF(famos.Misc.Resample2(Ch(s)));
            var ours = Call("spectrum", "Autocorrelation", Ch(s));
            Report("真实数据", Arr(reference.Data), Arr(ours.Data), 1.0);
        }

        return 0;
    }
}
