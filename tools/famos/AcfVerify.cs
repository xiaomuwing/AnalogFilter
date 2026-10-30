using System;
using System.Diagnostics;
using System.IO;
using AnalogFilter;

/// <summary>
/// 验证 Spectrum.Autocorrelation 与 FAMOS 的一致性，并检查边界行为与性能。
/// </summary>
internal static class AcfVerify
{
    private static double[] ReadBin(string path)
    {
        byte[] raw = File.ReadAllBytes(path);
        double[] r = new double[raw.Length / 8];
        Buffer.BlockCopy(raw, 0, r, 0, r.Length * 8);
        return r;
    }

    private static int Main(string[] args)
    {
        string dir = args[0];

        // ---- 1. 真实数据对比 FAMOS ----
        double[] x = ReadBin(Path.Combine(dir, "acf_input.bin"));
        double[] famos = ReadBin(Path.Combine(dir, "acf_famos.bin"));

        var sw = Stopwatch.StartNew();
        double[] ours = Spectrum.Autocorrelation(x);
        sw.Stop();

        Console.WriteLine($"=== 真实数据（n={x.Length}）===");
        Console.WriteLine($"  本库输出长度 {ours.Length}，FAMOS {famos.Length}");
        if (ours.Length == famos.Length)
        {
            double maxDiff = 0, sumSq = 0;
            int at = 0;
            for (int i = 0; i < ours.Length; i++)
            {
                double d = Math.Abs(ours[i] - famos[i]);
                if (d > maxDiff) { maxDiff = d; at = i; }
                sumSq += d * d;
            }
            Console.WriteLine($"  最大绝对差 {maxDiff:E3}（lag={at}）   RMS {Math.Sqrt(sumSq / ours.Length):E3}");
            Console.WriteLine($"  前 5 个: 本库 {ours[0]:G9} {ours[1]:G9} {ours[2]:G9} {ours[3]:G9} {ours[4]:G9}");
            Console.WriteLine($"          FAMOS {famos[0]:G9} {famos[1]:G9} {famos[2]:G9} {famos[3]:G9} {famos[4]:G9}");
        }
        Console.WriteLine($"  耗时 {sw.Elapsed.TotalMilliseconds:F2} ms");

        // ---- 2. 解析验证 ----
        Console.WriteLine();
        Console.WriteLine("=== 解析验证 ===");

        int n = 1024;
        double[] impulse = new double[n];
        impulse[0] = 1.0;
        double[] a1 = Spectrum.Autocorrelation(impulse);
        Console.Write("  冲激: ");
        for (int i = 0; i < 5; i++) Console.Write($"{a1[i]:G6} ");
        Console.WriteLine($"（期望 1 0 0 0 0）");

        double[] constant = new double[n];
        for (int i = 0; i < n; i++) constant[i] = 5.0;
        double[] a2 = Spectrum.Autocorrelation(constant);
        double maxConst = 0;
        for (int i = 0; i < n; i++) maxConst = Math.Max(maxConst, Math.Abs(a2[i] - 1.0));
        Console.WriteLine($"  常数 5: 与 1 的最大偏差 {maxConst:E3}（期望 0）");

        // 单频：ACF[k] = cos(2π·cycles·k/N)
        Console.WriteLine("  单频正弦（期望 ACF[0]=1, ACF[3]=cos(2π·cycles·3/N)）:");
        foreach (int cycles in new[] { 4, 8, 16, 32 })
        {
            double[] s = new double[n];
            for (int i = 0; i < n; i++) s[i] = Math.Cos(2 * Math.PI * cycles * i / n);
            double[] a = Spectrum.Autocorrelation(s);
            double expect = Math.Cos(2 * Math.PI * cycles * 3 / n);
            Console.WriteLine($"    {cycles,2} 周期: ACF[3]={a[3]:F9}  理论={expect:F9}  差={Math.Abs(a[3] - expect):E2}");
        }

        // 长度截断
        Console.WriteLine();
        Console.WriteLine("=== 长度处理（FAMOS：取最长 2 的幂前缀）===");
        foreach (int len in new[] { 64, 100, 128, 1000, 1024 })
        {
            double[] s = new double[len];
            for (int i = 0; i < len; i++) s[i] = Math.Sin(0.1 * i);
            Console.WriteLine($"  输入 {len,-5} -> 输出 {Spectrum.Autocorrelation(s).Length}");
        }

        // 全零
        Console.WriteLine();
        double[] zeros = new double[n];
        double[] a3 = Spectrum.Autocorrelation(zeros);
        Console.WriteLine($"=== 全零输入: 输出长度 {a3.Length}，[0]={a3[0]}（不抛异常）===");
        return 0;
    }
}
