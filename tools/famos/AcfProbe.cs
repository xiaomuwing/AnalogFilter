using System;
using System.Globalization;
using ImcCoreLib;
using ImcFamosLib;

/// <summary>
/// 用受控输入探测 FAMOS 的 ACF 行为：输出长度、归一化方式、是否去均值。
/// </summary>
internal static class AcfProbe
{
    private static double[] ToArray(object data)
    {
        Array a = (Array)data;
        int lo = a.GetLowerBound(0);
        double[] r = new double[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = Convert.ToDouble(a.GetValue(lo + i));
        return r;
    }

    private static double[] Run(Famos famos, double[] x, double xDelta = 0.001)
    {
        var ch = new DChannel();
        ch.Data = (double[])x.Clone();
        ch.xDelta = xDelta;
        try
        {
            return ToArray(famos.Spec.ACF(ch).Data);
        }
        catch (Exception ex)
        {
            Console.WriteLine("    ACF 抛异常: " + ex.GetType().Name + ": " + ex.Message);
            return null;
        }
    }

    private static void Show(string label, double[] y, int count = 8)
    {
        if (y == null) { Console.WriteLine($"  {label}: <失败>"); return; }
        Console.Write($"  {label}: n={y.Length}  ");
        for (int i = 0; i < Math.Min(count, y.Length); i++)
            Console.Write($"{y[i]:G9} ");
        Console.WriteLine();
    }

    private static int Main()
    {
        var famos = new Famos();
        const int N = 1024;

        Console.WriteLine("=== 1. 长度与基本行为 ===");
        double[] impulse = new double[N];
        impulse[0] = 1.0;
        Show("冲激 [1,0,0,...]", Run(famos, impulse));

        double[] constant = new double[N];
        for (int i = 0; i < N; i++) constant[i] = 5.0;
        Show("常数 5", Run(famos, constant));

        Console.WriteLine();
        Console.WriteLine("=== 2. 去均值？(常数 vs 零均值) ===");
        double[] shifted = new double[N];
        for (int i = 0; i < N; i++) shifted[i] = 100.0 + Math.Sin(2 * Math.PI * 8 * i / N);
        Show("100 + sin(8 周期)", Run(famos, shifted), 5);

        Console.WriteLine();
        Console.WriteLine("=== 3. 归一化：单位方差白噪声 ===");
        var rnd = new Random(12345);
        double[] noise = new double[N];
        for (int i = 0; i < N; i++) noise[i] = rnd.NextDouble() - 0.5;
        double[] acfNoise = Run(famos, noise);
        if (acfNoise != null)
        {
            double mean = 0, var = 0;
            for (int i = 0; i < N; i++) mean += noise[i];
            mean /= N;
            for (int i = 0; i < N; i++) var += (noise[i] - mean) * (noise[i] - mean);
            var /= N;
            Console.WriteLine($"  输入: 均值={mean:G9} 方差(1/N)={var:G9}");
            Show("ACF", acfNoise, 6);
            Console.WriteLine($"  比值 ACF[0]/方差 = {acfNoise[0] / var:G9}");
            Console.WriteLine($"  比值 ACF[0]/(N*方差) = {acfNoise[0] / (N * var):G9}");
        }

        Console.WriteLine();
        Console.WriteLine("=== 4. 输出长度随输入长度变化 ===");
        foreach (int n in new[] { 64, 100, 128, 256, 1000 })
        {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) x[i] = Math.Sin(2 * Math.PI * 8 * i / n);
            double[] y = Run(famos, x);
            Console.WriteLine($"  输入 n={n,-5} -> 输出 n={(y == null ? "失败" : y.Length.ToString())}");
        }

        Console.WriteLine();
        Console.WriteLine("=== 5. 单频正弦：验证是否 FFT 路径（周期延拓假设）===");
        for (int cycles = 4; cycles <= 32; cycles *= 2)
        {
            double[] x = new double[N];
            for (int i = 0; i < N; i++) x[i] = Math.Cos(2 * Math.PI * cycles * i / N);
            double[] y = Run(famos, x);
            if (y == null) continue;
            double expect = 0.5 * Math.Cos(2 * Math.PI * cycles * 3 / N);
            Console.WriteLine($"  {cycles,2} 周期: ACF[0]={y[0]:G8}  ACF[3]={y[3]:G8}  " +
                              $"理论cos×0.5={expect:G8}");
        }

        return 0;
    }
}
