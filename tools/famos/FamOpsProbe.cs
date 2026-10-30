using System;
using ImcCoreLib;
using ImcFamosLib;

/// <summary>
/// 探测 FAMOS 的 Basics.*、Misc.Resample2、Spec.FFTa 的确切语义。
/// </summary>
internal static class FamOpsProbe
{
    private static double[] Arr(object data)
    {
        Array a = (Array)data;
        int lo = a.GetLowerBound(0);
        double[] r = new double[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = Convert.ToDouble(a.GetValue(lo + i));
        return r;
    }

    private static DChannel Ch(double[] x, double xDelta = 0.001)
    {
        var ch = new DChannel();
        ch.Data = (double[])x.Clone();
        ch.xDelta = xDelta;
        return ch;
    }

    private static void Show(string label, double[] y, int count = 7)
    {
        if (y == null) { Console.WriteLine($"  {label,-12}: <失败>"); return; }
        Console.Write($"  {label,-12}: n={y.Length,-7} ");
        for (int i = 0; i < Math.Min(count, y.Length); i++) Console.Write($"{y[i]:G8} ");
        Console.WriteLine();
    }

    private static int Main()
    {
        var famos = new Famos();

        Console.WriteLine("=== 1. Basics.* 语义（输入 1..8，xDelta=0.001）===");
        double[] x = { 1, 2, 3, 4, 5, 6, 7, 8 };
        Show("输入", x);
        Show("Sin", Arr(famos.Basics.Sin(Ch(x)).Data));
        Show("Cos", Arr(famos.Basics.Cos(Ch(x)).Data));
        Show("Tan", Arr(famos.Basics.Tan(Ch(x)).Data));
        Show("Abs", Arr(famos.Basics.Abs(Ch(new double[] { -2, -1, 0, 1, 2 })).Data));
        Show("Exp", Arr(famos.Basics.Exp(Ch(x)).Data));
        Show("Ln", Arr(famos.Basics.Ln(Ch(x)).Data));
        Show("Log", Arr(famos.Basics.Log(Ch(x)).Data));
        Show("Sqr", Arr(famos.Basics.Sqr(Ch(x)).Data));
        Show("Sqrt", Arr(famos.Basics.Sqrt(Ch(x)).Data));
        Show("Diff", Arr(famos.Basics.Diff(Ch(x)).Data), 8);
        Show("Int", Arr(famos.Basics.Int(Ch(x)).Data), 8);
        Console.WriteLine("  ln(1..8)   = 0 0.6931 1.0986 1.3863 1.6094 1.7918 1.9459 2.0794");
        Console.WriteLine("  log10(1..8)= 0 0.3010 0.4771 0.6021 0.6990 0.7782 0.8451 0.9031");

        Console.WriteLine();
        Console.WriteLine("=== 2. Diff/Int 是否乘除 xDelta ===");
        Show("Diff dt=0.001", Arr(famos.Basics.Diff(Ch(x, 0.001)).Data), 5);
        Show("Diff dt=1    ", Arr(famos.Basics.Diff(Ch(x, 1.0)).Data), 5);
        Show("Int  dt=0.001", Arr(famos.Basics.Int(Ch(x, 0.001)).Data), 5);
        Show("Int  dt=1    ", Arr(famos.Basics.Int(Ch(x, 1.0)).Data), 5);

        Console.WriteLine();
        Console.WriteLine("=== 3. Misc.Resample2：长度行为 ===");
        foreach (int n in new[] { 100, 1000, 1024, 1500, 3000 })
        {
            double[] s = new double[n];
            for (int i = 0; i < n; i++) s[i] = Math.Sin(2 * Math.PI * 5 * i / n);
            var r = famos.Misc.Resample2(Ch(s));
            Console.WriteLine($"  输入 n={n,-6} -> 输出 n={Arr(r.Data).Length,-6} xDelta={r.xDelta:G9}");
        }

        Console.WriteLine();
        Console.WriteLine("=== 4. Misc.Resample2：插值方式（斜坡输入便于识别）===");
        {
            int n = 100;
            double[] s = new double[n];
            for (int i = 0; i < n; i++) s[i] = i;
            var r = Arr(famos.Misc.Resample2(Ch(s, 1.0)).Data);
            Console.WriteLine($"  输入 0..99 (n=100) -> 输出 n={r.Length}");
            Console.Write("  前 8 个: ");
            for (int i = 0; i < 8; i++) Console.Write($"{r[i]:G8} ");
            Console.WriteLine();
            Console.Write("  后 4 个: ");
            for (int i = r.Length - 4; i < r.Length; i++) Console.Write($"{r[i]:G8} ");
            Console.WriteLine();
        }

        Console.WriteLine();
        Console.WriteLine("=== 5. Spec.FFTa：幅值标定 ===");
        {
            int n = 1024, cycles = 32;
            double[] s = new double[n];
            for (int i = 0; i < n; i++) s[i] = Math.Cos(2 * Math.PI * cycles * i / n);
            var y = Arr(famos.Spec.FFTa(Ch(s, 1.0 / n)).Data);
            Console.WriteLine($"  输入 cos({cycles} 周期, n={n}) -> 输出 n={y.Length}");
            int peak = 0;
            for (int i = 1; i < y.Length; i++)
            {
                if (y[i] > y[peak]) peak = i;
            }
            Console.WriteLine($"  峰值 index={peak}  值={y[peak]:G10}");
            Console.WriteLine($"  候选标定: |X|/N={1.0 / n:G6}  2|X|/N={2.0 / n:G6}  |X|={n / 2.0:G6}");
            Console.WriteLine($"  DC 项 y[0]={y[0]:G10}");
            Console.Write("  前 5 个: ");
            for (int i = 0; i < 5; i++) Console.Write($"{y[i]:G8} ");
            Console.WriteLine();
        }

        return 0;
    }
}
