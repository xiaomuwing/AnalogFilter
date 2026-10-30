using System;
using System.Diagnostics;
using System.IO;
using FastAcquisition.Analysis.Common;
using MathNet.Numerics.LinearAlgebra;

/// <summary>
/// 验证重写后的 SavitzkyGolayFilter：
/// 1) 与旧算法结果一致；2) 与 SciPy savgol_filter(mode='interp') 一致；3) 性能与分配；4) 参数校验。
/// </summary>
internal static class SgVerify
{
    /// <summary>旧实现（原样保留，用于对照）。</summary>
    private sealed class Legacy
    {
        private readonly int sidePoints;
        private Matrix<double> coefficients;

        internal Legacy(int sidePoints, int polynomialOrder)
        {
            this.sidePoints = sidePoints;
            double[,] a = new double[(sidePoints << 1) + 1, polynomialOrder + 1];
            for (int m = -sidePoints; m <= sidePoints; m++)
                for (int i = 0; i <= polynomialOrder; i++)
                    a[m + sidePoints, i] = Math.Pow(m, i);
            Matrix<double> s = Matrix<double>.Build.DenseOfArray(a);
            coefficients = s.Multiply(s.TransposeThisAndMultiply(s).Inverse()).Multiply(s.Transpose());
        }

        internal double[] Process(double[] samples)
        {
            int length = samples.Length;
            double[] output = new double[length];
            int frameSize = (sidePoints << 1) + 1;
            double[] frame = new double[frameSize];
            Array.Copy(samples, frame, frameSize);
            for (int i = 0; i < sidePoints; ++i)
                output[i] = coefficients.Column(i).DotProduct(Vector<double>.Build.DenseOfArray(frame));
            for (int n = sidePoints; n < length - sidePoints; ++n)
            {
                Array.ConstrainedCopy(samples, n - sidePoints, frame, 0, frameSize);
                output[n] = coefficients.Column(sidePoints).DotProduct(Vector<double>.Build.DenseOfArray(frame));
            }
            Array.ConstrainedCopy(samples, length - frameSize, frame, 0, frameSize);
            for (int i = 0; i < sidePoints; ++i)
                output[length - sidePoints + i] = coefficients.Column(sidePoints + 1 + i).DotProduct(Vector<double>.Build.Dense(frame));
            return output;
        }
    }

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
        double[] input = ReadBin(Path.Combine(dir, "sg_input.bin"));
        Console.WriteLine($"输入 n={input.Length}");
        Console.WriteLine();

        int[][] cases = { new[] { 10, 3 }, new[] { 25, 2 }, new[] { 50, 4 } };
        Console.WriteLine("=== 精度：新实现 vs 旧算法 vs SciPy savgol_filter(mode='interp') ===");
        foreach (int[] c in cases)
        {
            int side = c[0], order = c[1];
            double[] reference = ReadBin(Path.Combine(dir, $"sg_ref_{side}_{order}.bin"));

            double[] fresh = new SavitzkyGolayFilter(side, order).Process(input);
            double[] legacy = new Legacy(side, order).Process(input);

            double dNew = 0, dOld = 0;
            for (int i = 0; i < input.Length; i++)
            {
                dNew = Math.Max(dNew, Math.Abs(fresh[i] - reference[i]));
                dOld = Math.Max(dOld, Math.Abs(legacy[i] - reference[i]));
            }
            Console.WriteLine($"  side={side,2} order={order} 窗口={2 * side + 1,3}  新实现 vs SciPy: {dNew:E3}   旧算法 vs SciPy: {dOld:E3}");
        }

        Console.WriteLine();
        Console.WriteLine("=== 性能与分配（n=300000）===");
        foreach (int[] c in cases)
        {
            int side = c[0], order = c[1];
            var newFilter = new SavitzkyGolayFilter(side, order);
            var oldFilter = new Legacy(side, order);

            newFilter.Process(input);
            oldFilter.Process(input);   // 预热

            long a0 = GC.GetAllocatedBytesForCurrentThread();
            var sw = Stopwatch.StartNew();
            newFilter.Process(input);
            sw.Stop();
            double tNew = sw.Elapsed.TotalMilliseconds;
            long allocNew = GC.GetAllocatedBytesForCurrentThread() - a0;

            long a1 = GC.GetAllocatedBytesForCurrentThread();
            sw.Restart();
            oldFilter.Process(input);
            sw.Stop();
            double tOld = sw.Elapsed.TotalMilliseconds;
            long allocOld = GC.GetAllocatedBytesForCurrentThread() - a1;

            Console.WriteLine($"  窗口={2 * side + 1,3}: 旧 {tOld,7:F2} ms / {allocOld / 1024.0 / 1024.0,7:F1} MB   " +
                              $"新 {tNew,6:F2} ms / {allocNew / 1024.0,6:F1} B   加速 {tOld / Math.Max(tNew, 0.001),5:F1}×");
        }

        Console.WriteLine();
        Console.WriteLine("=== 参数校验（应全部抛出明确异常）===");
        Check("sidePoints = -1", () => new SavitzkyGolayFilter(-1, 2));
        Check("polynomialOrder = -1", () => new SavitzkyGolayFilter(3, -1));
        Check("polynomialOrder >= 窗口长度", () => new SavitzkyGolayFilter(5, 11));
        Check("samples = null", () => new SavitzkyGolayFilter(3, 2).Process(null));
        Check("样本长度 < 窗口长度", () => new SavitzkyGolayFilter(50, 2).Process(new double[10]));
        return 0;
    }

    private static void Check(string label, Action action)
    {
        try
        {
            action();
            Console.WriteLine($"  [未拦截] {label}");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"  [OK] {label,-24} -> {ex.GetType().Name}: {ex.Message}");
        }
    }
}
