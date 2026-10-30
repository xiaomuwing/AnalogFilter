using System;
using System.Diagnostics;
using FastAcquisition.Analysis.Common;

/// <summary>
/// 对比现有 SavitzkyGolayFilter 与优化实现的耗时、内存分配与结果一致性。
/// </summary>
internal static class SgBench
{
    /// <summary>优化版：设计期一次性算出全部权重行，运行期纯滑动点积、零分配。</summary>
    private sealed class Optimized
    {
        private readonly int sidePoints;
        private readonly int frameSize;
        private readonly double[] centerKernel;
        private readonly double[][] rows;

        internal Optimized(int sidePoints, int polynomialOrder)
        {
            this.sidePoints = sidePoints;
            frameSize = (sidePoints << 1) + 1;

            int cols = polynomialOrder + 1;
            double[,] a = new double[frameSize, cols];
            for (int m = -sidePoints; m <= sidePoints; m++)
                for (int i = 0; i < cols; i++)
                    a[m + sidePoints, i] = Math.Pow(m, i);

            // 正规方程 AᵀA 及其逆
            double[,] ata = new double[cols, cols];
            for (int i = 0; i < cols; i++)
                for (int j = 0; j < cols; j++)
                {
                    double s = 0;
                    for (int m = 0; m < frameSize; m++) s += a[m, i] * a[m, j];
                    ata[i, j] = s;
                }
            double[,] inv = Invert(ata, cols);

            // 帽子矩阵 H[k][m] = Σ_i Σ_j a[k,i] * inv[i,j] * a[m,j]
            rows = new double[frameSize][];
            for (int k = 0; k < frameSize; k++)
            {
                double[] w = new double[frameSize];
                for (int m = 0; m < frameSize; m++)
                {
                    double s = 0;
                    for (int i = 0; i < cols; i++)
                        for (int j = 0; j < cols; j++)
                            s += a[k, i] * inv[i, j] * a[m, j];
                    w[m] = s;
                }
                rows[k] = w;
            }
            centerKernel = rows[sidePoints];
        }

        private static double[,] Invert(double[,] m, int n)
        {
            double[,] inv = new double[n, n];
            for (int i = 0; i < n; i++) inv[i, i] = 1.0;
            for (int i = 0; i < n; i++)
            {
                int piv = i;
                for (int r = i + 1; r < n; r++) if (Math.Abs(m[r, i]) > Math.Abs(m[piv, i])) piv = r;
                if (piv != i)
                    for (int c = 0; c < n; c++)
                    {
                        double t = m[i, c]; m[i, c] = m[piv, c]; m[piv, c] = t;
                        t = inv[i, c]; inv[i, c] = inv[piv, c]; inv[piv, c] = t;
                    }
                double d = m[i, i];
                for (int c = 0; c < n; c++) { m[i, c] /= d; inv[i, c] /= d; }
                for (int r = 0; r < n; r++)
                {
                    if (r == i) continue;
                    double f = m[r, i];
                    if (f == 0) continue;
                    for (int c = 0; c < n; c++) { m[r, c] -= f * m[i, c]; inv[r, c] -= f * inv[i, c]; }
                }
            }
            return inv;
        }

        internal double[] Process(double[] samples)
        {
            int n = samples.Length;
            int half = sidePoints;
            double[] output = new double[n];

            // 左边界：用首个窗口 + 对应的边界权重行
            for (int k = 0; k < half; k++)
            {
                double s = 0;
                double[] w = rows[k];
                for (int j = 0; j < frameSize; j++) s += samples[j] * w[j];
                output[k] = s;
            }

            // 中间：固定核滑动，利用对称性把乘法减半
            for (int i = half; i < n - half; i++)
            {
                double s = centerKernel[half] * samples[i];
                for (int j = 1; j <= half; j++)
                    s += centerKernel[half + j] * (samples[i + j] + samples[i - j]);
                output[i] = s;
            }

            // 右边界
            int start = n - frameSize;
            for (int k = 0; k < half; k++)
            {
                double s = 0;
                double[] w = rows[half + 1 + k];
                for (int j = 0; j < frameSize; j++) s += samples[start + j] * w[j];
                output[n - half + k] = s;
            }
            return output;
        }
    }

    private static void Main()
    {
        int n = 300000;
        double[] data = new double[n];
        var rnd = new Random(12345);
        for (int i = 0; i < n; i++)
            data[i] = Math.Sin(2 * Math.PI * 0.01 * i)
                    + 0.3 * Math.Sin(2 * Math.PI * 0.17 * i)
                    + 0.05 * (rnd.NextDouble() - 0.5);

        int[][] cases = { new[] { 10, 3 }, new[] { 25, 2 }, new[] { 50, 4 } };
        foreach (int[] c in cases)
        {
            int side = c[0], order = c[1];
            Console.WriteLine($"=== sidePoints={side}（窗口 {2 * side + 1} 点）阶数={order}  n={n} ===");

            var sw = Stopwatch.StartNew();
            var oldFilter = new SavitzkyGolayFilter(side, order);
            sw.Stop();
            double tDesignOld = sw.Elapsed.TotalMilliseconds;

            long a0 = GC.GetAllocatedBytesForCurrentThread();
            sw.Restart();
            double[] y1 = oldFilter.Process(data);
            sw.Stop();
            double tOld = sw.Elapsed.TotalMilliseconds;
            long allocOld = GC.GetAllocatedBytesForCurrentThread() - a0;

            sw.Restart();
            var newFilter = new Optimized(side, order);
            sw.Stop();
            double tDesignNew = sw.Elapsed.TotalMilliseconds;

            long a1 = GC.GetAllocatedBytesForCurrentThread();
            sw.Restart();
            double[] y2 = newFilter.Process(data);
            sw.Stop();
            double tNew = sw.Elapsed.TotalMilliseconds;
            long allocNew = GC.GetAllocatedBytesForCurrentThread() - a1;

            double maxDiff = 0;
            for (int i = 0; i < n; i++) maxDiff = Math.Max(maxDiff, Math.Abs(y1[i] - y2[i]));

            Console.WriteLine($"  设计:  现实现 {tDesignOld,8:F3} ms   优化版 {tDesignNew,8:F3} ms");
            Console.WriteLine($"  处理:  现实现 {tOld,8:F2} ms   优化版 {tNew,8:F2} ms   加速 {tOld / Math.Max(tNew, 0.001),6:F1}×");
            Console.WriteLine($"  分配:  现实现 {allocOld / 1024.0 / 1024.0,8:F1} MB  优化版 {allocNew / 1024.0 / 1024.0,8:F2} MB");
            Console.WriteLine($"  结果最大差异: {maxDiff:E3}");
            Console.WriteLine();
        }
    }
}
