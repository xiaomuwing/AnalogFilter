using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using AnalogFilter;

/// <summary>
/// 临界阻尼的 -3dB 截止换算因子扫描：对每种频带，找出使 FAMOS 与本地实现最吻合的缩放因子。
/// 直接调用公开的 ZPK 原语构造滤波器，绕开 FilterDesign 内部的固定因子。
/// </summary>
internal static class CritSweep
{
    private static double[] ReadCsvChannel(string path, int columnIndex, int skipLines)
    {
        List<double> values = new List<double>();
        using (StreamReader reader = new StreamReader(path, Encoding.Default))
        {
            for (int i = 0; i < skipLines; i++) reader.ReadLine();
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
            result[i] = Convert.ToDouble(array.GetValue(lower + i));
        return result;
    }

    private static double[] Build(double[] signal, int order, BandType band, double w1, double w2, double mult)
    {
        ZpkModel proto = AnalogPrototypes.CriticallyDamped(order);
        // 双线性预畸变
        double p1 = 4.0 * Math.Tan(Math.PI * w1 / 2.0);
        double p2 = w2 > 0 ? 4.0 * Math.Tan(Math.PI * w2 / 2.0) : 0.0;
        ZpkModel analog;
        switch (band)
        {
            case BandType.LowPass:
                analog = AnalogTransforms.LowPassToLowPass(proto, p1 * mult);
                break;
            case BandType.HighPass:
                analog = AnalogTransforms.LowPassToHighPass(proto, p1 * mult);
                break;
            case BandType.BandPass:
                analog = AnalogTransforms.LowPassToBandPass(
                    proto, Math.Sqrt(p1 * p2) * mult, (p2 - p1) * mult);
                break;
            default:
                analog = AnalogTransforms.LowPassToBandStop(
                    proto, Math.Sqrt(p1 * p2) * mult, (p2 - p1) * mult);
                break;
        }
        ZpkModel digital = BilinearTransform.Zpk(analog, 2.0);
        double[][] sos = SecondOrderSections.Zpk2Sos(digital, false);
        return SosFiltering.SosFilt(sos, signal);
    }

    private static double Score(double[] reference, double[] candidate, double pp, double rms, int edge)
    {
        int start = edge;
        int stop = reference.Length - edge;
        double sumSq = 0.0;
        int count = 0;
        for (int i = start; i < stop; i++)
        {
            double d = reference[i] - candidate[i];
            sumSq += d * d;
            count++;
        }
        return Math.Sqrt(sumSq / count) / rms * 100.0;
    }

    private static int Main(string[] args)
    {
        string csvPath = args[0];
        double fs = 1000.0;
        double[] signal = ReadCsvChannel(csvPath, 2, 2);

        double min = double.MaxValue, max = double.MinValue, sumSq = 0.0;
        foreach (double v in signal)
        {
            min = Math.Min(min, v);
            max = Math.Max(max, v);
            sumSq += v * v;
        }
        double pp = max - min;
        double rms = Math.Sqrt(sumSq / signal.Length);

        var famos = new ImcFamosLib.Famos();

        (string Label, BandType Band, int PrototypeOrder, int FamosOrder, double F1, double F2)[] cases =
        {
            ("低通", BandType.LowPass, 4, 4, 50.0, 0.0),
            ("高通", BandType.HighPass, 4, 4, 50.0, 0.0),
            ("带通", BandType.BandPass, 2, 4, 20.0, 80.0),
            ("带阻", BandType.BandStop, 2, 4, 20.0, 80.0),
        };

        Console.WriteLine("=== 临界阻尼：-3dB 换算因子扫描（因子 1.0 = 无补偿）===");
        Console.WriteLine("输入 n=" + signal.Length + "  峰峰值=" + pp.ToString("F3") + "  RMS=" + rms.ToString("F3"));
        Console.WriteLine();

        foreach (var testCase in cases)
        {
            var ch = new ImcCoreLib.DChannel();
            ch.Data = (double[])signal.Clone();
            ch.xDelta = 1.0 / fs;

            ImcCoreLib.DChannel famosResult = testCase.Band switch
            {
                BandType.LowPass => famos.Filter.LowPass(ch, 3, 0.0, testCase.FamosOrder, testCase.F1),
                BandType.HighPass => famos.Filter.HighPass(ch, 3, 0.0, testCase.FamosOrder, testCase.F1),
                BandType.BandPass => famos.Filter.BandPass(ch, 3, 0.0, testCase.FamosOrder, testCase.F1, testCase.F2),
                _ => famos.Filter.BandStop(ch, 3, 0.0, testCase.FamosOrder, testCase.F1, testCase.F2)
            };
            double[] reference = ToArray(famosResult.Data);

            double w1 = testCase.F1 / (fs / 2.0);
            double w2 = testCase.F2 / (fs / 2.0);
            double s = Math.Sqrt(Math.Pow(2.0, 1.0 / testCase.PrototypeOrder) - 1.0);

            double bestMult = 0.0;
            double bestScore = double.MaxValue;
            for (int i = 0; i <= 400; i++)
            {
                double mult = 0.2 * Math.Pow(25.0, i / 400.0); // 0.2 → 5.0
                double score = Score(reference, Build(signal, testCase.PrototypeOrder, testCase.Band, w1, w2, mult), pp, rms, 2000);
                if (score < bestScore)
                {
                    bestScore = score;
                    bestMult = mult;
                }
            }

            Console.WriteLine($"[{testCase.Label}] 原型阶数 N={testCase.PrototypeOrder}  sqrt(2^(1/N)-1)={s:F6}");
            Console.WriteLine($"    最优因子 = {bestMult:F6}   rms误差 = {bestScore:F6}%");
            Console.WriteLine($"    参考：1/s={1.0 / s:F6}   s={s:F6}   1.0");
            foreach (double probe in new[] { 1.0 / s, s, 1.0, 1.0 / Math.Sqrt(s), Math.Sqrt(s) })
            {
                double sc = Score(reference, Build(signal, testCase.PrototypeOrder, testCase.Band, w1, w2, probe), pp, rms, 2000);
                Console.WriteLine($"      因子 {probe,9:F6} -> rms={sc,9:F4}%");
            }
            Console.WriteLine();
        }

        return 0;
    }
}

