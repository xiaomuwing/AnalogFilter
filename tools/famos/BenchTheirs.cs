using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;
using AnalogFilter;

/// <summary>
/// 基准：另一份 native 依赖实现（SciPy 原生库版），与 BenchOurs 使用完全相同的用例、信号与计时方式。
/// </summary>
internal static class BenchTheirs
{
    private sealed class Case
    {
        public int Id;
        public string Name;
        public string Kind;      // butter / cheby1 / cheby2 / bessel
        public BandType Band;
        public int Order;
        public double[] Wn;
        public double Rp;
        public double Rs;
        public BesselNorm Norm;
    }

    private static Case[] BuildCases()
    {
        return new[]
        {
            new Case { Id=0,  Name="butter4_lp",  Kind="butter", Band=BandType.Lowpass,  Order=4, Wn=new[]{0.2} },
            new Case { Id=1,  Name="butter4_hp",  Kind="butter", Band=BandType.Highpass, Order=4, Wn=new[]{0.2} },
            new Case { Id=2,  Name="butter6_bp",  Kind="butter", Band=BandType.Bandpass, Order=6, Wn=new[]{0.2,0.4} },
            new Case { Id=3,  Name="butter6_bs",  Kind="butter", Band=BandType.Bandstop, Order=6, Wn=new[]{0.2,0.4} },
            new Case { Id=4,  Name="cheby1_4_lp", Kind="cheby1", Band=BandType.Lowpass,  Order=4, Wn=new[]{0.2}, Rp=1 },
            new Case { Id=5,  Name="cheby1_4_hp", Kind="cheby1", Band=BandType.Highpass, Order=4, Wn=new[]{0.2}, Rp=1 },
            new Case { Id=6,  Name="cheby2_4_lp", Kind="cheby2", Band=BandType.Lowpass,  Order=4, Wn=new[]{0.2}, Rs=40 },
            new Case { Id=7,  Name="cheby2_4_bp", Kind="cheby2", Band=BandType.Bandpass, Order=4, Wn=new[]{0.2,0.4}, Rs=40 },
            new Case { Id=8,  Name="bessel4_lp_phase", Kind="bessel", Band=BandType.Lowpass, Order=4, Wn=new[]{0.2}, Norm=BesselNorm.Phase },
            new Case { Id=9,  Name="bessel6_lp_mag",   Kind="bessel", Band=BandType.Lowpass, Order=6, Wn=new[]{0.2}, Norm=BesselNorm.Mag },
            new Case { Id=10, Name="bessel4_bp_mag",   Kind="bessel", Band=BandType.Bandpass, Order=4, Wn=new[]{0.2,0.4}, Norm=BesselNorm.Mag },
        };
    }

    private static double[,] Design(Case c)
    {
        switch (c.Kind)
        {
            case "butter":
                return Signal.ButterSos(c.Order, c.Wn, c.Band);
            case "cheby1":
                return Signal.Cheby1Sos(c.Order, c.Rp, c.Wn, c.Band);
            case "cheby2":
                return Signal.Cheby2Sos(c.Order, c.Rs, c.Wn, c.Band);
            default:
                return Signal.BesselSos(c.Order, c.Wn, c.Band, false, c.Norm);
        }
    }

    private static double[] ReadCsv(string path, int column, int skip)
    {
        var list = new System.Collections.Generic.List<double>();
        using (var reader = new StreamReader(path, Encoding.Default))
        {
            for (int i = 0; i < skip; i++) reader.ReadLine();
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length == 0) continue;
                string[] parts = line.Split(',');
                if (parts.Length <= column) continue;
                list.Add(double.Parse(parts[column].Trim(), CultureInfo.InvariantCulture));
            }
        }
        return list.ToArray();
    }

    private static void WriteDoubles(string path, double[] values)
    {
        using (var fs = new FileStream(path, FileMode.Create, FileAccess.Write))
        using (var bw = new BinaryWriter(fs))
            for (int i = 0; i < values.Length; i++) bw.Write(values[i]);
    }

    private static int Main(string[] args)
    {
        string csv = args[0];
        string outDir = args[1];
        Directory.CreateDirectory(outDir);

        double[] signal = ReadCsv(csv, 2, 2);
        StringBuilder timing = new StringBuilder();
        timing.AppendLine("case,design_ms_200,filter_ms_5");

        foreach (Case c in BuildCases())
        {
            double[,] sos = Design(c);
            var sw = Stopwatch.StartNew();
            for (int i = 0; i < 200; i++) sos = Design(c);
            sw.Stop();
            double designMs = sw.Elapsed.TotalMilliseconds / 200.0;

            Signal.Sosfilt(sos, signal);
            sw.Restart();
            double[] y = null;
            for (int i = 0; i < 5; i++) y = Signal.Sosfilt(sos, signal);
            sw.Stop();
            double filterMs = sw.Elapsed.TotalMilliseconds / 5.0;

            int rows = sos.GetLength(0);
            double[] flat = new double[rows * 6];
            for (int s = 0; s < rows; s++)
                for (int k = 0; k < 6; k++)
                    flat[s * 6 + k] = sos[s, k];

            WriteDoubles(Path.Combine(outDir, $"theirs_{c.Id}_sos.bin"), flat);
            WriteDoubles(Path.Combine(outDir, $"theirs_{c.Id}_y.bin"), y);
            timing.AppendLine($"{c.Name},{designMs:F6},{filterMs:F6}");
            Console.WriteLine($"[{c.Id,2}] {c.Name,-18} 设计 {designMs * 1000,8:F2} us   滤波 {filterMs,7:F2} ms");
        }

        File.WriteAllText(Path.Combine(outDir, "timing_theirs.csv"), timing.ToString());
        Console.WriteLine("输出: " + outDir);
        return 0;
    }
}
