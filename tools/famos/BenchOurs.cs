using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;
using AnalogFilter;

/// <summary>
/// 基准：本仓库的托管 AnalogFilter（src\AnalogFilter）。
/// 输出 SOS 系数、滤波结果和耗时，供与 AnalogFilter_CC 对比。
/// </summary>
internal static class BenchOurs
{
    private sealed class Case
    {
        public int Id;
        public string Name;
        public FilterType Type;
        public BandType Band;
        public int Order;
        public double[] Wn;
        public double Rp;
        public double Rs;
        public BesselNormalization Norm;
    }

    private static Case[] BuildCases()
    {
        return new[]
        {
            new Case { Id=0,  Name="butter4_lp",  Type=FilterType.Butterworth, Band=BandType.LowPass,  Order=4, Wn=new[]{0.2} },
            new Case { Id=1,  Name="butter4_hp",  Type=FilterType.Butterworth, Band=BandType.HighPass, Order=4, Wn=new[]{0.2} },
            new Case { Id=2,  Name="butter6_bp",  Type=FilterType.Butterworth, Band=BandType.BandPass, Order=6, Wn=new[]{0.2,0.4} },
            new Case { Id=3,  Name="butter6_bs",  Type=FilterType.Butterworth, Band=BandType.BandStop, Order=6, Wn=new[]{0.2,0.4} },
            new Case { Id=4,  Name="cheby1_4_lp", Type=FilterType.ChebyshevI,  Band=BandType.LowPass,  Order=4, Wn=new[]{0.2}, Rp=1 },
            new Case { Id=5,  Name="cheby1_4_hp", Type=FilterType.ChebyshevI,  Band=BandType.HighPass, Order=4, Wn=new[]{0.2}, Rp=1 },
            new Case { Id=6,  Name="cheby2_4_lp", Type=FilterType.ChebyshevII, Band=BandType.LowPass,  Order=4, Wn=new[]{0.2}, Rs=40 },
            new Case { Id=7,  Name="cheby2_4_bp", Type=FilterType.ChebyshevII, Band=BandType.BandPass, Order=4, Wn=new[]{0.2,0.4}, Rs=40 },
            new Case { Id=8,  Name="bessel4_lp_phase", Type=FilterType.Bessel, Band=BandType.LowPass,  Order=4, Wn=new[]{0.2}, Norm=BesselNormalization.Phase },
            new Case { Id=9,  Name="bessel6_lp_mag",   Type=FilterType.Bessel, Band=BandType.LowPass,  Order=6, Wn=new[]{0.2}, Norm=BesselNormalization.Magnitude },
            new Case { Id=10, Name="bessel4_bp_mag",   Type=FilterType.Bessel, Band=BandType.BandPass, Order=4, Wn=new[]{0.2,0.4}, Norm=BesselNormalization.Magnitude },
        };
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
            // 设计（含预热，取 200 次平均）
            double[][] sos = FilterDesign.DesignSos(c.Order, c.Wn, c.Type, c.Band, c.Rp, c.Rs, false, null, c.Norm);
            var sw = Stopwatch.StartNew();
            for (int i = 0; i < 200; i++)
                sos = FilterDesign.DesignSos(c.Order, c.Wn, c.Type, c.Band, c.Rp, c.Rs, false, null, c.Norm);
            sw.Stop();
            double designMs = sw.Elapsed.TotalMilliseconds / 200.0;

            // 滤波
            SosFiltering.SosFilt(sos, signal);
            sw.Restart();
            double[] y = null;
            for (int i = 0; i < 5; i++) y = SosFiltering.SosFilt(sos, signal);
            sw.Stop();
            double filterMs = sw.Elapsed.TotalMilliseconds / 5.0;

            double[] flat = new double[sos.Length * 6];
            for (int s = 0; s < sos.Length; s++)
                for (int k = 0; k < 6; k++)
                    flat[s * 6 + k] = sos[s][k];

            WriteDoubles(Path.Combine(outDir, $"ours_{c.Id}_sos.bin"), flat);
            WriteDoubles(Path.Combine(outDir, $"ours_{c.Id}_y.bin"), y);
            timing.AppendLine($"{c.Name},{designMs:F6},{filterMs:F6}");
            Console.WriteLine($"[{c.Id,2}] {c.Name,-18} 设计 {designMs * 1000,8:F2} us   滤波 {filterMs,7:F2} ms");
        }

        File.WriteAllText(Path.Combine(outDir, "timing_ours.csv"), timing.ToString());
        Console.WriteLine("输出: " + outDir);
        return 0;
    }
}
