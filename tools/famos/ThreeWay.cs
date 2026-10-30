using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using AnalogFilter;

/// <summary>
/// 三方对比数据生成器：读取 11.csv，对同一组滤波参数分别用 FAMOS COM 与 AnalogFilter 计算，
/// 把输入信号与两路输出写成二进制文件，供 Python 侧与 SciPy 结果合并比较。
/// </summary>
internal static class ThreeWay
{
    private sealed class Config
    {
        public int Id;
        public string Band;          // LowPass / HighPass / BandPass / BandStop
        public int Characteristic;   // FAMOS: 0=Butter 1=Bessel 2=Chebyshev 3=CritDamped
        public double Parameter;
        public int Order;
        public double Freq1;
        public double Freq2;
        public string Name;
    }

    private static double[] ReadCsvChannel(string path, int columnIndex, int skipLines)
    {
        List<double> values = new List<double>();
        using (StreamReader reader = new StreamReader(path, Encoding.Default))
        {
            for (int i = 0; i < skipLines; i++)
            {
                reader.ReadLine();
            }

            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length == 0)
                {
                    continue;
                }

                string[] parts = line.Split(',');
                if (parts.Length <= columnIndex)
                {
                    continue;
                }

                values.Add(double.Parse(parts[columnIndex].Trim(), CultureInfo.InvariantCulture));
            }
        }
        return values.ToArray();
    }

    private static void WriteDoubles(string path, double[] values)
    {
        using (FileStream fs = new FileStream(path, FileMode.Create, FileAccess.Write))
        using (BinaryWriter bw = new BinaryWriter(fs))
        {
            for (int i = 0; i < values.Length; i++)
            {
                bw.Write(values[i]);
            }
        }
    }

    private static int Main(string[] args)
    {
        string csvPath = args[0];
        string outDir = args[1];
        double fs = 1000.0;

        Directory.CreateDirectory(outDir);

        double[] signal = ReadCsvChannel(csvPath, 2, 2);
        Console.WriteLine($"输入通道 Strain1: n={signal.Length}  fs={fs} Hz");
        WriteDoubles(Path.Combine(outDir, "signal.bin"), signal);

        List<Config> configs = new List<Config>();
        int id = 0;
        int[] cases = { -1 }; // 统一在下面构造

        // 三种 SciPy 有对应原型的特征 × 四种频带，阶数 4
        string[] bands = { "LowPass", "HighPass", "BandPass", "BandStop" };
        for (int characteristic = 0; characteristic <= 2; characteristic++)
        {
            foreach (string band in bands)
            {
                bool twoFreq = band == "BandPass" || band == "BandStop";
                configs.Add(new Config
                {
                    Id = id++,
                    Band = band,
                    Characteristic = characteristic,
                    Parameter = characteristic == 2 ? 1.0 : 0.0, // 切比雪夫纹波 1 dB
                    Order = 4,
                    Freq1 = twoFreq ? 20.0 : 50.0,
                    Freq2 = twoFreq ? 80.0 : 0.0,
                    Name = characteristic switch { 0 => "butter", 1 => "bessel", _ => "cheby1" }
                });
            }
        }

        // 临界阻尼：FAMOS 专有，SciPy/AnalogFilter 无对应原型，仅记录 FAMOS 输出
        configs.Add(new Config
        {
            Id = id++,
            Band = "LowPass",
            Characteristic = 3,
            Parameter = 0.0,
            Order = 4,
            Freq1 = 50.0,
            Freq2 = 0.0,
            Name = "critdamped"
        });

        // 阶数 2 的巴特沃斯低通，用于观察阶数影响
        configs.Add(new Config
        {
            Id = id++,
            Band = "LowPass",
            Characteristic = 0,
            Parameter = 0.0,
            Order = 2,
            Freq1 = 50.0,
            Freq2 = 0.0,
            Name = "butter"
        });

        var famos = new ImcFamosLib.Famos();

        StringBuilder manifest = new StringBuilder();
        foreach (Config cfg in configs)
        {
            var ch = new ImcCoreLib.DChannel();
            ch.Data = (double[])signal.Clone();
            ch.xDelta = 1.0 / fs;

            double[] famosOut = null;
            try
            {
                ImcCoreLib.DChannel result;
                switch (cfg.Band)
                {
                    case "LowPass":
                        result = famos.Filter.LowPass(ch, cfg.Characteristic, cfg.Parameter, cfg.Order, cfg.Freq1);
                        break;
                    case "HighPass":
                        result = famos.Filter.HighPass(ch, cfg.Characteristic, cfg.Parameter, cfg.Order, cfg.Freq1);
                        break;
                    case "BandPass":
                        result = famos.Filter.BandPass(ch, cfg.Characteristic, cfg.Parameter, cfg.Order, cfg.Freq1, cfg.Freq2);
                        break;
                    default:
                        result = famos.Filter.BandStop(ch, cfg.Characteristic, cfg.Parameter, cfg.Order, cfg.Freq1, cfg.Freq2);
                        break;
                }
                famosOut = ToArray(result.Data);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"  [{cfg.Id}] FAMOS 失败: {ex.Message}");
            }

            double[] analogOut = null;
            if (cfg.Characteristic <= 2)
            {
                try
                {
                    FilterType type = cfg.Characteristic switch
                    {
                        0 => FilterType.Butterworth,
                        1 => FilterType.Bessel,
                        _ => FilterType.ChebyshevI
                    };
                    BandType band = cfg.Band switch
                    {
                        "LowPass" => BandType.LowPass,
                        "HighPass" => BandType.HighPass,
                        "BandPass" => BandType.BandPass,
                        _ => BandType.BandStop
                    };
                    double[] wn = band == BandType.BandPass || band == BandType.BandStop
                        ? new[] { cfg.Freq1 / (fs / 2.0), cfg.Freq2 / (fs / 2.0) }
                        : new[] { cfg.Freq1 / (fs / 2.0) };

                    double[][] sos = FilterDesign.DesignSos(
                        cfg.Order, wn, type, band, cfg.Parameter == 0.0 ? 1.0 : cfg.Parameter, 40.0);
                    analogOut = SosFiltering.SosFilt(sos, signal);
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"  [{cfg.Id}] AnalogFilter 失败: {ex.Message}");
                }
            }

            if (famosOut != null) WriteDoubles(Path.Combine(outDir, $"famos_{cfg.Id}.bin"), famosOut);
            if (analogOut != null) WriteDoubles(Path.Combine(outDir, $"analog_{cfg.Id}.bin"), analogOut);

            manifest.AppendLine(string.Join("|",
                cfg.Id, cfg.Name, cfg.Band, cfg.Characteristic, cfg.Parameter, cfg.Order,
                cfg.Freq1.ToString("R", CultureInfo.InvariantCulture),
                cfg.Freq2.ToString("R", CultureInfo.InvariantCulture),
                famosOut == null ? "fail" : "ok",
                analogOut == null ? "none" : "ok"));

            Console.WriteLine($"[{cfg.Id,2}] {cfg.Name,-11} {cfg.Band,-9} N={cfg.Order} " +
                              $"f={cfg.Freq1}{(cfg.Freq2 > 0 ? "/" + cfg.Freq2 : "")}Hz  " +
                              $"FAMOS={(famosOut == null ? "fail" : "ok")}  Analog={(analogOut == null ? "none" : "ok")}");
        }

        File.WriteAllText(Path.Combine(outDir, "manifest.txt"), manifest.ToString());
        Console.WriteLine();
        Console.WriteLine($"输出目录: {outDir}");
        return 0;
    }

    private static double[] ToArray(object comData)
    {
        Array array = (Array)comData;
        int lower = array.GetLowerBound(0);
        double[] result = new double[array.Length];
        for (int i = 0; i < array.Length; i++)
        {
            result[i] = Convert.ToDouble(array.GetValue(lower + i));
        }
        return result;
    }
}
