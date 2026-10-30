using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// <c>Data\scipy_reference.txt</c> 中的一个用例（由 <c>tools\gen_reference.py</c> 生成）。
    /// </summary>
    internal sealed class ReferenceCase
    {
        internal ReferenceCase(string id, string[] spec)
        {
            Id = id;
            Spec = spec;
        }

        /// <summary>用例 id，例如 <c>digital_butter_lp_n4</c>。</summary>
        internal string Id { get; private set; }

        /// <summary><c>SPEC</c> 行的标记（不含 SPEC 本身）。</summary>
        internal string[] Spec { get; private set; }

        private readonly Dictionary<string, List<double[]>> _rows =
            new Dictionary<string, List<double[]>>(StringComparer.Ordinal);

        internal void Add(string key, double[] values)
        {
            List<double[]> list;
            if (!_rows.TryGetValue(key, out list))
            {
                list = new List<double[]>();
                _rows.Add(key, list);
            }
            list.Add(values);
        }

        internal bool Has(string key)
        {
            return _rows.ContainsKey(key);
        }

        /// <summary>取单行数据（缺项返回 null）。</summary>
        internal double[] Single(string key)
        {
            List<double[]> list;
            return _rows.TryGetValue(key, out list) ? list[0] : null;
        }

        /// <summary>取多行数据（如 SOS 每节一行）。</summary>
        internal List<double[]> Multi(string key)
        {
            List<double[]> list;
            return _rows.TryGetValue(key, out list) ? list : new List<double[]>();
        }
    }

    /// <summary>SciPy 对照数据的解析入口。</summary>
    internal static class ReferenceData
    {
        private static readonly Lazy<List<ReferenceCase>> Lazy = new Lazy<List<ReferenceCase>>(Load);

        internal static IReadOnlyList<ReferenceCase> All
        {
            get { return Lazy.Value; }
        }

        internal static ReferenceCase Get(string id)
        {
            foreach (ReferenceCase item in Lazy.Value)
            {
                if (item.Id == id)
                {
                    return item;
                }
            }
            throw new InvalidOperationException("对照数据中没有用例：" + id);
        }

        private static List<ReferenceCase> Load()
        {
            string path = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Data", "scipy_reference.txt");
            if (!File.Exists(path))
            {
                throw new FileNotFoundException(
                    "缺少 SciPy 对照数据。请运行 tools\\gen_reference.py 生成后重新构建测试项目。", path);
            }

            List<ReferenceCase> cases = new List<ReferenceCase>();
            ReferenceCase current = null;

            foreach (string raw in File.ReadAllLines(path))
            {
                string line = raw.Trim();
                if (line.Length == 0 || line[0] == '#')
                {
                    continue;
                }

                int space = line.IndexOf(' ');
                string key = space < 0 ? line : line.Substring(0, space);
                string rest = space < 0 ? string.Empty : line.Substring(space + 1).Trim();

                if (key == "CASE")
                {
                    current = new ReferenceCase(rest, null);
                    cases.Add(current);
                    continue;
                }
                if (current == null)
                {
                    throw new InvalidDataException("对照数据格式异常：第一个非注释行必须是 CASE。");
                }
                if (key == "SPEC")
                {
                    current = new ReferenceCase(
                        current.Id,
                        rest.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries));
                    cases[cases.Count - 1] = current;
                    continue;
                }

                current.Add(key, ParseDoubles(rest));
            }

            return cases;
        }

        private static double[] ParseDoubles(string text)
        {
            string[] parts = text.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            double[] values = new double[parts.Length];
            for (int i = 0; i < parts.Length; i++)
            {
                values[i] = double.Parse(parts[i], CultureInfo.InvariantCulture);
            }
            return values;
        }
    }
}
