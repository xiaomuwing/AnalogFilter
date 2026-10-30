using System;

namespace AnalogFilter.Tests
{
    /// <summary>
    /// 测试输入信号。公式必须与 <c>tools\gen_reference.py</c> 中的同名函数逐项一致，
    /// 否则 SciPy 参照值与本库输出没有可比性。
    /// </summary>
    internal static class TestSignals
    {
        /// <summary>设计用例的输入：50 点，两个正弦叠加。对应 Python 的 <c>test_signal()</c>。</summary>
        internal static double[] DesignInput()
        {
            double[] x = new double[50];
            for (int n = 0; n < x.Length; n++)
            {
                x[n] = Math.Sin(2.0 * Math.PI * 0.05 * n) + 0.5 * Math.Cos(2.0 * Math.PI * 0.31 * n);
            }
            return x;
        }

        /// <summary>零相位滤波用例的输入：200 点，三路叠加。对应 Python 的 <c>test_long_signal()</c>。</summary>
        internal static double[] ZeroPhaseInput()
        {
            double[] x = new double[200];
            for (int n = 0; n < x.Length; n++)
            {
                x[n] = Math.Sin(2.0 * Math.PI * 0.02 * n)
                     + 0.3 * Math.Cos(2.0 * Math.PI * 0.17 * n)
                     + 0.1 * Math.Sin(2.0 * Math.PI * 0.41 * n);
            }
            return x;
        }

        /// <summary>纯正弦，<paramref name="cyclesPerSample"/> 单位是 cycles/sample（不是 rad/sample）。</summary>
        internal static double[] Sine(double cyclesPerSample, int length)
        {
            double[] x = new double[length];
            for (int n = 0; n < length; n++)
            {
                x[n] = Math.Sin(2.0 * Math.PI * cyclesPerSample * n);
            }
            return x;
        }

        /// <summary>阶跃信号，用于验证临界阻尼的单调性。</summary>
        internal static double[] Step(int length)
        {
            double[] x = new double[length];
            for (int n = 0; n < length; n++)
            {
                x[n] = 1.0;
            }
            return x;
        }
    }
}
