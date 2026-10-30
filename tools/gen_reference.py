#!/usr/bin/env python3
"""生成 tests/AnalogFilter.Tests/Data/scipy_reference.txt（SciPy 对照数据）。

用法：
    python -m venv %TEMP%\\af_ref_venv
    %TEMP%\\af_ref_venv\\Scripts\\pip install numpy scipy
    %TEMP%\\af_ref_venv\\Scripts\\python tools\\gen_reference.py

本脚本是参考数据的唯一来源：测试项目只读取生成的文本文件，不依赖 Python。
基准版本见 BENCHMARK 行；数值按 IEEE-754 十七位有效数字写出，double.Parse 可无损还原。

文件格式（行式，全部用 Latin-1 安全的 ASCII）：

    CASE   <用例 id>
    SPEC   <kind> <参数...>            kind = digital | analog | besselp | acf
    W      <浮点...>                   频响取样点（FREQZ/FREQS 使用）
    SOS    <6 个浮点>                  每行一节，可重复
    POLE   <实部> <虚部>               可重复，按 (实部, 虚部) 升序
    ZERO   <实部> <虚部>               可重复，同上
    GAIN   <浮点>
    FREQZ  <实部> <虚部> ...           数字频响（与 W 一一对应）
    FREQS  <实部> <虚部> ...           模拟频响（与 W 一一对应）
    SOSFILT  <浮点...>                 sosfilt 输出
    Y_ZP_<pad>_<len>  <浮点...>        sosfiltfilt 输出，len 为 default 或具体值
    Y_FF_<pad>_<len>  <浮点...>        filtfilt 输出
    ZI_SOS <2 个浮点>                  每行一节
    ZI_LFILTER <浮点...>
    X      <浮点...>                   输入信号（acf 用例自带）
    ACF    <浮点...>
"""

import os
import sys

import numpy as np
from scipy import signal

OUT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "tests", "AnalogFilter.Tests", "Data", "scipy_reference.txt")

BTYPE = {"lp": "lowpass", "hp": "highpass", "bp": "bandpass", "bs": "bandstop"}
WN = {"lp": [0.2], "hp": [0.25], "bp": [0.1, 0.3], "bs": [0.1, 0.3]}
ANALOG_WN = {"lp": [100.0], "hp": [150.0], "bp": [100.0, 300.0], "bs": [100.0, 300.0]}


def f17(value):
    """与 C# double.Parse 往返一致的文本。"""
    value = float(value)
    if value == 0.0:
        return "0"
    if not np.isfinite(value):
        raise ValueError("参考值必须是有限数")
    return format(value, ".17g")


def fmt_row(key, values):
    return key + " " + " ".join(f17(v) for v in np.ravel(values))


def test_signal(n=50):
    """与 C# 侧 TestSignals.DesignInput 一致。"""
    t = np.arange(n)
    return np.sin(2.0 * np.pi * 0.05 * t) + 0.5 * np.cos(2.0 * np.pi * 0.31 * t)


def test_long_signal(n=200):
    """与 C# 侧 TestSignals.ZeroPhaseInput 一致。"""
    t = np.arange(n)
    return (np.sin(2.0 * np.pi * 0.02 * t)
            + 0.3 * np.cos(2.0 * np.pi * 0.17 * t)
            + 0.1 * np.sin(2.0 * np.pi * 0.41 * t))


def freq_grid(count=16, lo=0.01, hi=0.95):
    return np.linspace(lo, hi, count)


def sorted_complex(values):
    """按 (实部, 虚部) 升序，与 C# 侧排序规则一致。"""
    arr = np.asarray(values, dtype=complex)
    return sorted(arr, key=lambda v: (float(np.real(v)), float(np.imag(v))))


class Writer:
    def __init__(self):
        self.lines = []
        self.case_count = 0

    def case(self, cid, spec):
        self.case_count += 1
        self.lines.append("CASE " + cid)
        self.lines.append("SPEC " + spec)

    def row(self, key, values):
        flat = np.ravel(np.asarray(values, dtype=float))
        # 整行为空时也写键，便于 C# 侧判定“此用例含该项”
        self.lines.append(fmt_row(key, flat) if flat.size else key)

    def zeros_sos(self, sos):
        for sec in np.asarray(sos):
            self.row("SOS", sec)

    def pairs(self, key, complex_values):
        data = np.empty(len(complex_values) * 2)
        data[0::2] = np.real(complex_values)
        data[1::2] = np.imag(complex_values)
        self.row(key, data)

    def text(self, header):
        body = list(self.lines)
        return "\n".join(header + body) + "\n"


def freqz_sos(sos, w):
    try:
        _, h = signal.freqz_sos(sos, worN=w)
    except AttributeError:  # SciPy < 1.15
        _, h = signal.sosfreqz(sos, worN=w)
    return h


def scalar_wn(values):
    """SciPy 的单边带要求标量，双边带要求二元序列。"""
    return float(values[0]) if len(values) == 1 else list(values)


def design_digital(ftype, band, order, rp=1.0, rs=40.0):
    """返回 (zpk, sos, ba)。bessel 需要单独入口以传 norm。"""
    wn = scalar_wn(WN[band])
    if ftype == "bessel":
        zpk = signal.bessel(order, wn, btype=BTYPE[band], analog=False, output="zpk", norm="phase")
        sos = signal.bessel(order, wn, btype=BTYPE[band], analog=False, output="sos", norm="phase")
        ba = signal.bessel(order, wn, btype=BTYPE[band], analog=False, output="ba", norm="phase")
    else:
        zpk = signal.iirfilter(order, wn, rp=rp, rs=rs, btype=BTYPE[band],
                               analog=False, ftype=ftype, output="zpk")
        sos = signal.iirfilter(order, wn, rp=rp, rs=rs, btype=BTYPE[band],
                               analog=False, ftype=ftype, output="sos")
        ba = signal.iirfilter(order, wn, rp=rp, rs=rs, btype=BTYPE[band],
                              analog=False, ftype=ftype, output="ba")
    return zpk, sos, ba


def design_analog(ftype, band, order, rp=1.0, rs=40.0):
    wn = scalar_wn(ANALOG_WN[band])
    if ftype == "bessel":
        return signal.bessel(order, wn, btype=BTYPE[band], analog=True, output="zpk", norm="phase")
    return signal.iirfilter(order, wn, rp=rp, rs=rs, btype=BTYPE[band],
                            analog=True, ftype=ftype, output="zpk")


def emit_digital(w, cases, zero_phase=()):
    for (ftype, band, order) in cases:
        cid = "digital_%s_%s_n%d" % (ftype, band, order)
        _, sos, ba = design_digital(ftype, band, order)
        wgrid = freq_grid()
        x = test_signal()
        w.case(cid, "digital %s %s %d %s %s %s" % (
            ftype, band, order, " ".join(f17(v) for v in WN[band]), f17(1.0), f17(40.0)))
        w.row("W", wgrid)
        w.zeros_sos(sos)
        w.pairs("FREQZ", freqz_sos(sos, wgrid))
        w.row("SOSFILT", signal.sosfilt(sos, x))
        w.row("ZI_SOS", np.ravel(signal.sosfilt_zi(sos)))
        w.row("ZI_LFILTER", signal.lfilter_zi(ba[0], ba[1]))

        if cid in zero_phase:
            xz = test_long_signal()
            for padtype in ("odd", "even", "constant"):
                w.row("Y_ZP_%s_default" % padtype,
                      signal.sosfiltfilt(sos, xz, padtype=padtype))
                w.row("Y_ZP_%s_20" % padtype,
                      signal.sosfiltfilt(sos, xz, padtype=padtype, padlen=20))
            w.row("Y_FF_odd_20", signal.filtfilt(ba[0], ba[1], xz, padtype="odd", padlen=20))
            w.row("Y_FF_even_default", signal.filtfilt(ba[0], ba[1], xz, padtype="even"))


def emit_analog(cases):
    for (ftype, band, order) in cases:
        cid = "analog_%s_%s_n%d" % (ftype, band, order)
        z, p, k = design_analog(ftype, band, order)
        wgrid = np.linspace(10.0, 400.0, 12)
        w.case(cid, "analog %s %s %d %s %s %s" % (
            ftype, band, order, " ".join(f17(v) for v in ANALOG_WN[band]), f17(1.0), f17(40.0)))
        w.row("W", wgrid)
        w.pairs("ZERO", sorted_complex(z))
        w.pairs("POLE", sorted_complex(p))
        w.row("GAIN", [k])
        _, h = signal.freqs_zpk(z, p, k, worN=wgrid)
        w.pairs("FREQS", h)


def emit_bessel(w):
    """贝塞尔原型极点精度随阶数下降，单独固化以便回归。"""
    for order in range(2, 13):
        z, p, k = signal.besselap(order, norm="phase")
        w.case("besselp_n%d" % order, "besselp %d" % order)
        w.pairs("POLE", sorted_complex(p))
        w.row("GAIN", [k])


def autocorrelation_reference(samples, length):
    """循环自相关，与 Spectrum.Autocorrelation 的定义一致。"""
    x = np.asarray(samples[:length], dtype=float)
    spectrum = np.fft.fft(x)
    power = np.abs(spectrum) ** 2
    r = np.real(np.fft.ifft(power))
    if r[0] == 0.0:
        return np.zeros(length)
    return r / r[0]


def emit_acf(w):
    rng = np.random.default_rng(20261023)
    datasets = {
        "impulse_64": np.eye(1, 64, 0).ravel(),
        "constant_128": np.ones(128),
        "sine_512": np.sin(2.0 * np.pi * 3.0 * np.arange(512) / 512.0),
        "noise_1024": rng.standard_normal(1024),
        "chirp_256": np.sin(2.0 * np.pi * np.arange(256) ** 2 / 8192.0),
    }
    for name, data in datasets.items():
        n = len(data)
        w.case("acf_" + name, "acf %d" % n)
        w.row("X", data)
        w.row("ACF", autocorrelation_reference(data, n))


def main():
    global w
    w = Writer()

    digital_cases = [
        (ftype, band, order)
        for ftype in ("butter", "cheby1", "cheby2", "bessel")
        for band in ("lp", "hp", "bp", "bs")
        for order in (2, 4)
    ]
    digital_cases += [
        ("butter", "lp", 3), ("butter", "lp", 6), ("butter", "hp", 3),
        ("cheby1", "bp", 6), ("cheby2", "bs", 6), ("bessel", "bp", 4),
    ]
    zero_phase = {
        "digital_butter_lp_n4", "digital_cheby1_bp_n4",
        "digital_bessel_lp_n4", "digital_butter_hp_n4",
    }
    emit_digital(w, digital_cases, zero_phase)

    # 高阶数字设计：SOS 路径的意义就在于这些阶数
    for order in (10, 16):
        _, sos, _ = design_digital("butter", "lp", order)
        wgrid = freq_grid()
        w.case("highorder_butter_lp_n%d" % order,
               "digital butter lp %d %s %s %s" % (order, f17(0.2), f17(1.0), f17(40.0)))
        w.row("W", wgrid)
        w.zeros_sos(sos)
        w.pairs("FREQZ", freqz_sos(sos, wgrid))
        w.row("SOSFILT", signal.sosfilt(sos, test_signal()))

    _, sos, _ = design_digital("cheby1", "lp", 8)
    wgrid = freq_grid()
    w.case("highorder_cheby1_lp_n8", "digital cheby1 lp 8 %s %s %s" % (f17(0.2), f17(1.0), f17(40.0)))
    w.row("W", wgrid)
    w.zeros_sos(sos)
    w.pairs("FREQZ", freqz_sos(sos, wgrid))
    w.row("SOSFILT", signal.sosfilt(sos, test_signal()))

    analog_cases = [
        (ftype, band, order)
        for ftype in ("butter", "cheby1", "cheby2", "bessel")
        for band in ("lp", "hp", "bp", "bs")
        for order in (2, 4)
    ]
    emit_analog(analog_cases)
    emit_bessel(w)
    emit_acf(w)

    header = [
        "# AnalogFilter 的 SciPy 对照数据（自动生成，请勿手工编辑）",
        "# 生成：tools/gen_reference.py",
        "# 基准：scipy %s / numpy %s" % (__import__("scipy").__version__, np.__version__),
        "# 用例数：%d" % w.case_count,
    ]
    text = w.text(header)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)

    print("已写入 %s" % OUT_PATH)
    print("用例 %d 个，%.1f KB" % (w.case_count, len(text) / 1024.0))


if __name__ == "__main__":
    sys.exit(main())
