import sys, os, struct
import numpy as np
from scipy.signal import butter, cheby1, bessel, sosfilt

outdir = sys.argv[1]

def read_bin(path):
    if not os.path.exists(path):
        return None
    with open(path, 'rb') as f:
        raw = f.read()
    return np.frombuffer(raw, dtype='<f8').copy()

signal = read_bin(os.path.join(outdir, 'signal.bin'))
n = signal.size
fs = 1000.0
pp = signal.max() - signal.min()
rms = float(np.sqrt(np.mean(signal**2)))
EDGE = 2000
interior = slice(EDGE, n - EDGE)

print("输入: n=%d  fs=%.0f Hz  峰峰值=%.6g  RMS=%.6g" % (n, fs, pp, rms))
print()

def metrics(a, b):
    d = a - b
    di = d[interior]
    return (np.max(np.abs(di)) / pp * 100.0,
            float(np.sqrt(np.mean(di**2))) / rms * 100.0)

def design(name, band, order, f1, f2, param):
    wn = [f1/(fs/2), f2/(fs/2)] if band in ('BandPass', 'BandStop') else f1/(fs/2)
    bt = {'LowPass': 'lowpass', 'HighPass': 'highpass',
          'BandPass': 'bandpass', 'BandStop': 'bandstop'}[band]
    if name == 'butter':
        return butter(order, wn, btype=bt, output='sos')
    if name == 'cheby1':
        return cheby1(order, param if param > 0 else 1.0, wn, btype=bt, output='sos')
    if name == 'bessel':
        return bessel(order, wn, btype=bt, output='sos', norm='phase')
    raise ValueError(name)

manifest = [l.strip().split('|') for l in open(os.path.join(outdir, 'manifest.txt'), encoding='utf-8') if l.strip()]

rows = []
for spec in manifest:
    cid, name, band, char, param, order, f1, f2, famos_ok, analog_ok = spec
    cid = int(cid); param = float(param); order = int(order); f1 = float(f1); f2 = float(f2)
    famos = read_bin(os.path.join(outdir, 'famos_%d.bin' % cid))
    analog = read_bin(os.path.join(outdir, 'analog_%d.bin' % cid))

    if name == 'critdamped':
        rows.append((cid, name, band, order, f1, f2, None, None, None, 'FAMOS 专有，无对应实现'))
        continue

    sos = design(name, band, order, f1, f2, param)
    scipy = sosfilt(sos, signal)

    a_vs_s = metrics(analog, scipy)
    f_vs_s = metrics(famos, scipy)
    a_vs_f = metrics(analog, famos)

    rows.append((cid, name, band, order, f1, f2, a_vs_s, f_vs_s, a_vs_f, ''))

print("=" * 122)
print("%-3s %-11s %-9s %-3s %-9s | %-19s | %-19s | %-19s" %
      ("ID", "类型", "频带", "N", "截止Hz", "AnalogFilter vs SciPy", "FAMOS vs SciPy", "AnalogFilter vs FAMOS"))
print("=" * 122)
for cid, name, band, order, f1, f2, a_s, f_s, a_f, note in rows:
    ff = "%.0f/%.0f" % (f1, f2) if f2 > 0 else "%.0f" % f1
    if a_s is None:
        print("%-3d %-11s %-9s %-3d %-9s | %-19s | %-19s | %-19s  %s" %
              (cid, name, band, order, ff, "-", "-", "-", note))
        continue
    def fmt(p):
        return "%8.3f%% / %7.3f%%" % (p[0], p[1])
    print("%-3d %-11s %-9s %-3d %-9s | %-19s | %-19s | %-19s" %
          (cid, name, band, order, ff, fmt(a_s), fmt(f_s), fmt(a_f)))

print("=" * 122)
print("说明: 数值为 最大绝对偏差 / 峰峰值 %%   与   RMS 偏差 / 信号RMS %%（已跳过两端各 2000 点，避开启动瞬态）")
print()

# ---- 贝塞尔归一化归属判定 ----
print("=== 贝塞尔：FAMOS 到底对应哪种归一化？（与 FAMOS 输出比较，取最接近者）===")
for spec in manifest:
    cid, name, band, char, param, order, f1, f2, _, _ = spec
    if name != 'bessel':
        continue
    cid = int(cid); order = int(order); f1 = float(f1); f2 = float(f2)
    famos = read_bin(os.path.join(outdir, 'famos_%d.bin' % cid))
    wn = [f1/(fs/2), f2/(fs/2)] if band in ('BandPass', 'BandStop') else f1/(fs/2)
    bt = {'LowPass': 'lowpass', 'HighPass': 'highpass',
          'BandPass': 'bandpass', 'BandStop': 'bandstop'}[band]
    scored = []
    for norm in ('phase', 'delay', 'mag'):
        y = sosfilt(bessel(order, wn, btype=bt, output='sos', norm=norm), signal)
        scored.append((metrics(famos, y)[0], norm))
    scored.sort()
    ff = "%.0f/%.0f" % (f1, f2) if f2 > 0 else "%.0f" % f1
    detail = "   ".join("%s=%.3f%%" % (nm, v) for v, nm in scored)
    print("  [%2d] %-9s %-3d阶 %-9s  最接近: %-6s   (%s)" % (cid, band, order, ff, scored[0][1], detail))
