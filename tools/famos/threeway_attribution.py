import sys, os
import numpy as np
from scipy.signal import butter, cheby1, cheby2, bessel, sosfilt, buttap, cheb1ap, cheb2ap, besselap
from scipy.signal import lp2bp_zpk, lp2bs_zpk, bilinear_zpk, zpk2sos

outdir = sys.argv[1]
def read_bin(p):
    if not os.path.exists(p): return None
    return np.frombuffer(open(p,'rb').read(), dtype='<f8').copy()

signal = read_bin(os.path.join(outdir,'signal.bin'))
n = signal.size; fs = 1000.0
pp = signal.max()-signal.min(); rms = float(np.sqrt(np.mean(signal**2)))
EDGE = 2000; interior = slice(EDGE, n-EDGE)

def score(a, b):
    d = (a-b)[interior]
    return np.max(np.abs(d))/pp*100.0, float(np.sqrt(np.mean(d**2)))/rms*100.0

def show(label, famos, y):
    m = score(famos, y)
    print("    %-46s max=%7.3f%%  rms=%7.3f%%" % (label, m[0], m[1]))
    return m[0]

print("=== A. 切比雪夫低通 N=4 50Hz：FAMOS 的 parameter 到底是什么？===")
famos = read_bin(os.path.join(outdir,'famos_8.bin'))
best=(1e9,None)
for rp in [0.1,0.2,0.5,1.0,2.0,3.0]:
    show("cheby1 rp=%.1f dB" % rp, famos, sosfilt(cheby1(4,rp,0.1,output='sos'), signal))
for rs in [3,6,10,20,40]:
    show("cheby2 rs=%.0f dB" % rs, famos, sosfilt(cheby2(4,rs,0.1,output='sos'), signal))
print("    参考: butter (rp=0 的极限)")
show("butter", famos, sosfilt(butter(4,0.1,output='sos'), signal))

print()
print("=== B. 巴特沃斯带通 20/80Hz：FAMOS 的 order 与中心频率约定？===")
famos_bp = read_bin(os.path.join(outdir,'famos_2.bin'))
famos_bs = read_bin(os.path.join(outdir,'famos_3.bin'))
w1, w2 = 20.0, 80.0
nw1, nw2 = w1/(fs/2), w2/(fs/2)
show("scipy butter(4, bandpass)  [8 极点]", famos_bp, sosfilt(butter(4,[nw1,nw2],btype='bandpass',output='sos'), signal))
show("scipy butter(2, bandpass)  [4 极点]", famos_bp, sosfilt(butter(2,[nw1,nw2],btype='bandpass',output='sos'), signal))

# 手工：不预畸变，直接用模拟频率做 lp2bp，再双线性
def manual_bp(order, f1, f2, prewarp=True):
    z,p,k = buttap(order)
    if prewarp:
        ww1 = 4*np.tan(np.pi*(f1/(fs/2))/2); ww2 = 4*np.tan(np.pi*(f2/(fs/2))/2)
    else:
        ww1 = 2*np.pi*f1; ww2 = 2*np.pi*f2   # rad/s
    wo = np.sqrt(ww1*ww2); bw = ww2-ww1
    z,p,k = lp2bp_zpk(z,p,k,wo,bw)
    if prewarp:
        z,p,k = bilinear_zpk(z,p,k,2.0)
    else:
        z,p,k = bilinear_zpk(z,p,k,2*np.pi*fs)
    return zpk2sos(z,p,k)

show("手工 lp2bp(几何中心, 预畸变)", famos_bp, sosfilt(manual_bp(4,w1,w2,True), signal))
show("手工 lp2bp(几何中心, 不预畸变)", famos_bp, sosfilt(manual_bp(4,w1,w2,False), signal))
show("手工 lp2bp(几何中心, 不预畸变, 阶2)", famos_bp, sosfilt(manual_bp(2,w1,w2,False), signal))

print()
print("=== C. 巴特沃斯带阻 20/80Hz ===")
show("scipy butter(4, bandstop)  [8 极点]", famos_bs, sosfilt(butter(4,[nw1,nw2],btype='bandstop',output='sos'), signal))
show("scipy butter(2, bandstop)  [4 极点]", famos_bs, sosfilt(butter(2,[nw1,nw2],btype='bandstop',output='sos'), signal))

print()
print("=== D. 贝塞尔带通/带阻：哪种归一化 + 哪种阶数解释 ===")
fbp = read_bin(os.path.join(outdir,'famos_6.bin'))
fbs = read_bin(os.path.join(outdir,'famos_7.bin'))
for norm in ('phase','delay','mag'):
    show("bessel(4,bp,norm=%s) [8 极点]" % norm, fbp, sosfilt(bessel(4,[nw1,nw2],btype='bandpass',output='sos',norm=norm), signal))
for norm in ('phase','delay','mag'):
    show("bessel(2,bp,norm=%s) [4 极点]" % norm, fbp, sosfilt(bessel(2,[nw1,nw2],btype='bandpass',output='sos',norm=norm), signal))
for norm in ('phase','delay','mag'):
    show("bessel(4,bs,norm=%s) [8 极点]" % norm, fbs, sosfilt(bessel(4,[nw1,nw2],btype='bandstop',output='sos',norm=norm), signal))
