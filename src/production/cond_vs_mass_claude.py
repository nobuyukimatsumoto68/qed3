# cond_vs_mass_claude.py
# Condensate vs physical mass, Nf2, all L at the SAME masses m=0.1,0.2,0.3,0.4 (R=1).
# norm=Nt*4pi; sigma_PS=2Re(etadag_xi)/norm+2 (contact +2);
# sigma_FS=Re(etadag_xi - xidag_1mDdag_eta)/norm+2-m (contact +2-m). Config-jackknife.
# L1: gsq1.5 hb1.0; L2: gsq3.0 hb1.0; L4: gsq6.0 (resc-shifted hb per mass). kmin=20.
import glob, math, re
import numpy as np, h5py, pickle
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Nt=128; norm=Nt*4.0*math.pi; KMIN=20
MASSES=[0.1,0.2,0.3,0.4]
# L -> (gsq, hb-per-mass dict).  L1/L2 uniform hb 1.000000; L4 resc-shifted.
L4HB={0.1:"0.425902-1.025902",0.2:"0.451804-1.051804",0.3:"0.477706-1.077706",0.4:"0.503608-1.103608"}
gsqL={1:1.5,2:3.0,4:6.0}
lcol={1:"tab:red",2:"tab:blue",4:"tab:green"}; lmk={1:"o",2:"s",4:"^"}

def ens_dir(L,mre):
    g=gsqL[L]; ms="%.6f"%mre
    hb=L4HB[round(mre,1)] if L==4 else "1.000000"
    return "data_Nf2_gsq%.6fat0.200000nu01.000000mRe%smIm0.000000nt128L%d_hb%s_vmRe%svmIm0.000000/corr_condensate_eo_nhits1/"%(g,ms,L,hb,ms)
def kof(fn): return int(fn.rsplit('/',1)[-1].split('.')[1])
def find(L,mre):
    return [f for f in sorted(glob.glob(ens_dir(L,mre)+"corr.*.h0.h5")) if kof(f)>=KMIN]
def jk(a):
    n=len(a)
    if n<2: return np.nan,np.nan
    m=a.mean(); sub=(n*m-a)/(n-1); return m, math.sqrt((n-1)/n*((sub-m)**2).sum())

rows=[]
for L in [1,2,4]:
    for mre in MASSES:
        fs=find(L,mre)
        if len(fs)<2:
            rows.append((L,mre,len(fs),np.nan,np.nan,np.nan,np.nan)); continue
        xi=[];xid=[]
        for fn in fs:
            with h5py.File(fn,'r') as f:
                xi.append(f['h0/condensate/etadag_xi/real'][0]); xid.append(f['h0/condensate/xidag_1mDdag_eta/real'][0])
        xi=np.array(xi); xid=np.array(xid)
        ps=2*xi/norm+2.0; fsv=(xi-xid)/norm+2.0-mre
        m1,e1=jk(ps); m2,e2=jk(fsv)
        rows.append((L,mre,len(fs),m1,e1,m2,e2))
pickle.dump(rows,open("/tmp/cond_rows.pkl","wb"))

fig,axs=plt.subplots(1,2,figsize=(12,5))
for L in [1,2,4]:
    d=[r for r in rows if r[0]==L and not np.isnan(r[3])]
    if not d: continue
    ms=[r[1] for r in d]
    axs[0].errorbar(ms,[r[3] for r in d],[r[4] for r in d],marker=lmk[L],color=lcol[L],capsize=3,label="L%d g%.1f"%(L,gsqL[L]))
    axs[1].errorbar(ms,[-r[5] for r in d],[r[6] for r in d],marker=lmk[L],color=lcol[L],capsize=3,label="L%d g%.1f"%(L,gsqL[L]))
axs[0].set_ylabel(r"$\sigma_{PS}$"); axs[1].set_ylabel(r"$-\sigma_{FS}$")
for ax in axs:
    ax.set_xlabel(r"$m$ (physical, R=1)"); ax.grid(alpha=0.3); ax.legend(fontsize=9); ax.axhline(0,color="gray",lw=0.5)
axs[0].set_title("Nf2 PS condensate (matched m)"); axs[1].set_title("Nf2 FS condensate (matched m)")
fig.tight_layout(); fig.savefig("cond_vs_mass_claude.png",dpi=150)
print("wrote cond_vs_mass_claude.png")
for L,mre,n,m1,e1,m2,e2 in rows:
    s = "%.4f(%.4f) %.4f(%.4f)"%(m1,e1,-m2,e2) if not np.isnan(m1) else "-- (ncfg=%d)"%n
    print("L%d m=%.1f ncfg=%d : %s"%(L,mre,n,s))
