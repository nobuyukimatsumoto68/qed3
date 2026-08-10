# effmass_full_claude.py
# Effective-mass analysis of the redo massless fermion correlators, l=0..3, config-jackknife.
#  - VECTOR tp (s3): PHYSICAL full = -C_conn + (disc_DCsub - plateau), ipynb subtraction scheme
#    (Eq.12 PhysRevD.110.054501 DC-subtract + residual plateau over dt in [16,Nt/2]); conn-only overlaid.
#    Disc is VECTOR-only (h0/disc/ylm/...); matched conn/disc configs.
#  - AXIAL tp (s3): conn-only (no axial disc).   - SCALAR PS/FS: conn-only.
# cosh effmass vs dt (timeslice separation); jackknife error band.  Stable window dt~[8,16] shaded.
# One png per (channel-group, L, gsq); Nf=2,4,6 overlaid.  kmin=20.  L1/L2 (L4 disc too sparse).
import glob, math, re
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

Nt=128; at=0.2; KMIN=20; inv4pi=1.0/(4.0*math.pi)
GS={1:[0.5,1.0,1.5],2:[1.0,2.0,3.0],3:[1.5,3.0,4.5],4:[2.0,4.0,6.0]}; HB={1:"1.000000",2:"1.000000",3:"0.400000-1.000000",4:"0.400000-1.000000"}
NFS=[2,4,6]; nfcol={2:"tab:red",4:"tab:blue",6:"tab:green"}
LSET=[0,1,2,3]
PWIN=np.arange(16,Nt//2+1)          # disc plateau window (dt)
WIN=(8,16)                          # stable effmass window (shaded)
dtp=np.arange(1,Nt//2)

def esn(nf,L,g): return "data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s_vmRe0.000000vmIm0.000000/"%(nf,g,L,HB[L])
def kof(fn): return int(fn.rsplit('/',1)[-1].split('.')[1])
def conn_files(nf,L,g):
    fs=sorted(glob.glob(esn(nf,L,g)+"corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if kof(f)>=KMIN]
def disc_files(nf,L,g):
    fs=sorted(glob.glob(esn(nf,L,g)+"corr_ylm_disc_tb2/corr.*.h0.h5"))
    return [f for f in fs if kof(f)>=KMIN]
def matched(nf,L,g):
    cfm={kof(f):f for f in conn_files(nf,L,g)}; dfm={kof(f):f for f in disc_files(nf,L,g)}
    ks=sorted(set(cfm)&set(dfm)); return [cfm[k] for k in ks],[dfm[k] for k in ks]
def load(fs,key):
    o=[]
    for fn in fs:
        with h5py.File(fn,'r') as f: o.append(f[key+'/real'][()]+1j*f[key+'/imag'][()])
    return np.array(o)

def gl_conn(fs,grp,l):   # m-averaged conn (Vpp), per config
    return sum(load(fs,"h0/%s/s3/l%d/m%d/Vpp"%(grp,l,m)) for m in range(-l,l+1))/(2*l+1)
def gl_scalar(fs,l,fs_ch):
    vpp=sum(load(fs,"h0/scalar/l%d/m%d/Vpp"%(l,m)) for m in range(-l,l+1))/(2*l+1)
    if fs_ch=="PS": return 2.0*vpp.real
    vfs=sum(load(fs,"h0/scalar_fs/l%d/m%d/Vmm"%(l,m)) for m in range(-l,l+1))/(2*l+1)
    return (vpp+vfs).real

def two_point_pp(Jc):    # C(dt)=(1/Nt) sum_t0 J(t0) J((t0+dt)%Nt), per config; Jc=(H,Nt)
    G=np.zeros(Jc.shape,dtype=complex)
    for dt in range(Nt): G[:,dt]=np.mean(Jc*np.roll(Jc,-dt,axis=1),axis=1)
    return G
def disc_samples(dfs,l):
    Jm=[load(dfs,"h0/disc/ylm/s3/l%d/m%d/J"%(l,m)) for m in range(-l,l+1)]
    return inv4pi*sum(two_point_pp(J) for J in Jm)/(2*l+1)

def vec_full(nf,L,g,l):   # returns (conn_cfg, full_cfg) real (H,Nt) or None
    cfm,dfm=matched(nf,L,g)
    if len(cfm)<4: return None
    conn=-gl_conn(cfm,'ylm',l).real                      # -C_conn
    samp=disc_samples(dfm,l).real                        # C_disc
    sub=samp-samp.mean(axis=1,keepdims=True)             # Eq.12 DC-subtracted
    plat=sub[:,PWIN].mean(axis=1)                        # per-config residual plateau
    full=conn+(sub-plat[:,None])                         # -C_conn + (disc - plateau) -> decays
    return conn, full

def meff_acosh(C):
    m=np.full(Nt,np.nan)
    for t in range(1,Nt-1):
        d=2.0*C[t]
        if d==0: continue
        r=(C[t-1]+C[t+1])/d
        if r>1.0: m[t]=np.arccosh(r)
    return m
def effmass_jk(samp):   # samp (H,Nt) real correlator -> mean,err of cosh effmass (1/a_t)
    H=samp.shape[0]; jk=(samp.sum(0)-samp)/(H-1)
    me=np.array([meff_acosh(jk[i]) for i in range(H)])/at
    mean=np.nanmean(me,0); err=np.sqrt(np.maximum((H-1)*np.nanmean((me-mean)**2,0),0.0))
    return mean,err

# ---------- VECTOR full vs conn ----------
for L in [1,2,3,4]:
    for g in GS[L]:
        fig,axs=plt.subplots(2,2,figsize=(13,9))
        drew=False
        for l,ax in zip(LSET,axs.ravel()):
            for nf in NFS:
                r=vec_full(nf,L,g,l)
                if r is None: continue
                conn,full=r
                cm,ce=effmass_jk(conn.astype(float)); fm,fe=effmass_jk(full.astype(float))
                ax.errorbar(dtp,cm[dtp],yerr=ce[dtp],fmt='o',ms=2.5,mfc='none',capsize=1.5,color=nfcol[nf],alpha=0.5,label="Nf%d conn"%nf)
                ax.errorbar(dtp,fm[dtp],yerr=fe[dtp],fmt='o',ms=2.5,capsize=1.5,color=nfcol[nf],label="Nf%d full"%nf)
                drew=True
            ax.axvspan(WIN[0],WIN[1],color="gray",alpha=0.1)
            ax.set_title("vector tp l=%d"%l); ax.set_xlim(0,30); ax.set_ylim(0,6)
            ax.set_xlabel("dt"); ax.set_ylabel(r"$m_\mathrm{eff}$ ($1/a_t$)"); ax.grid(alpha=0.3); ax.legend(fontsize=6,ncol=2)
        fig.suptitle("VECTOR tp: full (conn+disc, DC+plateau sub) vs conn -- L%d g%.1f (window dt[8,16])"%(L,g))
        fig.tight_layout()
        if drew: fig.savefig("figs/effmass_vecfull_L%d_g%.1f_claude.png"%(L,g),dpi=140)
        plt.close(fig)
        print("vec L%d g%.1f done"%(L,g))

# ---------- AXIAL (conn) and SCALAR (conn), l=0..3 ----------
for grp,name,func in [("axial","axial tp (s3) conn",None),("scalar","scalar PS/FS conn",None)]:
    for L in [1,2,3,4]:
        for g in GS[L]:
            fig,axs=plt.subplots(2,2,figsize=(13,9)); drew=False
            for l,ax in zip(LSET,axs.ravel()):
                for nf in NFS:
                    fs=conn_files(nf,L,g)
                    if len(fs)<4: continue
                    if grp=="axial":
                        C=2.0*gl_conn(fs,'ylm_axial',l).real
                        m,e=effmass_jk(C.astype(float))
                        ax.errorbar(dtp,m[dtp],yerr=e[dtp],fmt='o',ms=2.5,capsize=1.5,color=nfcol[nf],label="Nf%d"%nf); drew=True
                    else:
                        for ch,mk in [("PS",'o'),("FS",'^')]:
                            C=gl_scalar(fs,l,ch)
                            m,e=effmass_jk(C.astype(float))
                            ax.errorbar(dtp,m[dtp],yerr=e[dtp],fmt=mk,ms=2.5,capsize=1.5,color=nfcol[nf],mfc='none' if ch=="FS" else None,label="Nf%d %s"%(nf,ch)); drew=True
                ax.axvspan(WIN[0],WIN[1],color="gray",alpha=0.1)
                ax.set_title("%s l=%d"%(name,l)); ax.set_xlim(0,30); ax.set_ylim(0,6)
                ax.set_xlabel("dt"); ax.set_ylabel(r"$m_\mathrm{eff}$ ($1/a_t$)"); ax.grid(alpha=0.3); ax.legend(fontsize=6,ncol=2)
            fig.suptitle("%s -- L%d g%.1f"%(name,L,g)); fig.tight_layout()
            tag="axial" if grp=="axial" else "scalar"
            if drew: fig.savefig("figs/effmass_%s_L%d_g%.1f_claude.png"%(tag,L,g),dpi=140)
            plt.close(fig)
            print("%s L%d g%.1f done"%(tag,L,g))
