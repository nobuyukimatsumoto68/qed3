import marimo

__generated_with = "0.23.6"
app = marimo.App()


@app.cell
def _():
    import os
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import scipy as sp
    import scipy.linalg
    import scipy.optimize
    import glob
    # # import lsqfit
    # from lib import Jackknife
    # from lib import simple_mean
    # from lib import format_print
    # plt.rcParams['figure.dpi'] = 300
    from natsort import natsorted

    return glob, math, natsorted, np, plt, sp


@app.cell
def _(math, np):
    class Jackknife:

        def __init__(self, len_data, binsize):
            self.binsize = binsize
            self.nbins = math.floor(len_data / self.binsize)
            self.N = self.binsize * self.nbins
            self.jack_avg = []
            self.est = 0
            self.var_est = 0

        def set(self, func, list_of_data):
            for i in range(self.nbins):
                self.jack_avg.append(func(i, self.binsize, list_of_data))

        def do_it(self):
            for i in range(0, self.nbins):
                self.est = self.est + self.jack_avg[i]
            self.est = self.est / self.nbins
            for i in range(0, self.nbins):
                self.var_est = self.var_est + (self.jack_avg[i] - self.est) ** 2
            self.var_est = self.var_est / self.nbins
            self.var_est = self.var_est * (self.nbins - 1)

        def mean(self):
            return self.est

        def var(self):
            return self.var_est

        def err(self):
            return np.sqrt(self.var_est)

    def simple_mean(i, binsize, np_data):
        resmpld = np.delete(np_data, np.s_[i * binsize:(i + 1) * binsize], axis=0)
        return np.mean(resmpld, axis=0)

    def jack_avg(i, binsize, np_data):
        return np_data[i]

    def take_abs(i, binsize, np_data):
        resmpld = np.delete(np_data, np.s_[i * binsize:(i + 1) * binsize], axis=0)
        mean = np.mean(resmpld, axis=0)
        return np.sqrt(mean.T[1] ** 2 + mean.T[2] ** 2)  # mean = mean[1:int(Nt/2)] + np.flip( mean[int(Nt/2):Nt-1], axis=0 )

    def take_real(i, binsize, np_data):
        resmpld = np.delete(np_data, np.s_[i * binsize:(i + 1) * binsize], axis=0)
        mean = np.mean(resmpld, axis=0)
        return mean.T[1]

    def take_imag(i, binsize, np_data):
        resmpld = np.delete(np_data, np.s_[i * binsize:(i + 1) * binsize], axis=0)  # mean = mean[1:int(Nt/2)] + np.flip( mean[int(Nt/2):Nt-1], axis=0 )
        mean = np.mean(resmpld, axis=0)
        return mean.T[2]

    def format_print(cen, err):
        for i in range(-50, 50):
            if 10 ** (-i + 1) >= err > 10 ** (-i):  # mean = mean[1:int(Nt/2)] + np.flip( mean[int(Nt/2):Nt-1], axis=0 )
                tmp = err * 10 ** (i + 1)
                return '{num:.{width}f}'.format(num=cen, width=i + 1) + '(' + str(round(tmp)) + ')'

    def format_print_w_exact(exact, cen, err):
        for i in range(-50, 50):
            if 10 ** (-i + 1) >= err > 10 ** (-i):
                tmp = err * 10 ** (i + 1)
                return '{num:.{width}f}'.format(num=cen, width=i + 1) + '(' + str(round(tmp)) + ')' + ': ' + '{num:.{width}f}'.format(num=(exact - cen) / err, width=i + 1) + ' sigma'

    return Jackknife, format_print_w_exact, jack_avg, simple_mean, take_real


@app.cell
def _():
    # pwd
    return


@app.cell
def _():
    # dir1='/projectnb/qfe/nmatsum/qed3/src/both_3d/'
    # dir1='/media/nobu/baracuda_14/qed3/qed3/src/both_3d/'
    dir1='/mnt/baracuda_14/qed3/qed3/src/both_3d/'
    return (dir1,)


@app.cell
def _():
    # results=[]
    return


@app.cell
def _():
    # Nf=2

    # # gsq="8.000000"
    # gsq="4.000000"
    # # gsq="1.000000"
    # # gsq="0.500000"
    # # gsq="0.100000"

    # # nu0="0.800000"
    # # nu0="2.000000"
    # nu0="2.000000"

    # nu1="3.000000" # nu0 # "1.000000"
    # # nu1="2.000000"
    # # nu1="1.000000"

    # at="0.200000"
    # Nt=96
    # L=1
    return


@app.cell
def _():
    # Nf=6

    # gsq="4.000000"
    gsq="8.000000"
    # gsq="2.800000"
    # gsq="0.500000"

    # gsq="1.000000"
    # gsq="0.500000"
    # gsq="0.100000"To

    # nu0="2.000000"
    # nu0="1.500000"
    nu0="1.000000"

    nu1="2.750000" # "1.000000"
    # nu1="2.000000"


    at="0.200000"
    Nt=128
    L=1

    at="0.200000"
    return L, Nt, at, gsq, nu0, nu1


@app.cell
def _(L, Nt, at, gsq, nu1):
    # label="data_at"+str(at)+"nu0"+str(nu0)+"nu1"+str(nu1)+"nt"+str(Nt)+"L"+str(L)
    # label="data_at"+str(at)+"nu0"+str(nu0)+"nu1"+str(nu1)+"nt"+str(Nt)+"L"+str(L)


    # label="data_Nf"+str(Nf)+"_gsq"+str(gsq)+"at"+str(at)+"nu0"+str(nu0)+"nt"+str(Nt)+"L"+str(L)

    label="data_gsq"+str(gsq)+"at"+str(at)+"nu1"+str(nu1)+"nt"+str(Nt)+"L"+str(L)
    return (label,)


@app.cell
def _(label):
    label
    return


@app.cell
def _(dir1, label):
    path = dir1+label
    path
    return (path,)


@app.cell
def _():
    binsize=1
    return (binsize,)


@app.cell
def _():
    ell=0
    em=0
    nhits = 1
    return ell, em, nhits


@app.cell
def _(
    Jackknife,
    Nt,
    binsize,
    ell,
    em,
    glob,
    natsorted,
    nhits,
    np,
    path,
    plt,
    take_real,
):
    ninit = 20
    nt0s = 1
    imax = 503
    corrs_aa = np.full((nhits, Nt, 3), 0.0)
    corrs_33 = np.full((nhits, Nt, 3), 0.0)
    for ab in range(1, 4):
        _files_ = sorted(glob.glob(path + '/meson_' + str(ell) + '_' + str(em) + '_a' + str(ab) + '_h*_corr_v2.*'), key=sorted)
        _files = np.array(_files_)
    corrs_aa = np.full((len(_files[ninit * nt0s:imax]), Nt, 3), 0.0)
    corrs_33 = np.full((len(_files[ninit * nt0s:imax]), Nt, 3), 0.0)
    for ab in range(1, 4):
        _files_ = sorted(glob.glob(path + '/meson_' + str(ell) + '_' + str(em) + '_a' + str(ab) + '_h*_corr_v2.*'), key=sorted)
        _files = np.array(_files_)
        corrs_ = []
        for _file in natsorted(_files)[ninit * nt0s:imax]:
            corr = np.loadtxt(_file)
            corrs_.append(corr)
        corrs = np.array(corrs_)
        corrs_aa = corrs_aa + corrs
        if ab == 3:
            corrs_33 = corrs_33 + corrs
    plt.plot(np.mean(corrs_aa, axis=1).T[1])
    plt.plot(np.mean(corrs_33, axis=1).T[1])
    jk_aa = Jackknife(corrs_aa.shape[0], binsize)
    jk_33 = Jackknife(corrs_33.shape[0], binsize)
    jk_aa.set(take_real, corrs_aa)
    jk_33.set(take_real, corrs_33)
    jk_aa.do_it()
    jk_33.do_it()
    return jk_33, jk_aa


@app.cell
def _(Nf, at, gsq, nu0, nu1):
    # desc="at"+str(at)+"_nu0"+str(nu0)+"_nu1"+str(nu1)
    desc="Nf"+str(Nf)+"_gsq"+str(gsq)+"_at"+str(at)+"_nu0"+str(nu0)+"_nu1"+str(nu1)
    return (desc,)


@app.cell
def _(Nt, at, np):
    xxx = float(at) * np.arange(Nt)
    return (xxx,)


@app.cell
def _(desc, jk_33, jk_aa, plt, xxx):
    plt.errorbar( xxx, jk_33.mean(), jk_33.err(), label="33" )
    plt.errorbar( xxx, -jk_aa.mean(), jk_aa.err(), label="m_aa" )
    plt.errorbar( xxx, jk_aa.mean(), jk_aa.err(), label="aa" )

    # plt.plot( xxx, 0.01*np.exp(-2.0*xxx) )
    # plt.plot( xxx, 0.01*np.exp(-3.0*xxx) )

    plt.yscale("log")
    plt.ylim(1.0e-13, 1.0e0)

    plt.legend()

    plt.savefig(desc+"JA_corr.pdf")
    return


@app.cell
def _(jk_33):
    jk_ell1em1_33 = jk_33
    return (jk_ell1em1_33,)


@app.cell
def _(jk_33):
    jk_ell1em0_33 = jk_33
    return (jk_ell1em0_33,)


@app.cell
def _(jk_aa):
    jk_ell0em0_aa = jk_aa
    return (jk_ell0em0_aa,)


@app.cell
def _(desc, jk_ell0em0_aa, jk_ell1em0_33, jk_ell1em1_33, plt, xxx):
    plt.errorbar( xxx, jk_ell1em1_33.mean(), jk_ell1em1_33.err(), label="11_33" )
    plt.errorbar( xxx, jk_ell1em0_33.mean(), jk_ell1em0_33.err(), label="10_33" )
    plt.errorbar( xxx, -jk_ell0em0_aa.mean()/3.0, jk_ell0em0_aa.err()/3.0, label="00_aa" )

    # plt.plot( xxx, 0.01*np.exp(-2.0*xxx) )
    # plt.plot( xxx, 0.01*np.exp(-3.0*xxx) )

    plt.yscale("log")
    plt.ylim(1.0e-13, 1.0e0)

    plt.legend()
    plt.savefig(desc+"JA_corr_factor_check.pdf")
    return


@app.cell
def _(jk_ell0em0_aa, jk_ell1em0_33, jk_ell1em1_33, plt, xxx):
    plt.errorbar( xxx, -jk_ell1em1_33.mean()/jk_ell0em0_aa.mean(), jk_ell1em0_33.err(), label="10_33" )
    return


@app.cell
def _(jk_ell0em0_aa, jk_ell1em0_33, plt, xxx):
    plt.errorbar( xxx, -jk_ell1em0_33.mean()/jk_ell0em0_aa.mean(), jk_ell1em0_33.err(), label="10_33" )
    return


@app.cell
def _(jks_):
    jks_[0][1][47]
    return


@app.cell
def _():
    # 1.3157697610054052e-09
    return


@app.cell
def _(jks_):
    -jks_[0][1][47]
    return


@app.cell
def _(Nt, at, desc, jks_, np, plt):
    # plt.errorbar(np.arange(1,Nt/2)*float(at), jks_[0][1], jks_[0][2], capsize=2.0, label='$\\Gamma=1$')

    plt.errorbar(np.arange(Nt)*float(at), np.abs(jks_[0][1])+jks_[0][1][47], jks_[0][2], capsize=2.0, label='$\\Gamma=1$')
    # plt.errorbar(np.arange(Nt)*float(at), jks_[1][1], jks_[1][2], capsize=2.0, label='$\\Gamma=\\sigma_1$')
    # plt.errorbar(np.arange(Nt)*float(at), jks_[2][1], jks_[2][2], capsize=2.0, label='$\\Gamma=\\sigma_2$')
    # plt.errorbar(np.arange(Nt)*float(at), jks_[3][1], jks_[3][2], capsize=2.0, label='$\\Gamma=\\sigma_3$')
    plt.yscale("log")
    plt.title(desc)

    # plt.plot( np.arange(Nt)*float(at), free.T[1], label='free', ls='dashed' )

    x = np.arange(1,Nt/2)*float(at)
    # plt.plot(x, 0.008*np.exp(-2.4*x), label="exp(-2.6x)")
    plt.plot(x, 0.02*np.exp(-2.0*x), label="exp(-2.0x)")

    # plt.ylim(1.0e-4, 1.0e-5)
    plt.legend()

    plt.savefig( desc+"_wall.pdf", bbox_inches="tight" )
    return


@app.cell
def _(Nt, at, np):
    # fitm = int(5.0 / 0.2)
    # fitM = int(8.0 / 0.2)

    fitm = int(1.2 / 0.2)
    fitM = int(5.0 / 0.2)


    def meff( cor ):
        avg = (cor[1:int(Nt/2)] + np.flip(cor[int(Nt/2)+1:]))/2
        tmp = np.arccosh(avg/cor[int(Nt/2)])
        return (-tmp[1:]+tmp[:-1])/float(at)
        # return -np.log(avg[1:]/avg[:-1])/float(at)

    xx = np.arange(1,int(Nt/2)-1)*float(at)

    # def fitter( x, m1, m2, A ):
    #     return np.full( x.shape[0], m1)

    def fitter( x, m1, m2, A ):
        return m1 + A * np.exp(-m2*x)

    # def meff( cor ):
    #     tmp = np.arccosh(l1_1/l1_1[int(Nt/2)])
    #     return (-tmp[1:]+tmp[:-1])/float(at)
    return fitM, fitm, fitter, meff, xx


@app.cell
def _(jk, meff):
    jk_M = []

    for l1_1 in jk.jack_avg:
        jk_M.append( meff(l1_1) )
    return (jk_M,)


@app.cell
def _(Jackknife, jack_avg, jk_M):
    jk_meff = Jackknife( len(jk_M), 1 )
    jk_meff.set( jack_avg, jk_M )
    jk_meff.do_it()
    return (jk_meff,)


@app.cell
def _(fitM, fitm, fitter, jk, jk_meff, meff, np, sp, xx):
    jk_mass = []
    for l1_1_1 in jk.jack_avg:
        yy = meff(l1_1_1)
        opt = sp.optimize.curve_fit(fitter, xx[fitm:fitM], yy[fitm:fitM], sigma=jk_meff.err()[fitm:fitM], p0=(np.sqrt(2.0), 2, 0.01), full_output=True, absolute_sigma=True)
        jk_mass.append(opt[0])
    return jk_mass, opt, yy


@app.cell
def _(Jackknife, binsize, jack_avg, jk_mass):
    jk_fit = Jackknife( len(jk_mass), binsize )
    jk_fit.set( jack_avg, jk_mass )
    jk_fit.do_it()
    return (jk_fit,)


@app.cell
def _(jk_fit):
    jk_fit.mean()
    return


@app.cell
def _(
    at,
    desc,
    ell,
    em,
    fitM,
    fitm,
    fitter,
    format_print_w_exact,
    jk_fit,
    jk_meff,
    plt,
    xx,
):
    plt.errorbar( xx, jk_meff.mean(), jk_meff.err(), label=str(fitm*float(at))+":"+str(fitM*float(at)) )

    plt.plot( xx, fitter(xx, jk_fit.mean()[0], jk_fit.mean()[1], jk_fit.mean()[2] ) )

    plt.ylim(0.0,5.0)

    plt.title(format_print_w_exact(2.0+ell, jk_fit.mean()[0], jk_fit.err()[0] ))
    plt.legend()

    plt.savefig( desc+"_fit"+str(ell)+str(em)+".pdf", bbox_inches="tight" )
    return


@app.cell
def _(opt):
    opt
    return


@app.cell
def _(nu0, results, yy):
    results.append( [nu0, yy] )
    return


@app.cell
def _(Nf, at, gsq):
    desc_1 = 'Nf' + str(Nf) + '_gsq' + str(gsq) + '_at' + str(at)
    return (desc_1,)


@app.cell
def _(desc_1, plt, results):
    for _result in results:
        plt.plot(_result[1], label='nu0=' + _result[0])
    plt.hlines(1.0, 0, 40)
    plt.legend()
    plt.ylim(0.5, 5.5)
    plt.savefig(desc_1 + '_meff_wall.pdf', bbox_inches='tight')
    return


@app.cell
def _(desc_1, plt, results):
    for _result in results:
        plt.plot(_result[1], label='nu0=' + _result[0])
    plt.hlines(1.0, 0, 40)
    plt.legend()
    plt.savefig(desc_1 + '_meff_wall.pdf', bbox_inches='tight')
    return


@app.cell
def _(results):
    results
    return


@app.cell
def _():
    # 1 0.138394
    # 2 0.073868
    # 4 2.994744726728921e-01
    return


@app.cell
def _(plt):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica",
        'font.size'  : 20
    })
    return


@app.cell
def _(np):
    def f(t, nmax, T):
        _res = 0.0
        for n in range(1, nmax):
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * t)
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * (T - t))
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * (t + T))
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * (2 * T - t))
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * (t + 2 * T))
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * (3 * T - t))
        return _res / (8.0 * np.pi)

    return (f,)


@app.cell
def _():
    Nt_1 = 96
    return


@app.cell
def _():
    gsq_1 = 0.5
    beta = 1.0 / gsq_1
    return (beta,)


@app.cell
def _():
    # gsq=2.00000
    # beta=1.0/gsq

    # key_1='gsq2.000000at0.500000nt24L1__corr'
    # key_2='gsq2.000000at0.500000nt24L2__corr'
    # # key_4='gsq0.500000at0.375000nt32L4_'

    # l1_1 = np.loadtxt("corr_"+key_1+".dat")
    # l1_2 = np.loadtxt("corr_"+key_2+".dat")
    # # l1_4 = np.loadtxt("corr_"+key_4+".dat")
    return


@app.cell
def _():
    # plt.yscale("log")

    # xx=np.arange(24)*0.5000
    # plt.errorbar( xx, l1_1.T[0], l1_1.T[1], capsize=5, fmt='x-', label="$L=1$" )
    # plt.errorbar( xx, l1_2.T[0], l1_2.T[1], capsize=5, fmt='v-', label="$L=2$" )
    # # plt.errorbar( xx, l1_4.T[0], l1_4.T[1], capsize=5, fmt='o-', label="$L=4$" )

    # xx=np.linspace(0, 12, 5000)

    # yy= (1.0/beta)*f(xx, 4000, 12) # *fixme*
    # plt.plot( xx, yy, ls='dashed', c='black', label='exact' )

    # yy= (1.0/beta)*f(xx, 2, 12) # *fixme*
    # plt.plot( xx, yy, ls='dashed', c='red', label='exact' )

    # plt.xlabel("$t$")
    # plt.title("$1/g^2 = 20, a_t=0.05, N_t=96$")

    # plt.xlim(0, 12)
    # plt.ylim(1.0e-5, 1.0e2)
    # plt.legend()

    # # plt.savefig("jtjt.pdf")
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.500000nt24L1__corr'
    _key_2 = 'gsq0.500000at0.500000nt24L2__corr'
    _key_4 = 'gsq0.500000at0.500000nt24L4__corr'
    l1_1_2 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2 = np.loadtxt('corr_' + _key_2 + '.dat')
    l1_4 = np.loadtxt('corr_' + _key_4 + '.dat')
    return l1_1_2, l1_2


@app.cell
def _(beta, f, l1_1_2, l1_2, np, plt):
    plt.yscale('log')
    xx_1 = np.arange(24) * 0.5
    plt.errorbar(xx_1, l1_1_2.T[0], l1_1_2.T[1], capsize=5, fmt='x-', label='$L=1$')
    plt.errorbar(xx_1, l1_2.T[0], l1_2.T[1], capsize=5, fmt='v-', label='$L=2$')
    xx_1 = np.linspace(0, 12, 5000)
    # plt.errorbar( xx, l1_4.T[0], l1_4.T[1], capsize=5, fmt='o-', label="$L=4$" )
    yy_1 = 1.0 / beta * f(xx_1, 4000, 12)
    plt.plot(xx_1, yy_1, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 2, a_t=0.5, N_t=24$')  # *fixme*
    plt.xlim(0, 12)
    plt.ylim(1e-05, 100.0)
    # yy= (1.0/beta)*f(xx, 2, 12) # *fixme*
    # plt.plot( xx, yy, ls='dashed', c='red', label='exact' )
    plt.legend()
    plt.savefig('jtjt_Nt24.pdf')
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.500000nt24L1__meff'
    _key_2 = 'gsq0.500000at0.500000nt24L2__meff'
    l1_1_3 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2_1 = np.loadtxt('corr_' + _key_2 + '.dat')
    return l1_1_3, l1_2_1


@app.cell
def _(l1_1_3, l1_2_1, np, plt):
    xx_2 = np.arange(24) * 0.5
    plt.errorbar(xx_2, l1_1_3.T[0], l1_1_3.T[1])
    plt.errorbar(xx_2, l1_2_1.T[0], l1_2_1.T[1])
    plt.hlines(np.sqrt(2.0), 0, 6, ls='dashed')
    plt.xlim(0, 6)
    plt.ylim(-2, 7)
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.375000nt32L1__corr'
    _key_2 = 'gsq0.500000at0.375000nt32L2__corr'
    l1_1_4 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2_2 = np.loadtxt('corr_' + _key_2 + '.dat')
    return l1_1_4, l1_2_2


@app.cell
def _(beta, f, l1_1_4, l1_2_2, np, plt):
    plt.yscale('log')
    xx_3 = np.arange(32) * 0.375
    plt.errorbar(xx_3, l1_1_4.T[0], l1_1_4.T[1], capsize=5, fmt='x-', label='$L=1$')
    plt.errorbar(xx_3, l1_2_2.T[0], l1_2_2.T[1], capsize=5, fmt='v-', label='$L=2$')
    xx_3 = np.linspace(0, 12, 5000)
    # plt.errorbar( xx, l1_4.T[0], l1_4.T[1], capsize=5, fmt='o-', label="$L=4$" )
    yy_2 = 1.0 / beta * f(xx_3, 4000, 12)
    plt.plot(xx_3, yy_2, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 2, a_t=0.375, N_t=32$')  # *fixme*
    plt.xlim(0, 12)
    plt.ylim(1e-05, 100.0)
    # yy= (1.0/beta)*f(xx, 2, 12) # *fixme*
    # plt.plot( xx, yy, ls='dashed', c='red', label='exact' )
    plt.legend()
    plt.savefig('jtjt_Nt32.pdf')
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.375000nt32L1__meff'
    _key_2 = 'gsq0.500000at0.375000nt32L2__meff'
    l1_1_5 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2_3 = np.loadtxt('corr_' + _key_2 + '.dat')
    return l1_1_5, l1_2_3


@app.cell
def _(l1_1_5, l1_2_3, np, plt):
    xx_4 = np.arange(32) * 0.375
    plt.errorbar(xx_4, l1_1_5.T[0], l1_1_5.T[1])
    plt.errorbar(xx_4, l1_2_3.T[0], l1_2_3.T[1])
    plt.hlines(np.sqrt(2.0), 0, 6, ls='dashed')
    plt.xlim(0, 6)
    plt.ylim(-2, 7)
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.250000nt48L1__corr'
    _key_2 = 'gsq0.500000at0.250000nt48L2__corr'
    _key_4 = 'gsq0.500000at0.250000nt48L4__corr'
    l1_1_6 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2_4 = np.loadtxt('corr_' + _key_2 + '.dat')
    l1_4_1 = np.loadtxt('corr_' + _key_4 + '.dat')
    return l1_1_6, l1_2_4, l1_4_1


@app.cell
def _():
    binsize_1 = 100
    return


@app.cell
def _(beta, f, l1_1_6, l1_2_4, l1_4_1, np, plt):
    plt.yscale('log')
    xx_5 = np.arange(48) * 0.25
    plt.errorbar(xx_5, l1_1_6.T[0], l1_1_6.T[1], capsize=5, fmt='x-', label='$L=1$')
    plt.errorbar(xx_5, l1_2_4.T[0], l1_2_4.T[1], capsize=5, fmt='v-', label='$L=2$')
    plt.errorbar(xx_5, l1_4_1.T[0], l1_4_1.T[1], capsize=5, fmt='o-', label='$L=4$')
    xx_5 = np.linspace(0, 12, 5000)
    yy_3 = 1.0 / beta * f(xx_5, 4000, 12)
    plt.plot(xx_5, yy_3, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 20, a_t=0.05, N_t=96$')  # *fixme*
    plt.xlim(0, 12)
    plt.ylim(1e-05, 100.0)
    # plt.savefig("jtjt.pdf")
    plt.legend()
    return


@app.cell
def _(l1_4_1):
    l1_4_1
    return


@app.cell
def _(l1_4_1):
    l1_4_1
    return


@app.cell
def _(beta, f, l1_1_6, l1_2_4, l1_4_1, np, plt):
    plt.yscale('log')
    xx_6 = np.arange(48) * 0.25
    plt.errorbar(xx_6, l1_1_6.T[0], l1_1_6.T[1], capsize=5, fmt='x-', label='$L=1$')
    plt.errorbar(xx_6, l1_2_4.T[0], l1_2_4.T[1], capsize=5, fmt='v-', label='$L=2$')
    plt.errorbar(xx_6, l1_4_1.T[0], l1_4_1.T[1], capsize=5, fmt='o-', label='$L=4$')
    xx_6 = np.linspace(0, 12, 5000)
    yy_4 = 1.0 / beta * f(xx_6, 4000, 12)
    plt.plot(xx_6, yy_4, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 20, a_t=0.05, N_t=96$')  # *fixme*
    plt.xlim(0, 12)
    plt.ylim(1e-05, 100.0)
    # plt.savefig("jtjt.pdf")
    plt.legend()
    return


@app.cell
def _():
    binsize_2 = 2000
    return


@app.cell
def _(beta, f, l1_1_6, l1_2_4, l1_4_1, np, plt):
    plt.yscale('log')
    xx_7 = np.arange(48) * 0.25
    plt.errorbar(xx_7, l1_1_6.T[0], l1_1_6.T[1], capsize=5, fmt='x-', label='$L=1$')
    plt.errorbar(xx_7, l1_2_4.T[0], l1_2_4.T[1], capsize=5, fmt='v-', label='$L=2$')
    plt.errorbar(xx_7, l1_4_1.T[0], l1_4_1.T[1], capsize=5, fmt='o-', label='$L=4$')
    xx_7 = np.linspace(0, 12, 5000)
    yy_5 = 1.0 / beta * f(xx_7, 4000, 12)
    plt.plot(xx_7, yy_5, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 20, a_t=0.05, N_t=96$')  # *fixme*
    plt.xlim(0, 12)
    plt.ylim(1e-05, 100.0)
    # plt.savefig("jtjt.pdf")
    plt.legend()
    return


@app.cell
def _(np):
    _key_1 = 'gsq0.500000at0.250000nt48L1__meff'
    _key_2 = 'gsq0.500000at0.250000nt48L2__meff'
    l1_1_7 = np.loadtxt('corr_' + _key_1 + '.dat')
    l1_2_5 = np.loadtxt('corr_' + _key_2 + '.dat')
    return l1_1_7, l1_2_5


@app.cell
def _(l1_1_7, l1_2_5, np, plt):
    xx_8 = np.arange(48) * 0.25
    plt.errorbar(xx_8, l1_1_7.T[0], l1_1_7.T[1])
    plt.errorbar(xx_8, l1_2_5.T[0], l1_2_5.T[1])
    plt.hlines(np.sqrt(2.0), 0, 6, ls='dashed')
    plt.xlim(0, 6)
    plt.ylim(-2, 7)
    return


@app.cell
def _():
    at_1 = 0.25
    return (at_1,)


@app.cell
def _(at_1, l1_1_7, np, plt):
    meff_1 = -np.log(l1_1_7.T[0][1:] / l1_1_7.T[0][:-1])
    plt.plot(meff_1[:24] / at_1)
    plt.hlines(np.sqrt(2.0), 0, 24, ls='dashed')
    return


@app.cell
def _(at_1, l1_2_5, np, plt):
    meff_2 = -np.log(l1_2_5.T[0][1:] / l1_2_5.T[0][:-1])
    plt.plot(meff_2[:24] / at_1)
    plt.hlines(np.sqrt(2.0), 0, 24, ls='dashed')
    return


@app.cell
def _():
    96
    return


@app.cell
def _():
    at_2 = 0.05
    beta_1 = 20
    key_2_1 = 'gsq0.050000at0.100000nt48L2_'
    key_2_2 = 'gsq0.050000at0.075000nt64L2_'
    key_2_3 = 'gsq0.050000at0.050000nt96L2_'
    return beta_1, key_2_1, key_2_2, key_2_3


@app.cell
def _(key_2_1, key_2_2, key_2_3, np):
    l1_1_8 = np.loadtxt('corr_' + key_2_1 + '.dat')
    l1_2_6 = np.loadtxt('corr_' + key_2_2 + '.dat')
    # l1_3 = np.loadtxt("corr_"+key_3+".dat")
    # l2 = np.loadtxt("corr_"+key+"nt96L2_short.dat")
    # l4 = np.loadtxt("corr_"+key+"nt96L4_short.dat")
    l1_4_2 = np.loadtxt('corr_' + key_2_3 + '.dat')
    return l1_1_8, l1_2_6, l1_4_2


@app.cell
def _(beta_1, f, l1_1_8, l1_2_6, l1_4_2, np, plt):
    plt.yscale('log')
    xx_9 = np.arange(48) * 0.1
    plt.errorbar(xx_9, l1_1_8.T[0], l1_1_8.T[1], capsize=5, fmt='x-', label='$L=1$')
    xx_9 = np.arange(64) * 0.075
    plt.errorbar(xx_9, l1_2_6.T[0], l1_2_6.T[1], capsize=5, fmt='v-', label='$L=2$')
    xx_9 = np.arange(96) * 0.05
    plt.errorbar(xx_9, l1_4_2.T[0], l1_4_2.T[1], capsize=5, fmt='o-', label='$L=4$')
    xx_9 = np.linspace(0, 4.8, 5000)
    yy_6 = 1.0 / beta_1 * f(xx_9, 4000, 4.8)
    plt.plot(xx_9, yy_6, ls='dashed', c='black', label='exact')
    plt.xlabel('$t$')
    plt.title('$1/g^2 = 20, a_t=0.05, N_t=96$')  # *fixme*
    plt.xlim(0, 5)
    plt.ylim(0.0001, 100.0)
    # plt.savefig("jtjt.pdf")
    plt.legend()
    return


@app.cell
def _(np):
    l1 = np.loadtxt("corr_beta5.000000at0.100000nt96L1.dat")
    l2 = np.loadtxt("corr_beta5.000000at0.100000nt96L2.dat")
    l4 = np.loadtxt("corr_beta5.000000at0.100000nt96L4.dat")
    return l1, l2, l4


@app.cell
def _(f, l1, l2, l4, np, plt):
    plt.yscale('log')
    xx_10 = np.arange(96) * 0.2
    plt.errorbar(xx_10, l1.T[0], l1.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_10, l2.T[0], l2.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_10, l4.T[0], l4.T[1], capsize=5, label='L=1')
    xx_10 = np.linspace(0, 3.0, 50)
    yy_7 = 1.0 / 5.0 * f(xx_10, 100)
    plt.plot(xx_10, yy_7)
    # yy= 0.02*np.exp( -np.sqrt(2.0) * xx * 2.4 )
    # yy= (1.0/300.0)*f(xx, 6)
    # plt.plot( xx, yy )
    # yy= (1.0/800.0)*f(xx, 11)
    # yy= (1.0/800.0)*f(xx, 20)
    plt.xlim(0, 4)
    return


@app.cell
def _(np):
    l1_3 = np.loadtxt('corr_beta20.000000at0.100000nt96L1.dat')
    l2_1 = np.loadtxt('corr_beta20.000000at0.100000nt96L2.dat')
    l4_1 = np.loadtxt('corr_beta20.000000at0.100000nt96L4.dat')
    return l1_3, l2_1, l4_1


@app.cell
def _(f, l1_3, l2_1, l4_1, np, plt):
    plt.yscale('log')
    xx_11 = np.arange(96) * 0.1
    plt.errorbar(xx_11, l1_3.T[0], l1_3.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_11, l2_1.T[0], l2_1.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_11, l4_1.T[0], l4_1.T[1], capsize=5, label='L=1')
    xx_11 = np.linspace(0, 3.0, 50)
    yy_8 = 1.0 / 20.0 * f(xx_11, 100)
    plt.plot(xx_11, yy_8)
    # yy= 0.02*np.exp( -np.sqrt(2.0) * xx * 2.4 )
    # yy= (1.0/20.0)*f(xx, 6)
    # plt.plot( xx, yy )
    # yy= (1.0/20.0)*f(xx, 11)
    # yy= (1.0/20.0)*f(xx, 20)
    plt.xlim(0, 4)
    return


@app.cell
def _(np):
    l1_5 = np.loadtxt('corr_beta50.000000at0.100000nt96L1.dat')
    l2_2 = np.loadtxt('corr_beta50.000000at0.100000nt96L2.dat')
    l4_2 = np.loadtxt('corr_beta50.000000at0.100000nt96L4.dat')
    return l1_5, l2_2, l4_2


@app.cell
def _(f, l1_5, l2_2, l4_2, np, plt):
    plt.yscale('log')
    xx_12 = np.arange(96) * 0.1
    plt.errorbar(xx_12, l1_5.T[0], l1_5.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_12, l2_2.T[0], l2_2.T[1], capsize=5, label='L=2')
    plt.errorbar(xx_12, l4_2.T[0], l4_2.T[1], capsize=5, label='L=4')
    xx_12 = np.linspace(0, 3.0, 50)
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('$G(t)$')
    # yy= 0.02*np.exp( -np.sqrt(2.0) * xx * 2.4 )
    yy_9 = 1.0 / 50.0 * f(xx_12, 100)
    # yy= (1.0/50.0)*f(xx, 6)
    # plt.plot( xx, yy )
    plt.plot(xx_12, yy_9)
    # yy= (1.0/50.0)*f(xx, 11)
    # yy= (1.0/50.0)*f(xx, 20)
    plt.xlim(0, 4)
    return


@app.cell
def _(np):
    l1_6 = np.loadtxt('corr_beta50.000000at0.100000nt96L1.dat')
    l2_3 = np.loadtxt('corr_beta50.000000at0.100000nt96L2.dat')
    l4_3 = np.loadtxt('corr_beta50.000000at0.100000nt96L4.dat')
    return l1_6, l2_3, l4_3


@app.cell
def _(f, l1_6, l2_3, l4_3, np, plt):
    plt.yscale('log')
    xx_13 = np.arange(96) * 0.1
    plt.errorbar(xx_13, l1_6.T[0], l1_6.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_13, l2_3.T[0], l2_3.T[1], capsize=5, label='L=1')
    plt.errorbar(xx_13, l4_3.T[0], l4_3.T[1], capsize=5, label='L=1')
    xx_13 = np.linspace(0, 3.0, 50)
    yy_10 = 1.0 / 50.0 * f(xx_13, 6)
    plt.plot(xx_13, yy_10)
    yy_10 = 1.0 / 50.0 * f(xx_13, 11)
    # yy= 0.02*np.exp( -np.sqrt(2.0) * xx * 2.4 )
    plt.plot(xx_13, yy_10)
    yy_10 = 1.0 / 50.0 * f(xx_13, 20)
    plt.plot(xx_13, yy_10)
    yy_10 = 1.0 / 50.0 * f(xx_13, 100)
    plt.plot(xx_13, yy_10)
    plt.xlim(0, 4)
    return


@app.cell
def _(f, np, plt, res_list):
    plt.yscale('log')
    for _res in res_list:
        Nt_2 = _res[1].shape[0]
        at_3 = float(_res[0].split('at')[1].split('nt')[0])
        xx_14 = np.arange(int(Nt_2 / 2)) * at_3
        yy_11 = _res[1][:int(Nt_2 / 2)]
        _dy = _res[2][:int(Nt_2 / 2)]
        plt.errorbar(xx_14, yy_11, _dy, capsize=5, label=_res[0])
    xx_14 = np.linspace(0, 3.0, 50)
    yy_11 = 0.0006 * f(xx_14, 10)
    plt.plot(xx_14, yy_11)
    plt.legend()
    plt.ylim(1e-05, 1.0)
    return


@app.cell
def _(l1_6):
    l1_6
    return


@app.cell
def _():
    res_list = []
    binsize_3 = 10
    return binsize_3, res_list


@app.cell
def _(Jackknife, binsize_3, glob, np, res_list, simple_mean):
    path_1 = './beta50.000000at0.138394nt32L1_v2/'
    _files_ = sorted(glob.glob(path_1 + 'plaq_ss_t_*.dat'), key=len)
    _files = np.array(_files_)
    _data_ = []
    for _file in _files:
        _data_.append(np.loadtxt(_file))
    _data = np.array(_data_)
    Nt_3 = _data.shape[1]
    jk = Jackknife(_data.shape[0], binsize_3)
    jk.set(simple_mean, _data)
    jk.do_it()
    res_list.append([path_1, jk.mean(), jk.err()])
    return (jk,)


@app.cell
def _(Jackknife, binsize_3, glob, np, res_list, simple_mean):
    path_2 = './beta50.000000at0.069197nt32L1_v2/'
    _files_ = sorted(glob.glob(path_2 + 'plaq_ss_t_*.dat'), key=len)
    _files = np.array(_files_)
    _data_ = []
    for _file in _files:
        _data_.append(np.loadtxt(_file))
    _data = np.array(_data_)
    Nt_4 = _data.shape[1]
    jk_1 = Jackknife(_data.shape[0], binsize_3)
    jk_1.set(simple_mean, _data)
    jk_1.do_it()
    res_list.append([path_2, jk_1.mean(), jk_1.err()])
    return


@app.cell
def _(Jackknife, binsize_3, glob, np, res_list, simple_mean):
    path_3 = './beta50.000000at0.036934nt32L2_v2/'
    _files_ = sorted(glob.glob(path_3 + 'plaq_ss_t_*.dat'), key=len)
    _files = np.array(_files_)
    _data_ = []
    for _file in _files:
        _data_.append(np.loadtxt(_file))
    _data = np.array(_data_)
    Nt_5 = _data.shape[1]
    jk_2 = Jackknife(_data.shape[0], binsize_3)
    jk_2.set(simple_mean, _data)
    jk_2.do_it()
    res_list.append([path_3, jk_2.mean(), jk_2.err()])
    return (Nt_5,)


@app.cell
def _(res_list):
    res_list
    return


@app.cell
def _(np):
    def f_1(t, nmax):
        _res = 0.0
        for n in range(1, nmax):
            _res = _res + np.sqrt(n * (n + 1.0)) * (2.0 * n + 1.0) * np.exp(-np.sqrt(n * (n + 1.0)) * t)
        return _res

    return (f_1,)


@app.cell
def _(Nt_5):
    Nt_5
    return


@app.cell
def _(f_1, np, plt, res_list):
    plt.yscale('log')
    for _res in res_list:
        Nt_6 = _res[1].shape[0]
        at_4 = float(_res[0].split('at')[1].split('nt')[0])
        xx_15 = np.arange(int(Nt_6 / 2)) * at_4
        yy_12 = _res[1][:int(Nt_6 / 2)]
        _dy = _res[2][:int(Nt_6 / 2)]
        plt.errorbar(xx_15, yy_12, _dy, capsize=5, label=_res[0])
    xx_15 = np.linspace(0, 3.0, 50)
    yy_12 = 0.0006 * f_1(xx_15, 10)
    plt.plot(xx_15, yy_12)
    plt.legend()
    plt.ylim(1e-05, 1.0)
    return


@app.cell
def _():
    50.0/8.0
    return


if __name__ == "__main__":
    app.run()
