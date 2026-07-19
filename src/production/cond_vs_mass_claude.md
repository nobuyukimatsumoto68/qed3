# Condensate vs mass (redo massive, Nf2, e/o+spin dilution)

Convention (jj_ylm_condensate notebook): norm = Nt*4pi = 128*4pi;
sigma_PS = 2*Re(etadag_xi)/norm + 2  (contact-sub +2);
sigma_FS = Re(etadag_xi - xidag_1mDdag_eta)/norm + 2 - m  (contact-sub +2-m).
m = physical mass (R=1), SAME convention across L (measure weighting only enters the local
operator m_L = m*A_y/abar_s; the ensemble label m is identical across L).
Config-jackknife; NOTE thin statistics (3-6 configs/point, stride 20).

| L | gsq | m | ncfg | sigma_PS | -sigma_FS |
|---|-----|---|------|----------|-----------|
| 1 | 1.5 | 0.10 | 6 | 0.0736(0.0001) | 0.0089(0.0000) |
| 1 | 1.5 | 0.20 | 6 | 0.1434(0.0001) | 0.0243(0.0000) |
| 1 | 1.5 | 0.30 | 6 | 0.2096(0.0004) | 0.0459(0.0001) |
| 1 | 1.5 | 0.40 | 6 | 0.2730(0.0006) | 0.0733(0.0001) |
| 2 | 3.0 | 0.10 | 4 | 0.0397(0.0002) | 0.0502(0.0000) |
| 2 | 3.0 | 0.20 | 4 | 0.0790(0.0003) | 0.1024(0.0000) |
| 2 | 3.0 | 0.30 | 4 | 0.1163(0.0004) | 0.1565(0.0001) |
| 2 | 3.0 | 0.40 | 4 | 0.1534(0.0006) | 0.2126(0.0001) |
| 4 | 6.0 | 0.10 | 3 | 0.0197(0.0001) | 0.0741(0.0000) |
| 4 | 6.0 | 0.50 | 3 | 0.0967(0.0004) | 0.3756(0.0001) |
| 4 | 6.0 | 1.00 | 3 | 0.1906(0.0008) | 0.7637(0.0002) |
| 4 | 6.0 | 1.50 | 3 | 0.2791(0.0010) | 1.1631(0.0002) |
