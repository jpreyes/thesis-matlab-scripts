# thesis-matlab-scripts
## APPENDIX D – MATLAB SCRIPTS DEVELOPED FOR THIS RESEARCH

The main MATLAB scripts developed for ground motion modeling, spectral acceleration analysis, and statistical evaluation of GMPEs are available in the following public repository:

🔗 **GitHub Repository:** [https://github.com/jpreyes/thesis-matlab-scripts](https://github.com/jpreyes/thesis-matlab-scripts)

### Included Scripts:

- `GMPE_MBR17.m` – Main script to compute spectral acceleration values (Sa) using the GMPE proposed by Montalva et al. (2017), with full vector support and log-log interpolation.
- `CH_MBR17_CoefModel.m` – Core implementation of the GMPE model, returning Sa and associated statistical components (σ, ϕ, τ) across periods.
- `CH_MBR17_PGA1000.m` – Auxiliary function to compute PGA at Vs30 = 1000 m/s, required for nonlinear site amplification when Vs30 < vlin.
- `GMPEFit.m` – Fits a smoothing spline to empirical data of ln(Sa) versus Rhyp, and generates residual plots for model validation.
- `generateSyntheticSeismicSignal1.m` – Script to generate a synthetic seismic acceleration time-history with Gaussian pulse envelope and white noise.
- `signal_analysis_roof_directions.m` – Reads and analyzes roof-level seismic responses in X, Y, and RZ directions using spectrograms, periodograms, and Yule-Walker estimators to assess modal energy evolution and damage effects.

All scripts are distributed under the **MIT License**.  
The repository also includes example datasets, coefficient tables, and plotting routines for visualization and model comparison.

├── GMPE/
│   ├── GMPE_MBR17.m
│   ├── CH_MBR17_CoefModel.m
│   └── ...
├── SignalAnalysis/
│   ├── signal_analysis_roof_directions.m
│   └── generateSyntheticSeismicSignal1.m
├── Models/
│   ├── Illapel_RC_Model_ETABS.s2k
│   ├── Illapel_RC_Model_Readme.md
│   └── ...
└── README.md
