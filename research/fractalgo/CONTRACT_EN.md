# FractalGo · computational contract

Bilingual software and documentation; manuscript and supplement distributed separately.

Command: `python run.py topology`.

Finite assertion: fixed-graph excursion connectivity is not determined by marginal Gaussian bands or permutation-invariant scalar summaries. Test the exact 80-node witness; a smooth conditional GP versus independent marginals on 51/101/201 nodes with 1024 draws; and conditional Cassini residual excursions on 1696 reconstructed epochs in 20 passes, with 512 predictive draws for each noise model.

GP graph: fixed ordered path, strictly positive cell prior masses, posterior-to-prior ratios, thresholds 0.5/1/2. Uniform and exp(-x) priors are explicit sensitivity scenarios. Same GP draws are restricted to nested grids. Node masses use trapezoidal weights; this differs from exact segment integration used for the earlier quantile analysis. No continuous-limit claim.

Cassini graph: chronological adjacency within pass; gaps exceeding twice that pass's median cadence break edges. Same graph for observed and replicated residuals. Signed thresholds 0/1/2 whitened-noise standard deviations. The inherited optical template, weights and nuisance basis remain conditioned upon. The four-kernel mixture and all new diagnostics are exploratory, not prospectively selected. Posterior predictive tail fractions are not uniformly calibrated p-values or force detections. No new station data, no physical Yukawa fit, no complete hierarchical SBC.

Retain every node ID, original row, timestamp, edge and input hash. No capture/deletion, collision displacement, invented dimensions or independent likelihood factors from derived summaries. The prototype does not estimate a fractal dimension or validate the entire FractalGo construct.
