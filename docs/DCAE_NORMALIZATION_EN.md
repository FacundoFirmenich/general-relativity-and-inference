# DCAE: normalization and thresholds

**Canonical v0.6 notation.** Here $R_j=1/\mathcal R_j$: $R_j$ is relative asymmetry and $\mathcal R_j$ is relative access. Therefore $I^\uparrow_j=[1-\mathcal R_j]_+$ and $\Delta_j=-\ln\mathcal R_j$. The earlier normalization is preserved exactly. See [v0.6 integration](DCAE_V06_EN.md).

**Thresholds imposed by normalization.** Defining $R=\mathcal A(s)/\mathcal A(s_0)>0$, the logarithm cancels in $I=\max\{0,1-R^{-1}\}$. For $0<c<1$, $I=c$ is equivalent to $R=(1-c)^{-1}$ and $\Delta=-\ln(1-c)$. In particular, $I=0.90$ corresponds exactly to $R=10$ and $\Delta=\ln10$, whereas $I=0.95$ corresponds to $R=20$. These values follow from the definition; they are not independent evidence of an information barrier or physical transition. Both $R$ and $\Delta$ must also be reported with their baseline, uncertainties and covariances; proximity of $I$ to one does not measure confidence or statistical significance.

## Normalization, base conversion and constructed thresholds

For $R=\mathcal A(s)/\mathcal A(s_0)>0$, the identities

$$\Delta=\ln R,\qquad I=[1-R^{-1}]_+,\qquad I=c\Longleftrightarrow R=(1-c)^{-1},\quad 0<c<1$$

follow by substitution, without an additional dynamical hypothesis. The values $\ln10$ at 90% and $\ln20$ at 95% are coordinates of chosen thresholds. They are not empirical constants discovered by the index. The axis $s$ need not be time, and this normalization does not establish an exponential relaxation law.

For a base $b>1$, $\Delta_b=\log_bR$ and $\Delta=(\ln b)\Delta_b$. The consistent transformation is $I=1-b^{-\max(0,\Delta_b)}$. Replacing $\ln$ by $\log_b$ while retaining a base-$e$ exponential without converting its exponent changes the index. Obtaining 2.302585 in natural coordinates and 1 in decimal coordinates for $R=10$ does not identify a different physical mechanism.

On the branch $R>1$, $dI/dR=R^{-2}$ and $dI/d\Delta=1-I$. The display scale flattens near one: $R=10,20,100$ yields $I=0.90,0.95,0.99$. This transformation is nevertheless bijective on that branch; in exact arithmetic, with the probability law correctly transformed, flattening alone does not remove information. Loss occurs when all $R\leq1$ are clipped to zero, through rounding or numerical saturation, or when only a threshold is retained. These operations must be distinguished from decreasing sensitivity of a physical instrument. Retain the signed gap and source scales; in extreme regimes report $\ln R$ to avoid overflow, without inventing a finite $R$.

Uncertainty propagates from the joint law of the four scales or their logarithms, including shared-baseline covariances. $I=0.90$ means neither a 90% probability nor a 90% loss of Fisher information. If $P(R>10\mid Y)$ is computed, it is a model- and data-dependent posterior probability, distinct from the descriptor value. Since $I$ is a deterministic function of $R$, treating them as independent observations would duplicate evidence. A physical claim needs an observable, null model and validation independent of the normalization identity.

The requirement applies symmetrically to DCAE and the mantissa and Cassini audit: proximity to $\ln10$ neither proves an artifact nor validates a physical barrier when it appears in our index. The tool `research/dcae/normalization_audit.py` checks the logarithm-free form, thresholds, base conversion, positive-branch inverse, clipping and sensitivity. Its outputs are numerical mathematical controls, not cosmological measurements.

