# One-sided moment-prior calculations: manuscript comparison

This historical check used package revision `dfc8343` and compared it with Table 1 (p. 5), Sections 2–3.2 and Appendix B of Pawel and Held, [Bayes Factor Group Sequential Designs, arXiv:2601.02851v2](https://arxiv.org/pdf/2601.02851v2) (22 June 2026).

I agree with the relevant mathematics under the stated normal/canonical sampling model. The added fixed-sample and sequential calculations implement that model; I found no statistical discrepancy. This is a check of the added moment-prior calculations, not a replication of every application or the paper's t-test approximation.

For `alternative = "greater"`, set

$$
z=(\widehat\theta-\theta_0)/\sigma,\quad r=\tau/\sigma,
\quad q=z/\sqrt{1+1/r^2},\quad A=(1+r^2)^{3/2}.
$$

Here `se` is $\sigma$, `psd` is the normal-kernel scale $\tau$, and `null` is $\theta_0$. Multiplying the Table 1 numerator and denominator by $e^{-q^2/2}$ gives

$$
\mathrm{BF}_{01}
=\frac{A}{2e^{q^2/2}(1+q^2)\Phi(q)+q\sqrt{2/\pi}}
=\frac{Ae^{-q^2/2}}{2\{(1+q^2)\Phi(q)+q\phi(q)\}},
$$

which is the formula in `nmbf01`. It corresponds to the normalized prior $2\delta^2 N(\delta\mid0,\tau^2)/\tau^2$ on $\delta=\theta-\theta_0>0$; the null remains a point. The `less` alternative reflects $z$.

The predictive distribution agrees with Eq. (3.1): mean $(\mathrm{dpm}-\theta_0)/\sigma_i$, covariance $\sqrt{I_{\min}/I_{\max}}+\mathrm{dpsd}^2/(\sigma_i\sigma_j)$, where $I_i=1/\sigma_i^2$. The design prior remains untruncated. Stopping probabilities require continuation at every earlier look; inconclusive trials still consume the final planned sample size when calculating expected sample size and variance.

The differences are computational and API extensions:

- The printed expression is unsuitable for direct floating-point evaluation in extreme tails. The implementation uses log arithmetic, bounded scale ratios, a small-$q$ expansion and a continued fraction for the negative tail. These evaluate the same integral; they do not change the prior or evidence threshold.
- Numerical root-finding is retained, with finite-range bracketing. Nonzero null values and the `less` direction extend the paper's displayed positive, zero-null formula by translation/reflection.
- Multivariate-normal integration uses a finite deterministic grid. It approximates the same stopping probabilities; accuracy depends on integration settings. Normal/canonical design assumptions remain essential.

Independent verification ([R script](scripts/verify_one_sided_moment_manuscript.R), run against the installed branch package):

| Check | Result |
| --- | --- |
| Table 1 and independent marginal integral versus `nmbf01` | 240 grid rows (60 unique r,z pairs for the printed expression); maximum absolute log error 2.9e-13 for the printed expression and 2.6e-13 for the implementation |
| Fixed-sample power versus independent boundaries/normal tails | 96 fixed-power and 48 one-look probabilities; maximum absolute error 1.8e-12 |
| Two-look stopping probabilities and sample-size moments versus conditional-normal integration | 8 two-look designs, ngrid=50000: probability error <=1.6e-6; expected-size error <=5.3e-11; variance error <=2.1e-10 |

At `dfc8343`, 1,076 default assertions passed and `R CMD check --no-manual --ignore-vignettes` returned Status OK. That run omitted extended tests and vignette and manual regeneration.

The same check recorded a search limitation. For `k=.1`, target power `.016`, `usd=1`, `psd=.3`, `dpm=dpsd=.03`, and `alternative="greater"`, `nnmbf01` returned **405**, although **18** already met the target. The power curve is nonmonotone, so the numerical root search can return a later crossing. The call `nbf01seq(..., looks=1, search="exhaustive", nrange=c(2,1000), type="moment")` returned 18.

This is a search-policy limitation, not evidence against Table 1 or Eq. (3.1).

To run the check, install this branch's `package` directory, then run `Rscript simulations/scripts/verify_one_sided_moment_manuscript.R path/to/results`. The script writes detailed CSV/RDS results and reports its loaded package path. The PDF checked has SHA256 `ddc13db23f215fa7cb21efc988e78f7277584e2163fa7ad50c0993f3666a9880`.
