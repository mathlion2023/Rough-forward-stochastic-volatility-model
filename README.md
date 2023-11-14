# Rough-forward-stochastic-volatility-model
Matlab codes for the paper: Implied Roughness in the Term Structure of Oil Markets Volatility
The data file is not shared.
Any difficulties in understanding the functionality of these codes, feel free to reach out to mesias@sun.ac.za

The rough forward stochastic volatility is expressed as
\[\begin{eqnarray}
dF(t,T, \mathbf{V_t})&=&F(t,T, \mathbf{V_t})\sum_{i=1}^n \sigma_i(t,T, \mathbf{V_t}) \left(\rho_{i}dW_i(t) + \sqrt{1-\rho^2_{i}}dW_i^V(t)\right)\label{Forwardheston} \\
 \mathbf{V_t}^i- \mathbf{V_0}^i&=&\frac{\kappa_i}{\Gamma(\alpha)}\int_0^t(t-s)^{\alpha-1}(\theta_i- \mathbf{V_s})dt +\frac{1}{\Gamma(\alpha)}\int_0^t(t-s)^{\alpha-1}\sigma_V^i\sqrt{ \mathbf{V_s}^i}dW^V_i,\label{RoughBM}
\end{eqnarray}
\]
where, $\alpha=H+\frac{1}{2}\in \left(\frac{1}{2},1\right]$ 
