# An apparent means of obtaining bathymetry from mass-conserved Manning's equation 



## Intro

Manning's McFli has an apparent contradiction. The system is such that unknowns vary solely in time or space, whereas observations vary along both time and space. Durand (2010) showed that the spatially varying unknown can be estimated  independently of the time-varying unknown, but GaMo and I have shown that $Q_t$ is unobtainable even if $A_0$ is known perfectly. What gives?! 



The following is an attempt to elucidate the differences in these approaches, and the conclusions reached by each one. 



### Commonalities across approaches

Both approaches have the following aspects in common.

- Apply Manning's equation at multiple times and locations
- Transform to linearize desired unknowns
- Apply conservation of mass
- Further manipulations to infer or show impossibility of inference



### Problem setup

Assume we have measurements of slope, width, and partial cross-section area for T times and N locations within a mass-conserved river segment. Then we can write


$$
\begin{aligned}
Q_{it} &= \frac{1}{n} W_{it}^{-2/3} A_{it}^{5/3}S_{it}^{1/2} \\
\implies Q^{3/5}_{it} & = A_{it}\Big(\frac{1}{n} W_{it}^{-2/3} S_{it}^{1/2} \Big)^{3/5}
\end{aligned}
$$
For simplicity let $x_{it} = \Big(\frac{1}{n} W_{it}^{-2/3} S_{it}^{1/2} \Big)^{3/5}$, so that $Q^{3/5}_{it} = A_{it}x_{it}$. 



### 2. Linearize

The linearization to isolate $Q$ is necessarily different from that to isolate $A_0$. This is at the heart of the difference in inferrability of these two quantities. 

From the above, we have
$$
Q^{3/5}_{it} = A_{it}x_{it} = (A_{0,i} + \delta A_{it}) x_{it}
$$


Note that this equation is linear in both $A_{0,i}$ and $Q_{it}$. 



### Mark's way of asserting mass conservation

If we take one index off of the Q term, then voila! we are mass-conserved.
$$
Q^{3/5}_{t} = A_{it}x_{it} = (A_{0,i} + \delta A_{it})x_{it}
$$
Alright, then
$$
\begin{aligned}
Q^{3/5}_{t} = (A_{0,i} + \delta A_{it})x_{it} \\
\implies Q^{3/5}_t / x_{it} - \delta A_{it} = A_{0,i}
\end{aligned}
$$


Now average over time. 
$$
Q^{3/5}_t / x_{it} - \delta A_{it} = A_{0,i} \\
\implies \bar{Q}(\frac{1 - \delta_{i \cdot}}{})
$$








Now, steady-state mass conservation requires that $Q_{it} = Q_{i't}$ for all $i, i'$ in $1, 2, ..., N$. This is equivalent to  
$$
\sum_{i = 1}^N\omega_i A_{it} x_{it} = 0
$$
for any set of $\omega_i$'s such that $\sum_{i = 1}^N \omega_i = 0$. In matrix notation, this is


$$
(\mathbf{A} \circ \mathbf{X}) \mathbf{\omega} = \mathbf{0}
$$
where $\circ$ denotes the elementwise product. We can now split this into the unobserved and observed parts of area:
$$
((\mathbf{A_0} + \mathbf{\delta A}) \circ \mathbf{X}) \mathbf{\omega} = (\mathbf{A_0} \circ \mathbf{X}) \omega + (\mathbf{\delta A} \circ \mathbf{X}) \mathbf{\omega} =  \mathbf{0} \\
\implies (\mathbf{A_0} \circ \mathbf{X}) \omega = (-\delta \mathbf{A} \circ \mathbf{X}) \omega
$$
A little math shows that the following representations of the LHS are equivalent:
$$
(\mathbf{A}_0 \circ \mathbf{X}) \omega = \mathbf{X} \mathbf{\Omega} \mathbf{a}_0 = \mathbf{X} (\mathbf{\omega} \circ \mathbf{a_0})
$$
where $\mathbf{\Omega}$ is the diagonal matrix whose diagonal entries are the elements of $\mathbf{\omega}$ and $\mathbf{a_0}$ is the the vector of $A_0$'s that is repeated to form the matrix $\mathbf{A_0}$, i.e. 


$$
\mathbf{A_0} = 

\begin{bmatrix} 

\mathbf{a_0} \\
\mathbf{a_0} \\
\vdots \\
\mathbf{a_0} \\

\end{bmatrix}

= 

\begin{bmatrix} 

A_{0,1} & A_{0,2} & ... & A_{0,N} \\
A_{0,1} & A_{0,2} & ... & A_{0,N} \\
\vdots \\
A_{0,1} & A_{0,2} & ... & A_{0,N} \\


\end{bmatrix}
$$


### Inference procedure

Given the above, the vector $\mathbf{a_0}$ can be inferred by solving the following linear model:


$$
\mathbf{X} \mathbf{\Omega} \mathbf{a}_0 = 1
$$


