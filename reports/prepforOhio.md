Questions for the working group following Pepsi Challenge 2.0 release:



- What is overall expected error statitsics for discharge estimation?
- Which cases perform better or worse?
- In which cases are statistical assumptions most closely followed?
  - $A_0$ log-normality
  - Error autocorrelation in time and space



Recommendations for new cases:

- Wide range of flow conditions, in order to assure that sample statistics closely approximate global (in time) parameters.
  - sample mean q is approximately global mean q
  - Maximize entropy of data



What is needed:

- Quick way to scope data interactively
- Ditto for method validation
- Assessment of errors in water-balance estimates





Bayes in info-theory form:


$$
\begin{aligned}
H(\theta | \mathbf{y}) &= H(\theta, \mathbf{y}) - H(\mathbf{y}) \\
&= H(\theta) - I(\theta; \mathbf{y}) \\
&= H(\mathbf{y} | \theta) + H(\theta)  - H(\mathbf{y})
\end{aligned}
$$
Now thinking about conditional entropy for various QOI's:

- $\bar{q}$
- $\dot{q}_t$
- $\tilde{n}$



As an example, take the "optimistic" case, and assume multivariate normal.


$$
\begin{aligned}
H(\bar{q} | \mathbf{y}) &= H(\bar{q}, \mathbf{y}) - H(\mathbf{y}) \\
&= H(\mathbf{y} | \theta) + H(\theta)  - H(\mathbf{y})
\end{aligned}     
$$


