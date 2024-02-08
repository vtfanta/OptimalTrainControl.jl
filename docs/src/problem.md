# Problem Statement

The optimal train control problem which is solved by this package is mathematically formulated as
```math
\mathrm{min}\ J=\mathrm{min} \int\limits_0^X \left( \frac{u+|u|}{2} + \rho\frac{u-|u|}{2} \right)\mathrm{d}x
```
subject to the constraints
```math
\begin{gather*}
t' = \frac1v,\quad v' = \frac{u-r(v)+g(x)}{v},\\
t(0) = 0,\quad t(X) = T,\\
v(0) = v_i,\quad v(X) = v_f,\\
\overline{U}(v) \leq u \leq \underbar{U}(v).
\end{gather*}
```
In other words, calculate control signal $u$ subject to the dynamical system (the constraints on the first line above)
and the boundary constraints given by initial speed $v_i$, final speed $v_f$ and total journey time $T$ such that the
criterion $J$ is minimised. The control signal $u$ is constrained itself by the train traction characteristics
$\overline{U}(v) \leq u \leq \underbar{U}(v)$. 

The resistance function $r(v)$ which combines the influences of mechanic and
aerodynamic resistances is modeled as $r(v) = A + Bv + Cv^2$ with $A,B,C\geq0$. The speed parameters
$v_i$ and $v_f$ have to be strictly positive (it is highly recommended to use values $v\geq1$) since
the state equations contain singularity at $v=0$. The control limits have to be monotone ($\overline{U}$ non-increasing,
$\underbar{U}$ non-decreasing) with
```math
\lim_{v\to\infty}\underbar{U}(v)=\lim_{v\to\infty}\overline{U}(v) = 0.
```

It can be shown that the optimal control is a switching control between a small number of modes: Maximum
acceleration, holding constant speed, coasting (no braking and no traction) and maximum braking.

The output of the implemented solution are the trajectories of $t$, $v$ and $u$ as well as the location of switching points between the control modes.