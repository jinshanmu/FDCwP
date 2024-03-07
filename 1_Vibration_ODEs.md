# 1 Vibration ODEs

## 1.1 Finite Difference Discretization

### 1.1.1 A Basic Model for Vibrations

Consider the vibrating mechanical system
$$u''+\omega^2 u=0, \quad u(0)=I, \quad u'(0)=0, \quad t \in[0, T] \quad \quad \quad (1.1)$$
with constant amplitude $I$ and angular frequency $\omega$.
The corresponsing period is $P=2\pi/\omega$, and the corresponding frequency is $f=\omega/(2\pi)$.
In this system, $u(t)$ represents displacement, $u'(t)$ velocity, and $u''(t)$ acceleration.
The exact solution of $(1.1)$ is $u(t)=I \cos (\omega t)$.

### 1.1.2 A Centered Finite Difference Scheme
Four steps are to follow to formulate a finite difference method for the problem $(1.1)$.

**Step 1: Discretizing the domain**

The domain is discretized by introducing a uniformly partitioned time mesh.
The points in the mesh are $t_n=n \Delta t, n=0,1, \ldots ,N_t$, where $\Delta t=T/N_t$ is the constant length of the time steps.
A mesh function $u^n$ for $n=0,1, \ldots, N_t$ is introduced to approximate the exact solution $u(t_n)$ at the mesh points:
$$u(t_n) \sim u^n. \quad \quad \quad (1.2)$$
The mesh function $u^n$ will be computed from algebraic equations derived from the differential equation problem. 

**Step 2: Fulfilling the equation at discrete time points**

The ODE is to be satisfied at each mesh point, then the solutions must be found for
$$u''(t_n)+\omega^2 u(t_n)=0, \quad n=0,1, \ldots, N_t. \quad \quad \quad (1.3)$$

**Step 3: Replacing derivatives with finite differences**

By Taylor expansion, 
$$u(t+\Delta t)=u(t)+u'(t) \Delta t +\frac{1}{2} u''(t)(\Delta t)^2+\frac{1}{6} u'''(t)(\Delta t)^3+O((\Delta t)^4),$$
$$u(t-\Delta t)=u(t)-u'(t) \Delta t +\frac{1}{2} u''(t)(\Delta t)^2-\frac{1}{6} u'''(t)(\Delta t)^3+O((\Delta t)^4).$$
Adding these two expansions and solving for $u''(t)$, we have
$$u''(t)=\frac{u(t+\Delta t)-2 u(t)+u(t-\Delta t)}{(\Delta t)^2}+O((\Delta t)^2).$$
Replacing $t$ by $t_n=n \Delta t$ and using $(1.2)$ yields the common second-order accurate centered approximation to the second-order derivative
$$u''(t_n) \sim  \frac{u^{n+1}-2 u^n+u^{n-1}}{\Delta t^2}. \quad \quad \quad (1.4)$$
Inserting $(1.2)$ and $(1.4)$ into $(1.3)$, we have
$$\frac{u^{n+1}-2 u^n+u^{n-1}}{\Delta t^2}+\omega^2 u^n=0, \quad n=0,1, \ldots, N_t. \quad \quad \quad (1.5)$$
Additionally, we need to replace the derivative in the initial condition by a finite difference.
Utilizing the same two Taylor expansions and solving for $u'(t)$ yields
$$u'(t)=\frac{u(t+\Delta t)-u(t-\Delta t)}{2 \Delta t}+O((\Delta t)^2).$$
Substituting $t$ by $t_n=n \Delta t$ and using $(1.2)$ produces the second-order accurate centered approximation to the first-order derivative
$$u'(t_n) \sim \frac{u^{n+1}-u^{n-1}}{2 \Delta t}. \quad \quad \quad (1.6)$$
From the initial conditions $u(0)=I$ and $u'(0)=0$ we have
$$u^0=I \quad \text{and} \quad \frac{u^1-u^{-1}}{2 \Delta t}=0. \quad \quad \quad (1.7)$$

**Step 4: Formulating a recursive algorithm**

To formulate the computational algorithm, we assume that $u^{n-1}$ and $u^n$ have already been calculated, then $u^{n+1}$ is the unknown value solved as
$$u^{n+1}=2 u^n-u^{n-1}-\Delta t^2 \omega^2 u^n. \quad \quad \quad (1.8)$$
The algorithm is to apply $(1.8)$ successively for $n=0,1, \ldots, N_t-1$.
This numerical scheme sometimes goes under the name Störmer’s method, Verlet integration, or the Leapfrog method.
When computing the first step at $n=0$, the discretization of initial conditions $(1.7)$ is needed to resolve the values of $u^{-1}$ and $u^0$.
As a result, $u^1$ is calculated as
$$u^1=I-\frac{1}{2} \Delta t^2 \omega^2 I. \quad \quad \quad (1.9)$$
We introduce a group of compact operator notations for finite differences:
$$u'(t_n) \sim \frac{u^n-u^{n-1}}{\Delta t} \equiv [D_t^{-} u]^n$$
$$u'(t_n) \sim \frac{u^{n+1}-u^{n}}{\Delta t} \equiv [D_t^{+} u]^n$$
$$u'(t_n) \sim \frac{u^{n+\frac{1}{2}}-u^{n-\frac{1}{2}}}{\Delta t} \equiv [D_t u]^n$$
$$u'(t_n) \sim \frac{u^{n+1}-u^{n-1}}{2 \Delta t} \equiv [D_{2t} u]^n$$
$$u''(t_n) \sim \frac{u^{n+1}-2 u^n+u^{n-1}}{\Delta t^2} \equiv  [D_t(D_t u)]^n \equiv [D_t D_t u]^n$$
The operator is built of the symbol $D$, with the independent variable as its subscript and the type of difference as its superscript.
The superscripts $-$ and $+$ indicate backward and forward differences, respectively.
Square brackets are placed around the operator and the function it operates on to specify the mesh point that the finite difference approximates to by a superscript after the closing bracket.
Therefore, the scheme $(1.5)$ can be written as
$$[D_t D_t u+\omega^2 u=0]^n, \quad n=0,1, \ldots, N_t. \quad \quad \quad$$
The discretization of initial conditions can be expressed as
$$[u=I]^0 \quad \text{and} \quad [D_{2 t} u=0]^0.$$
In difference equations the square brackets are often placed around the whole equation to indicate at which mesh point the equation applies,
since each term must be approximated at the same point.

## 1.2 Implementation



