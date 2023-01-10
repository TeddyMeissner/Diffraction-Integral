# Fresnel Diffraction Integral

To describe the propagation of a light beam through a medium one begins with the general wave equation for a Nonlinear Optical Media derived from the Maxwell equations. In this paper we focus on solving the Helmholtz Paraxial equation, a slowly-varying envelope function for monochromatic waves and short-range approximation to the general wave equation. The solution to the Helmholtz Paraxial equation is known as the Fresnel Diffraction approximation. 

Below is a brief overview of the two methods while a more in depth overview can be found in the [writup](Fresnel_Diffraction_writeup.pdf) and [presentation](Fresnel_diffraction_presentation.pptx). 

### Findings: 
In the case of the Fresnel Diffraction Approximation we found the FFT based method introduces artificial periodic boundary conditions and the accuracy of the algorithm depends on propagation distance, wavelength, and observation plane discretization whereas the sinc based method relies on only how well the source field (initial condition) is approximated. 

We are looking to solve the Helmholtz Paraxial equation, 
$$(2ik\partial_z + \nabla_T^2)\psi(\boldsymbol{r}) = 0$$
where the solution is known as the Fresnel Diffraction Approximation or integral, 
$$U(X,Y) = \frac{-ike^{ikz}}{2\pi z}\int_{-\infty}^\infty\int_{-\infty}^\infty e^{\frac{ik}{2z}((X-x)^2 + (Y-y)^2)}u(x,y) dx dy$$

## Spectral Method
    
We first note we can write the Fresnel Diffraction Approximation in terms of a convolution. Let the Fresnel kernel be denoted as, 
$$h_F(x,y) = \frac{-ike^{ikz}}{2\pi z}e^{\frac{ik}{2z}(x^2 + y^2)}$$
where, 
$$U(X,Y) = (h_F\star u)(X,Y) = \frac{-ike^{ikz}}{2\pi z}\int_{-\infty}^\infty\int_{-\infty}^\infty e^{\frac{ik}{2z}((X-x)^2 + (Y-y)^2)}u(x,y) dx dy$$
Now recalling the convolution theorem, 
$$\mathcal{F}(h_F\star u)(\xi,\eta) = \hat{h_F}(\xi,\eta)\hat{u}(\xi,\eta)$$
where, 
$$\hat{f}(\xi,\eta) = \int_{-\infty}^\infty \int_{-\infty}^\infty f(x,y)e^{-i2\pi(x\xi + y\eta)}dx dy$$
We can analytically compute,
$$\hat{h_f}(\xi,\eta) = \frac{-ike^{ikz}}{2\pi z}\int_{-\infty}^\infty \int_{-\infty}^\infty e^{\frac{ik}{2z}(x^2 + y^2)}e^{-i2\pi(x\xi + y\eta)}dx dy$$
$$\hat{h_f}(\xi,\eta) = e^{ikz}\exp\left(\frac{-i2\pi^2z}{k}(\xi^2+\eta^2)\right)$$
Depending on the initial source $u(x,y)$ may need to be computed numerically thus we get, 
$$U(X,Y) = \mathcal{F}^{-1}[\hat{h_F}\hat{u}]$$
$$ = e^{ikz}\int_{-\infty}^\infty \int_{-\infty}^\infty \exp\left(\frac{-i2\pi^2z}{k}(\xi^2+\eta^2)\right)\hat{u}(\xi,\eta) e^{i2\pi(X\xi + Y\eta)}d\xi d\eta$$
As the spectral method approximates the inverse Fourier Transform by the composite trapezoidal rule, $U(X,Y)$ is computed as, 
$$U(X,Y) \approx e^{ikz}\frac{1}{L^2}\sum_{m = -L/2}^{L/2} \sum_{n = -L/2}^{L/2}  \exp\left(\frac{-i2\pi^2z}{k}(\xi_m^2+\eta_n^2)\right)\hat{u}(\xi_m,\eta_n) e^{i2\pi(X\xi_m + Y\eta_n)}d\xi d\eta$$

Where we can see artificial boundary conditions are imposed and waved reaching the boundary are reintroduced rather than dispersed through infinity. 

## Sinc Based Method

Let $f(x)$ be a band limited function (chopping off the initial condition where we want it to be true 0 outside a certain domain). 
Then $f$ can be represented exactly by,
$$f(x) = \sum_{n=-\infty}^\infty f_n \sinc\left( \frac{x - x_n}{\Delta x}\right)$$
The full derivation can be found in Shannon's Sampling Theorem. Thus we can rewrite our initial condition as, 
$$u(x,y) = \sum_{m=-\infty}^\infty\sum_{n=-\infty}^\infty u_{mn} \sinc\left( \frac{x - x_m}{\Delta x}\right) \sinc\left( \frac{y - y_n}{\Delta y}\right)
$$
The solution the Helmholtz Paraxial equation is given by the convolution with the Fresnel kernel $h_f$, 
$$U(X,Y) = \sum_{m=-\infty}^\infty\sum_{n=-\infty}^\infty u_{mn} \left(h_f \star (\sinc\left( \frac{x - x_m}{\Delta x}\right) \sinc\left( \frac{y - y_n}{\Delta y}\right))\right)(X,Y)$$

To generalize the right hand side (let us assume $\delta = \Delta x,\Delta y$) and let, 
$$\Phi(X,Y) = \int_{-\infty}^\infty\int_{-\infty}^\infty h_f(X - x,Y-y)\sinc(\frac{x}{\delta})\sinc(\frac{y}{\delta})dxdy$$

Where we use the convolution property of the Fourier Transform  and the fact that the Transform of a sinc function,
$$\mathcal{F}\left[\sinc(\frac{x}{\delta})\right](\xi) = \delta \rect(\delta \xi )$$

Giving us the inverse Fourier transform vanishes outside of the square $-\frac{1}{2}\leq \xi\delta \leq \frac{1}{2}$,$-\frac{1}{2}\leq \eta\delta \leq \frac{1}{2}$ thus let $W = \frac{1}{2\delta}$ then, 

$$\Phi(X,Y) = \delta^2 \int_{-W}^W\int_{-W}^W \hat{h}_f(\xi,\eta)e^{i2\pi(X\xi + Y\eta)}d\xi d\eta$$

Which we can be accurately computed on a chosen domain. We can now rewrite our solution as, 
$$U(X,Y) = \sum_{n = -\infty}^{\infty}\sum_{m = -\infty}^{\infty} u_{nm} \Phi(X - x_n,Y-y_m)$$
We can see clearly now the Sinc based method never assumes periodic boundary conditions. Thus the numerical solution will only depend on the fineness of the grid rather than the is propagation distance, wavelength, and observation plane discretization   \\
    
A numerically efficient method is derived in [1] which for conciseness of the paper will not be derived but is used in the creation of the numerical results. The method derived was writing the solution in terms of matrix multiplications, 

$$\boldsymbol{U}_{mn} = e^{ikz} \boldsymbol{w}_{mj}^x \boldsymbol{u}_{jl} (\boldsymbol{w}_{nl}^y)^T$$
but hello
$$
\begin{align*}
    \boldsymbol{U}_{mn} &= U(X_m,Y_n) &
    \boldsymbol{w}_{mj}^x &= \phi(X_m - x_j)\\
    \boldsymbol{u}_{jl} &= u(x_j,y_l) & 
    \boldsymbol{w}_{nl}^y &= \phi(Y_n - y_l)
\end{align*}
$$
Where, 
$$\phi(X) = \frac{\delta}{\pi}\sqrt{\frac{k}{2z}}\exp\left(i\frac{X^2k}{2z}\right)\left(C(\mu_2) - C(\mu_1) - iS(\mu_2) - iS(\mu_1))\right) $$
$$
\begin{align*}
    C(x) &= \int_0^x \cos(\mu^2)d\mu &
    S(x) &= \int_0^x \sin(\mu^2)d\mu \\
    \\ \mu_1 &= -\pi\sqrt{\frac{2z}{k}}W - \sqrt{\frac{k}{2z}}X &
    \mu_1 &= \pi\sqrt{\frac{2z}{k}}W - \sqrt{\frac{k}{2z}}X \\
\end{align*}
$$