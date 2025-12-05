
The method we try to use to calculate the crosssection is the born approximation.

The relevant chapters of the book *Quantum Mechanics, Volume 2*
are: 
Complement C<sub>VIII</sub> 
Chapter VIII, An elementary approach to the quantum theory of scattering by a potential
as well as Appendix I




Here are the most relevant formulas as screenshots:

Formula 52 is very relevant, as it is a general result of how the fourier transform of a spherially symmetric potential (or wavefunction $\Psi$) looks like:
![[ft of spherically symmetric.png]]


The next formula is the scattering crosssection in the Born approximation:
![[born approx.png]]

for a spherically symmetric potential V(r)


in the book, this was presented for the yukawa potential:
![[born integral yukawa.png]]

recognizing the yukawa potential

![[yukawa.png]]
in the formula above, we can replace by $V(r)$ 
to get the formula for our spherically symmetric potential $V(r)$
(also cancelling the 4$\pi$)

$f = -2 \mu/\hbar^2 1/K \int_0^\infty dr \, r \sin(K r) V(r)$

where:
![[K.png]]

and k = k0
![[group velocity.png]]



the crosssection $\sigma$ is then obtained by:
![[sigma.png]]
### yukawa potential formulas for comparing to numerical solution

yukawa potential solution
![[sigma yukawa.png]]

yukawa potential solution, crosssection at $\theta = 0$ 
![[yukawa potential solution.png]]

