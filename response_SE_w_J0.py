# matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import math as m

from qutip import *

# parameter definition
gamma = 1        # spontaneous emission rate of atoms
yij = 5*gamma        # atom dissipation rate induced by waveguide
kappa = 20*gamma      # cavity dissipation rate
N = 3            # number of cavity fock states
gc = 20*gamma        # QE-cavity coupling
p1 = 0               # phase factor of the waveguide-cavity junction to the first atom in waveguide
p2 = np.pi * 1.5         # phase factor between the atoms in waveguide
J0 = 8*yij            # interaction strength for dimerization interaction
p = 0.3*np.pi          # tunable parameter of dimerization interaction
dw = 0             #  atom-cavity detuning
dw0 = 0             #  atom-cavity detuning

Jp = J0*(1+np.cos(p)) 
Jn = J0*(1-np.cos(p))   

# intial state
psi0 = tensor(basis(N,0), basis(N,0), basis(2,0), basis(2,0), basis(2,0), basis(2,1));    # start with an excited atom

# operators
cR  = tensor(destroy(N), qeye(N), qeye(2), qeye(2), qeye(2), qeye(2));
cL  = tensor(qeye(N), destroy(N), qeye(2), qeye(2), qeye(2), qeye(2));
sm1 = tensor(qeye(N), qeye(N), destroy(2), qeye(2), qeye(2), qeye(2));
sm2 = tensor(qeye(N), qeye(N), qeye(2), destroy(2), qeye(2), qeye(2));
sm3 = tensor(qeye(N), qeye(N), qeye(2), qeye(2), destroy(2), qeye(2));
sm0 = tensor(qeye(N), qeye(N), qeye(2), qeye(2), qeye(2), destroy(2));



expca1 = complex( np.cos(p1), np.sin(p1) )     # phase factor for cavity-atom coupling in exp form
expca2 = complex( np.cos(p1+p2), np.sin(p1+p2) )
expca3 = complex( np.cos(p1+2*p2), np.sin(p1+2*p2) )
expca1c = complex( np.cos(-p1), np.sin(-p1) )    # conjugate
expca2c = complex( np.cos(-p1-p2), np.sin(-p1-p2) )
expca3c = complex( np.cos(-p1-2*p2), np.sin(-p1-2*p2) )

expaa1 = complex( np.cos(p2), np.sin(p2) )     # phase factor for atom-atom coupling in exp form
expaa2 = complex( np.cos(2*p2), np.sin(2*p2) )         
expaa1c = complex( np.cos(-p2), np.sin(-p2) )    # conjugate
expaa2c = complex( np.cos(-2*p2), np.sin(-2*p2) )

Laa123 = yij * expaa1 * ( spre(sm1) * spost(sm2.dag()) - spre(sm2.dag() * sm1)) \
    + yij * expaa1c * ( spre(sm2) * spost(sm1.dag()) - spost(sm1.dag() * sm2)) \
        + yij * expaa2 * ( spre(sm1) * spost(sm3.dag()) - spre(sm3.dag() * sm1)) \
             + yij * expaa2c * ( spre(sm3) * spost(sm1.dag()) - spost(sm1.dag() * sm3));    # the second line in Eq.(2)
Laa321 = yij * expaa1 * ( spre(sm3) * spost(sm2.dag()) - spre(sm2.dag() * sm3)) \
    + yij * expaa1c * ( spre(sm2) * spost(sm3.dag()) - spost(sm3.dag() * sm2)) \
        + yij * expaa2 * ( spre(sm3) * spost(sm1.dag()) - spre(sm1.dag() * sm3)) \
             + yij * expaa2c * ( spre(sm1) * spost(sm3.dag()) - spost(sm3.dag() * sm1));    # the second line in Eq.(2)
Laa23 = yij * expaa1 * ( spre(sm2) * spost(sm3.dag()) - spre(sm3.dag() * sm2)) \
    + yij * expaa1c * ( spre(sm3) * spost(sm2.dag()) - spost(sm2.dag() * sm3));    # the second line in Eq.(2)
Laa21 = yij * expaa1 * ( spre(sm2) * spost(sm1.dag()) - spre(sm1.dag() * sm2)) \
    + yij * expaa1c * ( spre(sm1) * spost(sm2.dag()) - spost(sm2.dag() * sm1));    # the second line in Eq.(2)
LacR = np.sqrt(kappa*yij) * expca1 * ( spre(cR) * spost(sm1.dag()) - spre(sm1.dag() * cR)) \
    + np.sqrt(kappa*yij) * expca1c * ( spre(sm1) * spost(cR.dag()) - spost(cR.dag() * sm1)) \
        + np.sqrt(kappa*yij) * expca2 * ( spre(cR) * spost(sm2.dag()) - spre(sm2.dag() * cR)) \
             + np.sqrt(kappa*yij) * expca2c * ( spre(sm2) * spost(cR.dag()) - spost(cR.dag() * sm2)) \
                  + np.sqrt(kappa*yij) * expca3 * ( spre(cR) * spost(sm3.dag()) - spre(sm3.dag() * cR)) \
                       + np.sqrt(kappa*yij) * expca3c * ( spre(sm3) * spost(cR.dag()) - spost(cR.dag() * sm3));    # the third line in Eq.(2)
LacL = np.sqrt(kappa*yij) * expca1 * ( spre(sm1) * spost(cL.dag()) - spre(cL.dag() * sm1)) \
    + np.sqrt(kappa*yij) * expca1c * ( spre(cL) * spost(sm1.dag()) - spost(sm1.dag() * cL)) \
         + np.sqrt(kappa*yij) * expca2 * ( spre(sm2) * spost(cL.dag()) - spre(cL.dag() * sm2)) \
              + np.sqrt(kappa*yij) * expca2c * ( spre(cL) * spost(sm2.dag()) - spost(sm2.dag() * cL)) \
                   + np.sqrt(kappa*yij) * expca3 * ( spre(sm3) * spost(cL.dag()) - spre(cL.dag() * sm3)) \
                        + np.sqrt(kappa*yij) * expca3c * ( spre(cL) * spost(sm3.dag()) - spost(sm3.dag() * cL));    # the fourth line in Eq.(2)
La1 = 2*yij * ( spre(sm1) * spost(sm1.dag()) - 0.5 * spost(sm1.dag() * sm1) - 0.5 * spre(sm1.dag() * sm1));  # atom decay
La2 = 2*yij * ( spre(sm2) * spost(sm2.dag()) - 0.5 * spost(sm2.dag() * sm2) - 0.5 * spre(sm2.dag() * sm2));
La3 = 2*yij * ( spre(sm3) * spost(sm3.dag()) - 0.5 * spost(sm3.dag() * sm3) - 0.5 * spre(sm3.dag() * sm3));
LcR = kappa * ( spre(cR) * spost(cR.dag()) - 0.5 * spost(cR.dag() * cR) - 0.5 * spre(cR.dag() * cR));  # cR decay
LcL = kappa * ( spre(cL) * spost(cL.dag()) - 0.5 * spost(cL.dag() * cL) - 0.5 * spre(cL.dag() * cL));  # cL decay


H = dw * cR.dag() * cR + dw * cL.dag() * cL + (dw + dw0) * sm1.dag() * sm1 + (dw + dw0) * sm2.dag() * sm2 \
    + (dw + dw0) * sm3.dag() * sm3 + gc * (cR.dag() * sm0 + cR * sm0.dag()) + gc * (cL.dag() * sm0 + cL * sm0.dag()) \
         + Jn*( sm1.dag()*sm2 + sm2.dag()*sm1 ) + Jp*( sm2.dag()*sm3 + sm3.dag()*sm2 )
L0 = liouvillian(H);
L = L0 + La1 + La2 + La3 + LcR + LcL + LacR + LacL + Laa123 + Laa321 + Laa23 + Laa21;    
   
# collapse operator that describes dissipation
c_ops = [ np.sqrt(gamma)*sm0, np.sqrt(gamma)*sm1, np.sqrt(gamma)*sm2, np.sqrt(gamma)*sm3 ]  # represents spontaneous emission

tlist = np.linspace(0, 0.5, 5000)

options = Options()
# options.nsteps = 20000
# options.max_step = (tlist[1]-tlist[0])/2

n = [sm0.dag()*sm0, cR.dag()*cR, cL.dag()*cL, sm1.dag()*sm1, sm2.dag()*sm2, sm3.dag()*sm3]     # operators
 
ms = mesolve(L, psi0, tlist, c_ops, n, options=options)

n0 = ms.expect[0]
ncR = ms.expect[1]
ncL = ms.expect[2]
n1 = ms.expect[3]
n2 = ms.expect[4]
n3 = ms.expect[5]



fig, ax = plt.subplots(figsize=(8,5))
ax.plot(tlist, n0, 'k', label="emitter")
ax.plot(tlist, ncR, 'b', label="CCW")
ax.plot(tlist, ncL, 'r', label="CW")
ax.legend()
# ax.set_xlim(0, tmax)
# ax.set_ylim(0, 1)
ax.set_xlabel('Time, $t$')
ax.set_ylabel('Occupation')

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(tlist, n1, label="atom 1")
ax.plot(tlist, n2, label="atom 2")
ax.plot(tlist, n3, label="atom 3")
ax.legend()
# ax.set_xlim(0, tmax)
# ax.set_ylim(0, 1)
ax.set_xlabel('Time, $t$')
ax.set_ylabel('Population')