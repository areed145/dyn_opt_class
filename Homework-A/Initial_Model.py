from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

m = GEKKO(remote=False)

m.time = np.linspace(0,600,601)

K = m.FV(value=5/8,lb=0,ub=1)
tau = m.Param(value=120)

Qd = np.zeros(601)
Qd[10:] = 80
Qd[200:] = 0
Q = m.Param(value=Qd)

T_ss = 23 # degC
Q_ss = 0  # %

T = m.Var(value=T_ss)

m.Equation(tau*T.dt()==-(T-T_ss)+K*(Q-Q_ss))

m.options.IMODE=4
m.solve()

plt.plot(m.time,Q.value)
plt.plot(m.time,T.value)

plt.show()
