"""
Plot SHMIP and KAN_L temperature timeseries
"""

import numpy as np
from matplotlib import pyplot as plt

fig, (ax1, ax2) = plt.subplots(figsize=(8, 4), ncols=2, sharey=True)

# Get sea-level temperature forcing
tt_temp = np.loadtxt('../glads/data/kan_l_melt/KAN_L_2014_temp_clipped.txt', delimiter=',')
tt_days = tt_temp[:, 0]
temp_sl = tt_temp[:, 1]

# SHMIP temp model
DT = 0.0075*390
temp_fun = lambda t: -16*np.cos(2*np.pi*t/86400/365) - 5 + DT

shmip_adj_fun = lambda t: -9.0684*np.cos(2*np.pi*t/86400/365) + 9.0684*(DT - 5)/16

tt = np.arange(0, 365*86400, 86400)
tt_a = tt/365/86400

KAN_temp = np.interp(tt, tt_days*86400, temp_sl, left=0, right=0)
KAN_temp[KAN_temp<0] = 0

SHMIP_temp = temp_fun(tt)
SHMIP_temp[SHMIP_temp<0] = 0

SHMIP_adj_temp = shmip_adj_fun(tt)
SHMIP_adj_temp[SHMIP_adj_temp<0] = 0
ax1.plot(tt_a, SHMIP_temp, label='SHMIP')
ax1.plot(tt_a, SHMIP_adj_temp, label='SHMIP adjusted')
ax1.grid()

ax2.plot(tt_a, KAN_temp)
ax2.grid()

ax1.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticks([0, 0.25, 0.5, 0.75, 1])

ax1.set_ylim([-1, 16])
ax2.set_ylim([-1, 16])

ax1.set_title('SHMIP')
ax2.set_title('KAN_L')

ax1.set_ylabel('Temp ($^\circ$C)')
ax1.set_xlabel('Time (a)')
ax2.set_xlabel('Time (a)')

fig.tight_layout()
fig.savefig('figures/aux/melt_forcing.png', dpi=600)

plt.show()
