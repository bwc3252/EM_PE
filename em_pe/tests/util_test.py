from __future__ import print_function
import em_pe

f = em_pe.utils.precompute_redshift()
print("redshift at 40.0 Mpc: z =", f(40.0))
