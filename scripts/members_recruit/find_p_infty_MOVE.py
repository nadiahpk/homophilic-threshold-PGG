
# find p_infty
# ---

# agrees with pdf writeup
pinfty_root_eqn = lambda p: X - Z + (W-X) * binom(n-1,tau-1) * p**(tau-1) + sum( binom(n-1,l) * p**l * (-1)**(l-tau) * ( binom(l-1, tau-2) * (X-W) + binom(l-1, tau-1) * (Z-Y)) for l in range(tau, n) ) 

# TODO - put the p_peak stuff here

# solve the two roots for pinfty
# ---

pinf1 = bisect( pinfty_root_eqn, 1e-6, p_peak )
pinf2 = bisect( pinfty_root_eqn, p_peak, 1-1e-6)


# plot where the p_infty lines are

plt.axhline(pinf1, ls='dotted', color='black')
plt.text(1, pinf1, r'$p_\infty$', va='center', ha='left', fontsize='large')

plt.axhline(pinf2, ls='dotted', color='black')
plt.text(1, pinf2, r'$p_\infty$', va='center', ha='left', fontsize='large')
