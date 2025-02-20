import matrixfuncs as mf
import numpy as np
import matplotlib.pyplot as plt


def est_M(sampled_vec):
  '''Estimate transition matrix M from sampled function values.
  
  Parameters
  ----------
  sampled_vec : ndarray of shape (dim, N)
    A matrix where each column represents function values over a time window.

  Returns
  -------
  M : ndarray of shape (dim, dim)
    The estimated transition matrix that best maps sampled_vec[:, :-1] to sampled_vec[:, 1:].
  '''
  dim = sampled_vec.shape[0]
  r = range(sampled_vec.shape[1] - dim)
  Ms = [np.linalg.solve(sampled_vec.T[k:k+dim], sampled_vec.T[k+1:k+1+dim]).T for k in r]
  M = np.mean(Ms, 0)
  return M


# Number of frequency components
num = 4

# Generate random function coefficients and frequencies
f_coeffs = np.random.randn(num)
f_freqs = 2*np.pi/(.1 + .9*np.random.rand(num))

# Define the reference function as a weighted sum of sines and cosines
f_ref = lambda t: f_coeffs @ np.concat([np.sin(f_freqs[:num//2]*t), np.cos(f_freqs[num//2:]*t)])

# Sample function values over a finite time range
t_sampled = np.linspace(0, 1, 5*num+1)
f_sampled = np.array(list(map(f_ref, t_sampled)))

# Construct sampled windows for transition matrix estimation
f_sampled_window = np.array([f_sampled[np.arange(i, i+2*num)] for i in range(len(t_sampled) - 2*num)]).T
M = est_M(f_sampled_window)
assert np.allclose(M @ f_sampled_window[:, :-1], f_sampled_window[:, 1:])


# Define evaluation times for function continuation
ts = np.linspace(-1, 2, 4000)
dt_sample = t_sampled[1] - t_sampled[0]
assert np.allclose(t_sampled[1:], t_sampled[:-1] + dt_sample)

# Define the function to be applied to the matrix
f = lambda x: x ** (ts / dt_sample)

# Convert the transition matrix into a functional array representation
arr = mf.FArray.from_matrix(M)

# Define initial left-hand side and right-hand side vectors
v0_lhs = (np.arange(f_sampled_window.shape[0])==0)*1.
v0_rhs = f_sampled_window[:, 0]

# Compute the function applied to the matrix and evaluate it
vals = v0_lhs @ arr @ v0_rhs  # Compute matrix function application
f_M = vals(f)  # Evaluate the function over the time steps


# Plot the computed function values
fig, (ax, bx) = plt.subplots(2, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
fig.suptitle(f'Smooth Continuation of Random Signal with DOF={2*num}')
fig.tight_layout(pad=1.)

# Plot the reconstructed function
ax.plot(ts, f_M, 'b-', label='Continuation of sampled function')
ax.plot(t_sampled, f_sampled, 'ro', label=f'Sampled at step size {dt_sample}')
ax.set_xlim(-1, 2)
ax.set_ylabel('Signal')
ax.legend()
ax.grid(True)

# Compute error between estimated function and reference function
err = f_M - np.array([f_ref(t) for t in ts])

# Plot the error
bx.plot(ts, err, 'b-', label='Error')
bx.plot(ts, np.zeros_like(ts), 'g-')
bx.set_ylim(np.array([-1.05, 1.05]) * np.max(abs(err)))
bx.set_xlabel('Time ($t$)')
bx.set_ylabel('Signal Error')
bx.legend()
bx.grid(True)

# Display and save the plot
plt.show(fig)
fig.savefig('many_frequencies.svg', dpi=200, bbox_inches='tight')
