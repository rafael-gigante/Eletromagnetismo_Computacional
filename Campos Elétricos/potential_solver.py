import numpy as np

def update_V(V, ρ):
    """
    Function that receives a potential V and aplies the Jacobi method to it
    and calculates the absolute difference between the new and old potentials.

    Returns:
    - `V_new`: Array with the new potential calculate by the Jacobi method,
    - 'ΔV': Float of the absolute difference between the potential received and the calculated.
    """
        
    # Total (cumulative) change in V during this update loop
    ΔV = 0

    # Dimension of the potential
    dim = len(V.shape)

    # Number of sites along a direction
    N = V.shape[0]

    # Step size
    dx = 2.0 / (N - 1)
    
    # Loop through all values (i,j) except the boundary, where the values are fixed
    if dim == 2:
        V_new = np.copy(V)
        for i in range(1, N-1):
            for j in range(1, N-1):
                if (abs(V[i, j]) != 1):
                    V_new[i, j] = (V[i+1, j] + V[i-1, j] +
                                V[i, j+1] + V[i, j-1] + ρ[i, j]*(dx**2)) / 4
                    ΔV += abs(V_new[i, j] - V[i, j])
    elif dim == 3:
        V_new = np.copy(V)
        for i in range(1, N-1):
            for j in range(1, N-1):
                for k in range(1, N-1):
                    V_new[i, j, k] = (V[i+1, j, k] + V[i-1, j, k] +
                                V[i, j+1, k] + V[i, j-1, k] +
                                V[i, j, k+1] + V[i, j, k-1] + ρ[i, j, k]*(dx**2)) / 6
                    ΔV += abs(V_new[i, j, k] - V[i, j, k])
    return V_new, ΔV


def laplace_calculate(V, ρ):
    """
    Function that aplies the relaxation method to a initial potential V
    until the new potential is within the tolerance range.

    Returns:
    - `V`: Array with final potential that is within the tolerance range
    """

    # Define the error tolerance
    N = V.shape[0]
    ε = 10**(-6) * N**2

    # Perform 10 iterations initially
    for _ in range(10):
        V_new, ΔV = update_V(V, ρ)
        V, ΔV = update_V(V_new, ρ)

    # Continue iterating until ΔV is less than ε
    while ΔV < ε:
        V_new, ΔV = update_V(V, ρ)
        V, ΔV = update_V(V_new, ρ)

    return V
