import numpy as np

def initialize_V_box(V1, V2, N):
    """
    Initialize an electric potential of a box with two opposite sides having potentials V1 and V2, 
    and the other sides are nonconducting, assuming the potential varies linearly on those sides. 
    The potential grid has dimensions N x N.

    Returns:
    - `V`: Array with the proper boundary for the specified problem
    - `ρ`: Array with the charge distribution for the specified problem
    """
    # Create an initial guess for the potential (N x N grid)
    V = np.zeros((N, N))
    ρ = np.zeros((N, N))

    # Boundary conditions:
    # For the sides with potentials
    for j in range(N):
        V[j, 0] = V1  # Left boundary
        V[j, N-1] = V2  # Right boundary

    # For nonconducting sides
    # Set linear variation between the potential sides
    ΔV = abs(V2 - V1)
    dV = ΔV / (N-1)
    for i in range(1, N-1):
        if V1 < V2:
            V[0, i] = V[0, i-1] + dV
            V[N-1, i] = V[N-1, i-1] + dV
        else:
            V[0, i] = V[0, i-1] - dV
            V[N-1, i] = V[N-1, i-1] - dV

    return V, ρ

def initialize_V_capacitor(V1, V2, N):
    """
    Initialize an electric potential of a capacitor with the plates having potentials V1 and V2.
    The potential grid has dimensions N x N.

    Returns:
    - `V`: Array with the proper boundary for the specified problem
    - `ρ`: Array with the charge distribution for the specified problem
    """
    # Create an initial guess for the potential (N x N grid)
    V = np.zeros((N, N))
    ρ = np.zeros((N, N))

    plate_begin = int(3*(N-1)/10)
    plate_end = int(7*(N-1)/10)
    first_third = int(N/3)
    second_third = int(2*N/3)
    # Position of the plates
    V[plate_begin:plate_end, first_third] = V1
    V[plate_begin:plate_end, second_third] = V2

    return V, ρ

def initialize_V_box_charge(N):
    """
    Initialize an electric potential of a cubic box and a point charge in the center.
    The potential grid and the charge grid has dimensions N x N x N.

    Returns:
    - `V`: Array with the proper boundary for the specified problem
    - `ρ`: Array with the charge distribution for the specified problem
    """
    # Create an initial guess for the potential (N x N x N grid)
    V = np.zeros((N, N, N))
    ρ = np.zeros((N, N, N))

    zero = int((N-1)/2) + 1
    ρ[zero, zero, zero] = 1
    return V, ρ