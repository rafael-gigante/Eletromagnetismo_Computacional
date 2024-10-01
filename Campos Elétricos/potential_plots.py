import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def equipotential_lines(V, Lx, Ly):
    """
    Plots the equipotential lines for a given potential V, for a distance
    Lx and Ly of these axis.

    """

    N = V.shape[0]
    
    # Create a grid for plotting
    x = np.linspace(-Lx/2, Lx/2, N)
    y = np.linspace(-Ly/2, Ly/2, N)

    # Plot equipotential lines with color fill
    plt.figure()
    levels = np.linspace(V.min(), V.max(), 500)
    plt.contourf(x, y, V, levels=levels, cmap=cm.jet, alpha=0.8)
    plt.title(r'Equipotential Lines')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    plt.colorbar(label='Potential (V)')
    plt.savefig("Campos Elétricos/figures/equipotential_lines1.png")
    plt.close()

    # Plot equipotential lines without color fill
    plt.figure()
    plt.contour(x, y, V, levels=np.linspace(V.min(), V.max(), 50) ,colors='black')
    plt.title(r'Equipotential Lines')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    plt.savefig("Campos Elétricos/figures/equipotential_lines2.png")
    plt.close()

def potential_surface(V, Lx, Ly):
    """
    Plots the potential surface for a given potential V, for a distance
    Lx and Ly of these axis.

    """

    N = V.shape[0]
    
    # Create a grid for plotting
    x = np.linspace(-Lx/2, Lx/2, N)
    y = np.linspace(-Ly/2, Ly/2, N)
    X, Y = np.meshgrid(x, y)

    # 3D surface plot of the electric potential
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Create the surface plot and keep a reference to the surface object
    surf = ax.plot_surface(X, Y, V, cmap=cm.jet, alpha=0.8)

    # Add color bar
    cbar = fig.colorbar(surf, location="left")
    cbar.set_label(r'Potential V(x,y)')

    ax.set_title(r'Electric Potential $V(x, y)$')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$V$')

    plt.savefig("Campos Elétricos/figures/3d_potential.png")
    plt.close()
    
def electric_field(V, Lx, Ly):
    """
    Plots the electric field for a given potential V, with the potential in the
    background, for a distance Lx and Ly of these axis.

    """

    N = V.shape[0]

    # Create the grid for x and y
    x = np.linspace(-Lx/2, Lx/2, N)
    y = np.linspace(-Ly/2, Ly/2, N)
    X, Y = np.meshgrid(x, y)

    # Step size for the grid
    dx = Lx / (N - 1)
    dy = Ly / (N - 1)

    # Compute the electric field (Ex, Ey) from the potential V using central differences
    Ex = np.zeros_like(V)
    Ey = np.zeros_like(V)

    # Central differences to compute the gradient of V
    Ex[:, 1:-1] = (-(V[:, 2:] - V[:, :-2]) / (2 * dx))  # ∂V/∂x
    Ey[1:-1, :] = (-(V[2:, :] - V[:-2, :]) / (2 * dy))  # ∂V/∂y

    # Plot the potential as a heatmap
    plt.figure()
    levels = np.linspace(V.min(), V.max(), 500)
    plt.contourf(X, Y, V, levels=levels, cmap=cm.jet, alpha=0.8)
    plt.colorbar(label='Potential (V)')

    # Plot the electric field vectors as a quiver plot
    n=2
    plt.quiver(X[::n, ::n], Y[::n, ::n], Ex[::n, ::n], Ey[::n, ::n], color='black')

    plt.title('Electric Field and Potential')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min() , y.max()])
    
    plt.savefig('Campos Elétricos/figures/electric_field.png')
