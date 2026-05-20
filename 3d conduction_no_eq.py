import numpy as np
import matplotlib.pyplot as plt

def solve_3d_conduction():
    print("\n3D Steady State Heat Conduction Solver (Fixed Boundaries)")
      
    # Inputs
    try:
        Lx = float(input("Enter length in x (Lx): "))
        Ly = float(input("Enter length in y (Ly): "))
        Lz = float(input("Enter length in z (Lz): "))
        nx_div = int(input("Enter divisions in x (nx_div): "))
        ny_div = int(input("Enter divisions in y (ny_div): "))
        nz_div = int(input("Enter divisions in z (nz_div): "))
        
        t_left = float(input("Enter x=0 (Left Face) temperature: "))
        t_right = float(input("Enter x=Lx (Right Face) temperature: "))
        t_bottom = float(input("Enter y=0 (Bottom Face) temperature: "))
        t_top = float(input("Enter y=Ly (Top Face) temperature: "))
        t_front = float(input("Enter z=0 (Front Face) temperature: "))
        t_back = float(input("Enter z=Lz (Back Face) temperature: "))
    except ValueError:
        print("Invalid input.")
        return

    nx = nx_div + 1
    ny = ny_div + 1
    nz = nz_div + 1
    total_nodes = nx * ny * nz
    
    dx = Lx / nx_div
    dy = Ly / ny_div
    dz = Lz / nz_div
    
    bx = 1.0 / dx**2
    by = 1.0 / dy**2
    bz = 1.0 / dz**2
    b_center = -2 * (bx + by + bz)
    
    A = np.zeros((total_nodes, total_nodes))
    B = np.zeros(total_nodes)
    def get_idx(i, j, k): return k * (nx * ny) + i * nx + j  # k is z, i is y, j is x
    
    for k in range(nz):
        for i in range(ny):
            for j in range(nx):
                idx = get_idx(i, j, k)
                
                # Check which faces this node belongs to
                faces_temp = []
                if j == 0: faces_temp.append(t_left)
                if j == nx - 1: faces_temp.append(t_right)
                if i == 0: faces_temp.append(t_bottom)
                if i == ny - 1: faces_temp.append(t_top)
                if k == 0: faces_temp.append(t_front)
                if k == nz - 1: faces_temp.append(t_back)
                
                if faces_temp:
                    # Average the temperatures if it's on an edge or corner (combines 2 or 3 faces)
                    bound_temp = sum(faces_temp) / len(faces_temp)
                    A[idx, idx] = 1; B[idx] = bound_temp
                else:
                    # Interior Nodes
                    A[idx, idx] = b_center
                    A[idx, get_idx(i, j-1, k)] = bx
                    A[idx, get_idx(i, j+1, k)] = bx
                    A[idx, get_idx(i-1, j, k)] = by
                    A[idx, get_idx(i+1, j, k)] = by
                    A[idx, get_idx(i, j, k-1)] = bz
                    A[idx, get_idx(i, j, k+1)] = bz
                    B[idx] = 0

    T_flat = np.linalg.solve(A, B)
    T = T_flat.reshape((nz, ny, nx))

    # Calculate mid-z slice
    x_vals, y_vals, z_vals = np.linspace(0, Lx, nx), np.linspace(0, Ly, ny), np.linspace(0, Lz, nz)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    
    print("\nFinal Results (All Nodes by Z-Slice):")
    for k in range(nz):
        print(f"Z layer = {k} (z = {z_vals[k]:.3f} m)")
        for i in range(ny):
            for j in range(nx):
                print(f"Node ({j},{i},{k}): {T[k, i, j]:.2f} K")
        print()

    # Plotting 1: Interactive 3D Volume (Scatter)
    fig2 = plt.figure("Interactive 3D Grid", figsize=(10, 8))
    ax = fig2.add_subplot(111, projection='3d')
    
    X_list, Y_list, Z_list, T_list = [], [], [], []
    for k in range(nz):
        for i in range(ny):
            for j in range(nx):
                x, y, z = x_vals[j], y_vals[i], z_vals[k]
                t_val = T[k, i, j]
                
                X_list.append(x)
                Y_list.append(y)
                Z_list.append(z)
                T_list.append(t_val)
                
                # Add text label for each 3D node (zorder=10 brings text to the very front)
                ax.text(x, y, z, f"({j},{i},{k})\n{t_val:.0f} K", 
                        size=8, color='black', ha='center', va='bottom', zorder=10)

    # Plot the 3D grid points (zorder=1 puts them in the back)
    scatter = ax.scatter(X_list, Y_list, Z_list, c=T_list, cmap='turbo', s=100, alpha=0.9, edgecolors='black', zorder=1)
    fig2.colorbar(scatter, ax=ax, label='Temperature (K)', pad=0.1)

    ax.set_title('3D Temperature Distribution', weight='bold', pad=15, fontname='Times New Roman')
    ax.set_xlabel('Position x (m)', weight='bold', fontname='Times New Roman')
    ax.set_ylabel('Position y (m)', weight='bold', fontname='Times New Roman')
    ax.set_zlabel('Position z (m)', weight='bold', fontname='Times New Roman')

    # Adjust perspective slightly
    ax.view_init(elev=20, azim=45)

    plt.show()

if __name__ == "__main__":
    solve_3d_conduction()
