# 3D Steady-State Heat Conduction Solver (Fixed Boundaries)

## Aim
To solve 3D steady-state heat conduction (Laplace's Equation) in a cuboid with constant temperatures imposed on all six faces (Dirichlet boundary conditions).

## Theory
The governing equation for 3D steady-state conduction is:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} + \frac{\partial^2 T}{\partial z^2} = 0$$

Discretizing the derivatives using second-order central difference schemes:

$$\frac{T_{i,j+1,k} - 2T_{i,j,k} + T_{i,j-1,k}}{\Delta x^2} + \frac{T_{i+1,j,k} - 2T_{i,j,k} + T_{i-1,j,k}}{\Delta y^2} + \frac{T_{i,j,k+1} - 2T_{i,j,k} + T_{i,j,k-1}}{\Delta z^2} = 0$$

For a general node $(i,j,k)$ where $k$ is $z$-index, $i$ is $y$-index, and $j$ is $x$-index, this simplifies to:
$$\beta_x (T_{i,j+1,k} + T_{i,j-1,k}) + \beta_y (T_{i+1,j,k} + T_{i-1,j,k}) + \beta_z (T_{i,j,k+1} + T_{i,j,k-1}) - 2(\beta_x + \beta_y + \beta_z) T_{i,j,k} = 0$$
Where $\beta_x = \frac{1}{\Delta x^2}$, $\beta_y = \frac{1}{\Delta y^2}$, and $\beta_z = \frac{1}{\Delta z^2}$.

Boundary nodes at edges or corners are evaluated by taking the average of their intersecting faces. The solver generates a large sparse system $AT = B$ and displays the temperature values via a 3D scatter grid.

## File Structure
- `3d conduction_no_eq.py` - Core implementation generating the equations, solving the matrix, printing nodal temperatures slice-by-slice, and plotting the interactive 3D scatter.
- `v0.py` - Initial prototype implementation.
- `output.txt` - Output showing temperature readings listed by Z-slices.
- `3D Grid.png` - Scatter visualization of the 3D temperature volume grid.

## How to Run
Ensure you have the required dependencies:
```bash
pip install numpy matplotlib
python "3d conduction_no_eq.py"
```
