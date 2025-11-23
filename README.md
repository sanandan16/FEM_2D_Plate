# 2D Finite Element Solver for Steady-State Heat Transfer on a Flat Plate

This repository implements a **2D Finite Element Method (FEM)** solver to compute the **steady-state temperature distribution** on a flat plate. The project demonstrates the complete FEM pipeline—mesh construction, stiffness matrix assembly, application of mixed boundary conditions, and numerical solution of the heat conduction equation.

This project forms part of my finite element methods coursework.
---

## 🔥 Problem Description

The solver computes the temperature field \( T(x, y) \) on a 2D plate by solving the steady heat conduction equation:

\[
\nabla \cdot (k \nabla T) = 0
\]

with a mix of boundary conditions:

- **Dirichlet BCs** (prescribed temperature)
- **Neumann BCs** (prescribed heat flux)


The implementation uses **linear quadrilateral elements**, though the structure allows extension to higher-order elements.

---

## ✨ Key Features

- Full FEM pipeline written from scratch
- Linear quadrilateral element stiffness assembly
- Multiple boundary conditions (Dirichlet, Neumann)
- Modular code for easy extension
- Contour and surface visualizations of the temperature field
- Example scripts demonstrating typical thermal boundary conditions

---

## 🧠 Numerical Approach

1. **Domain Discretization**
   - Flat Plate
   - Node/element connectivity generation

2. **Stiffness Matrix Assembly**
   - Element-level conductivity matrix:
     \[
     k_e = \int B^T k B \, dA
     \]
   - Assembly into global matrix \( K \)

3. **Boundary Condition Application**
   - Temperature (Dirichlet) nodes fixed
   - Flux boundaries incorporated via load vector

4. **Solve Linear System**
   \[
   K \mathbf{T} = \mathbf{F}
   \]

5. **Post-Processing**
   - Temperature contour maps
   - Mesh visualization
   - Node-level interpolation

---

## ▶️ How to Run

### **MATLAB**
Open MATLAB and run:

```matlab
run('src/fem2d_plate_main.m')
