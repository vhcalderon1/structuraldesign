# Structural Engineering Analysis Tools

_A collection of MATLAB scripts for structural analysis, design, and visualization, tailored for academic research and professional engineering applications._

<P align = "center">
   <img src="https://github.com/user-attachments/assets/b936b220-2c94-4efa-9e55-9f3e9ee604e6" width="500" height="350" />
   <img src="https://github.com/user-attachments/assets/816750b8-1b05-4da3-8638-70a00f6f0412" width="500" height="350" />

## References

Calderon, V. H. (2021). “Feasibility of protecting with seismic isolation a limited-ductility-wall social housing.” (In Spanish) Pontifical Catholic University of Peru Digital Repository. Available at: https://tesis.pucp.edu.pe/items/033969e7-0b01-494f-ac01-6d2167bd67d4

Originally created by: Victor Calderon (July - November 2019, April - November 2020, and September 2022)

Updated by: Jefferson de la Cuba (February 2025)

## What is this repository?

This repository contains eleven MATLAB scripts designed for structural engineering analysis, including:
Material behavior modeling: Bilinear stress-strain curves for steel and concrete.
Interaction diagrams: Axial load-moment (P-M) capacity for beams, walls, and T-sections.
Shear force analysis: Comparison of fixed-base vs. isolated-base structural systems.
Moment-curvature relationships: Nonlinear behavior of reinforced concrete sections.
Applied load validation: Overlay of design demands on interaction diagrams.
These tools are developed with academic rigor and industry standards (e.g., ACI 318) in mind, making them suitable for:
Research: Validating theoretical models.
Education: Teaching structural mechanics and reinforced concrete design.
Professional use: Preliminary design and code-compliant analysis.

## How does it work?

### Modular Workflow:
Each script is self-contained and annotated with input parameters (e.g., material properties, geometry).
Uses relative paths (../datasets/, ../outputs/) for cross-platform compatibility.

### Dynamic Analysis:
Iterative neutral-axis depth calculations for interaction diagrams.
Strain-dependent strength reduction factors (φ) per ACI guidelines.

### Visualization:
Publication-ready plots (Times New Roman fonts, axis formatting, gridlines).
Comparison of nominal vs. factored capacities and applied loads.

## Why use this repository?

### For Academic:
Transparent Algorithms: Fully documented code for peer review and reproducibility.
Teaching Resource: Demonstrate concepts like plastic hinge formation or shear distribution.

### For Engineers:
Code Compliance: Embedding ACI 318 assumptions (e.g., β₁ = 0.85, ε_cu = 0.003).
Rapid Prototyping: Adjust inputs (e.g., fpc, fy, reinforcement layout) to test design alternatives.

## License

[MIT](./LICENSE)