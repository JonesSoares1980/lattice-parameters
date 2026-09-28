Triclinic System Implementation via Linearized Reciprocal Metric Tensors

Overview
The triclinic crystal system features the lowest symmetry in crystallography, governed by six independent geometric variables: three axial lengths ((a), (b), (c)) and three interaxial angles ((alpha), (beta), (gamma)). Direct combinatorial brute-force calculation using standard trigonometric equations introduces severe non-linearity, which leads to numerical divergence and high computational overhead. 

To overcome this bottleneck, this routine introduces a deterministic algebraic linearization by mapping the discrete experimental Bragg reflections ((hkl) onto a Reciprocal Metric Tensor (G*) framework.

Mathematical Foundation

The relationship between the interplanar spacing (d{hkl}) and the reciprocal lattice constants in a triclinic cell is governed by the quadratic form:

1/d{hkl}^2 = h^2 Q{11} + k^2 Q{22} + l^2 Q{33} + 2hk Q{12} + 2hl Q{13} + 2kl Q{23}

Where the parameters (Q{ij}) represent the independent elements of the reciprocal metric tensor (G*):
* (x1 = Q{11})
* (x2 = Q{22})
* (x3 = Q{33})
* (x4 = 2Q{12})
* (x5 = 2Q{13})
* (x6 = 2Q{23})

Since there are 6 unknown variables, the algorithm partitions the pre-indexed data array using the binomial coefficient C(n,6). For a selected batch of (n = 15) peaks, the routine evaluates exactly 5,005 independent linear sub-matrices:

C(15,6) = 15!/(6!(15-6)!) = 5005

Algorithmic Workflow

The `triclinica()` function executes through four fundamental modular stages:

1. Bragg Conversion: Converts raw experimental diffraction angles (2theta) into interplanar distances (d{hkl}) via Bragg's Law, assigning them to their respective (hkl) plane vectors.
2. Combinatorial Combinations (C(n,6)): Generates all possible 6-by-6 combinations of equations using Python's `itertools.combinations`.
3. Linear Matrix Solving: For each combination, a linear equation system of the form (Ax = B) is constructed. The algorithm evaluates the matrix determinant via SymPy/NumPy to prevent linear dependence (discarding ill-conditioned or coplanar reflection subsets where (det(A) approx 0).
4. Tensor Inversion & Metric Extraction: 
   Once the reciprocal matrix \(G^*\) is resolved, the direct metric tensor (G) is obtained via matrix inversion:
   
   [G = (G*)^(-1)]
   
   The direct cell constants are then computed deterministically from the components of (G):
   
   a = sqrt{G{11}}, quad b = sqrt{G{22}}, quad c = sqrt{G{33}}
   [alpha = arccos\left(\frac{G_{23}}{b \cdot c}\right), \quad \beta = \arccos\left(\frac{G_{13}}{a \cdot c}\right), \quad \gamma = \arccos\left(\frac{G_{12}}{a \cdot b}\right)\]

5. **Statistical Aggregation:** Solutions that violate physical constraints (e.g., \(\vert{}\cos\vert{} \ge 1\) or negative tensor diagonals) are filtered out. The remaining successful solutions are aggregated to output the final mean values and the internal geometric consistency index (\(\sigma\)) via `np.std()`.
