# Data Inventory

This data is from epithelial monolayer inflation experiments using the **Monolayer Inflator (MOLI)** device (see Figures 3 & 4 of Ouzeri et al., 2025).

---

### 1. ConstantPressureData.dat

This file contains "Pressure-Controlled Ensemble" data (Figure 3d-g). It tracks the geometric of multiple domes subjected to a constant step pressure (e.g., 200 Pa) as they expand and approach a mechanical steady state.

**Columns:**
| Column | Description |
| :--- | :--- |
| `time` | Relative time from pressure application. |
| `p` | Applied hydrostatic pressure (Pa). |
| `height` | Apex height of the dome in microns. |
| `baseradius` | Radius of the micropatterned base in microns. |
| `curvatureRadius` | Calculated radius of curvature in microns. |
| `strain` | Dimensionless strain calculated from geometry. |
| `tension` | Resulting wall tension in mN/m. |
| `SourceFile` | Name of the raw data source file. |

---

### 2. ConstitutiveRelationData.dat

This dataset contains the steady-state mechanical response of the tissue (Figure 3j).

**Columns:**
| Column | Description |
| :--- | :--- |
| `Trial_ID` | Numerical identifier for each experimental repeat. |
| `Time_s` | Relative time in seconds. |
| `Strain` | Dimensionless tissue strain. |
| `Tension_mNm` | Tissue surface tension in mN/m. |
---

### 3. HystersisData.dat
This dataset corresponds to the cyclic loading experiments (Figure 4). Domes were subjected to triangular pressure waves at three distinct periods (20s, 266s, and 2000s). The resulting tension-strain loops reveal the rate-dependent visco-elastic nature and hysteresis of the epithelial tissue.

**Columns:**
| Column | Description |
| :--- | :--- |
| `Time` | Experimental time in seconds. |
| `Pressure` | Dynamic pressure applied during the cycle in Pa |
| `Strain` | Resulting tissue strain. |
| `Tension` | Resulting tissue tension in mN/m. |
| `Height` | Apex height during the cycle in microns. |
| `Base` | Radius of the dome base in microns. |
| `Radius` | Radius of curvature in microns. |
| `Period` | Loading rate label (**20s**, **266s**, or **2000s**). |
| `Dome_ID` | Identifier for specific domes (1–5) within each period group. |