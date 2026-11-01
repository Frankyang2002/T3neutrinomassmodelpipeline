# T3 16-model numerical comparison

**Exploratory, not certified physical neutrino predictions.**

14 models: improved provisional direct matching and SM+Weinberg running.
B(-1), C(-1): archived fixed-scale C5 imposed at 1 TeV on a separate improved SM background; mixed-scalar matching and consistency not established.
All values are read from saved outputs; no physics is recalculated.

| model | calculation_method | m1_ev | m2_ev | m3_ev | delta_m21_sq_ev2 | delta_m31_sq_ev2 |
| --- | --- | --- | --- | --- | --- | --- |
| A_-4 | ImprovedDirectMatchingProvisional | 5.1080841 | 5.4789171 | 5.9158598 | 3.9260096 | 8.9048743 |
| A_-2 | ImprovedDirectMatchingProvisional | 12.63485 | 13.552494 | 14.63379 | 24.030652 | 54.508394 |
| A_+0 | ImprovedDirectMatchingProvisional | 12.614917 | 13.53021 | 14.608597 | 23.930447 | 54.274953 |
| B_-3 | ImprovedDirectMatchingProvisional | 4.7470316 | 5.0911114 | 5.4964467 | 3.3851067 | 7.6766173 |
| B_-1 | ConditionalArchivedC5 | 0.0037044234 | 0.0082268928 | 0.041753348 | 5.3959013e-05 | 0.0017296193 |
| B_+1 | ImprovedDirectMatchingProvisional | 4.7706314 | 5.1162592 | 5.5233929 | 3.4171843 | 7.7489447 |
| C_-3 | ImprovedDirectMatchingProvisional | 2.7232727 | 2.9234388 | 3.1596748 | 1.1302801 | 2.5673307 |
| C_-1 | ConditionalArchivedC5 | 0.030874767 | 0.03212983 | 0.051463445 | 7.9074713e-05 | 0.0016952349 |
| C_+1 | ImprovedDirectMatchingProvisional | 2.734674 | 2.9355957 | 3.1727109 | 1.1392802 | 2.5876522 |
| D_-2 | ImprovedDirectMatchingProvisional | 12.726652 | 13.651749 | 14.741949 | 24.402577 | 55.357381 |
| D_+0 | ImprovedDirectMatchingProvisional | 12.732369 | 13.658228 | 14.74938 | 24.433975 | 55.431001 |
| D_+2 | ImprovedDirectMatchingProvisional | 5.1795976 | 5.5560218 | 5.9996142 | 4.0411469 | 9.1671392 |
| E_-4 | ImprovedDirectMatchingProvisional | 8.3400874 | 8.9386721 | 9.6428945 | 10.342802 | 23.428356 |
| E_-2 | ImprovedDirectMatchingProvisional | 18.66256 | 20.006066 | 21.587311 | 51.951512 | 117.72084 |
| E_+0 | ImprovedDirectMatchingProvisional | 18.643592 | 19.985334 | 21.564442 | 51.830047 | 117.44162 |
| E_+2 | ImprovedDirectMatchingProvisional | 8.2840179 | 8.8779553 | 9.5766132 | 10.193137 | 23.086568 |

## Absolute PMNS matrices

Entries are magnitudes |U_ij| (phases are not reported).

### A_-4 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### A_-2 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 1.7736243e-31 | 1.2784221e-16 | 1 |
| 2 | 1.3323434e-15 | 1 | 6.2487696e-17 |
| 3 | 1 | 1.5552836e-15 | 0 |

### A_+0 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### B_-3 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### B_-1 — ConditionalArchivedC5

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0.82093935 | 0.55128045 | 0.14882357 |
| 2 | 0.32149723 | 0.66163103 | 0.67740971 |
| 3 | 0.47190902 | 0.50826592 | 0.72039408 |

### B_+1 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 3.3509255e-31 | 3.669893e-16 | 1 |
| 2 | 1.3372024e-15 | 1 | 3.3461008e-16 |
| 3 | 1 | 1.5612845e-15 | 2.4751829e-31 |

### C_-3 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 1.5515354e-16 | 1 | 0 |
| 3 | 1 | 1.7804942e-16 | 0 |

### C_-1 — ConditionalArchivedC5

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0.82138361 | 0.55011774 | 0.15066332 |
| 2 | 0.31987004 | 0.66297578 | 0.67686504 |
| 3 | 0.4722416 | 0.50777317 | 0.72052362 |

### C_+1 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### D_-2 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 2.0496557e-16 | 1 |
| 2 | 0 | 1 | 1.6663896e-16 |
| 3 | 1 | 0 | 0 |

### D_+0 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### D_+2 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### E_-4 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### E_-2 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### E_+0 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 1 | 0 | 0 |

### E_+2 — ImprovedDirectMatchingProvisional

| i / j | 1 | 2 | 3 |
|---|---:|---:|---:|
| 1 | 0 | 0 | 1 |
| 2 | 1.5879267e-16 | 1 | 0 |
| 3 | 1 | 1.7807875e-16 | 0 |

