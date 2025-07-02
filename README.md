# near_field_simulator

This interactive MATLAB GUI lets you explore the near-field MIMO performance of custom phased-array geometries. You can sweep:

- **Vertical Offset** (distance between TX and RX planes)  
- **Horizontal Offset** (lateral shift between arrays)  
- **Number of TX/RX sub-arrays**  
- **Array spreading**  

and instantly view four heatmaps showing:  

1. **Spectral efficiency** (Capacity in bits/s/Hz)  
2. **Eigenvalue ratio** (dB)  
3. **Principal eigenvalue** (dB)  
4. **Secondary eigenvalue** (dB)  

Under the hood, it builds the near-field channel matrix from element positions, performs SVD, applies water-filling power allocation, and computes throughput.

---
