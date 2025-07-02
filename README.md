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
Usage:

1) Download *pa_nearfield_flexible_core.m* and *pa_nearfield_flexible_wrapper.m*
2) Run pa_nearfield_flexible_wrapper.m
3) Adjust settings on GUI to your liking

Please contact me via lukepeterwilson26@gmail.com for bugs or further inquiry
