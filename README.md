# Moer-zcx
## 简介
基于Moer-lite，侧重于volumetric渲染算法的实现

## Features
- [X] Volumetric Rendering (Regular tracking homegeneous/heterogeneous medium) 
- [X] Transmittance sampling with delta-tracking and transmittance estimation with ratio-tracking
- [X] MajorantGrid with a DDA-like tracker to accelerate the delta/ratio tracking
- [X] Path-space filtering

## Todo
- [ ] Progressive Photon Mapping (代码的bug我找不到)
  - [ ] 考虑用点光源去做测试

## Other focus
- [ ] Cuda-accelerated 