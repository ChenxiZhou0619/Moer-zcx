# Moer-lite：面向教学的蒙特卡洛路径追踪渲染框架

## 简介

Based on Moer-lite

## Features
- [X] Volumetric Rendering (Regular tracking homegeneous/heterogeneous medium) 
- [X] Transmittance sampling with delta-tracking and transmittance estimation with ratio-tracking
- [X] MajorantGrid with a DDA-like tracker to accelerate the delta/ratio tracking

## Todo
- [ ] Path-space filtering
  - [X] KD-Tree for neighbor points query (using nanoflann)