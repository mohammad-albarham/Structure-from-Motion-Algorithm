[![View Structure-from-Motion-Algorithm on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/179169-structure-from-motion-algorithm)

# Structure from Motion (SfM) Pipeline

The `run_sfm` function implements a **Structure from Motion (SfM)** pipeline for reconstructing 3D structures from image datasets. It supports multiple options for visualization, logging, and dataset control, making it suitable for educational, research, and development purposes.

---

## Function Overview

The `run_sfm` function processes a specific dataset and performs the following tasks:
1. **Feature extraction** using SIFT.
2. **Estimation of relative camera poses** using robust methods like RANSAC.
3. **Triangulation** of 3D points.
4. **Optimization and reconstruction** of absolute camera poses and 3D points.
5. **Visualization of 3D structures and camera poses.

---

## Function Signature

```matlab
function [] = run_sfm(dataset_num, varargin)
```
---

# How to Run the Code

Follow the steps below to successfully execute the `run_sfm` function and process a dataset for Structure from Motion (SfM):

---

## 1. **Prerequisites**
Before running the code, ensure you have the following:
- **MATLAB**: Installed and properly set up.
- **VLFeat Toolbox**: 
  - Download the toolbox from [VLFeat Website](http://www.vlfeat.org/).
  - Extract it and make sure its path is correctly added to your MATLAB workspace.
  - Example: Ensure the `vl_setup` script from VLFeat is accessible.

---

## 2. **Setup**
1. **Add the Required Paths**:
   - Ensure all necessary files and functions (e.g., `feature_extraction`, `estimate_E_robust`, `Cheirality_triangulate`) are in your MATLAB path.
   - Run the `vl_setup` script to initialize VLFeat.

2. **Prepare the Dataset**:
   - Ensure the dataset files (e.g., images and camera calibration data) are in the required format and location.
   - Each dataset must be supported by the `get_dataset_info` function, which extracts intrinsic camera parameters and image names.

---

## 3. **Basic Execution**
Run the `run_sfm` function in MATLAB using the following syntax:

```matlab
run_sfm(dataset_num);
```
---
## Sample of the results:

![image](https://github.com/user-attachments/assets/03c75c94-62bf-4293-b61b-8694f0680bed)
