# Gait Kinematics and Joint Coordinate Systems

Portfolio documentation for an academic movement-analysis project that reconstructs lower-limb kinematics from marker trajectories using joint coordinate systems and SVD-based rigid-body registration.

For a fuller description of the project, see the [project overview](https://ilalunghini.github.io/projects/gait-kinematics-jcs-analysis.html).

## Availability

The main analysis script is included as a workflow reference. Motion-capture data, derived outputs, and helper-function implementations are not yet included in this repository.

## Included code

- `src/gait_kinematics_jcs_analysis.m` — the main workflow, organized into data loading, static calibration, dynamic registration, joint-angle analysis, gait-cycle normalization, and angular-velocity analysis.

## Workflow functions

The original workflow uses the following helper functions. They are described here for methodological transparency but are not included in this repository.

| Function | Purpose |
| --- | --- |
| `computeRotationMatrix` | Estimates the best-fitting segment rotation at each frame using SVD. |
| `coordChange` | Expresses 3D marker coordinates in a different reference frame. |
| `CrossProduct` | Builds frame-wise cross-covariance matrices for rigid-body registration. |
| `Tinv` | Computes the inverse of a rigid homogeneous transformation. |
| `eventsnormalize` | Resamples gait-cycle signals to a common 0–100% time scale. |
| `plotCS` | Draws 3D coordinate systems for visual inspection. |
