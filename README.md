# Robot Kinematics & Trajectory Analysis Toolkit

## Project Overview

This MATLAB toolkit provides a comprehensive solution for analyzing and simulating a multi-degree-of-freedom robotic manipulator. The project encompasses forward kinematics, inverse kinematics, Jacobian analysis, workspace characterization, and dynamic trajectory planning. It enables detailed kinematic studies and 3D visualization of robot motion along predefined trajectories.

## Core Components

### 1. Forward Kinematics (MGD)

**Files:** `MGD_angle.m`, `MGD_ortho.m`

#### Overview
The forward kinematics functions compute the end-effector Cartesian position and orientation given joint angles. These functions form the foundation for all kinematic analyses.

#### Function: MGD (Geometric Direct Model)
- **Input:** Joint position vector `q` (joint angles in radians)
- **Output:** Homogeneous transformation matrix representing the end-effector pose relative to the base frame
- **Purpose:** Establishes the mapping from joint space to Cartesian space using Denavit-Hartenberg parameters

**Usage:**
```matlab
q = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6];  % Joint configuration
T = MGD(q);  % Returns 4x4 transformation matrix
```

#### Function: MGD_angle
- **Input:** Joint position vector `q`
- **Output:** Cartesian state including position (x, y, z) and orientation (roll, pitch, yaw)
- **Purpose:** Provides a user-friendly representation by extracting position and Euler angles from the transformation matrix
- **Utility:** Commonly used in trajectory planning and visualization scripts

**Usage:**
```matlab
q = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6];
[x, y, z, roll, pitch, yaw] = MGD_angle(q);  % Cartesian state
```

#### Variant: MGD_ortho
- **Specialization:** Optimized version for manipulators with orthogonal shoulder axes
- **Performance:** Faster computation for this specific kinematic configuration
- **Use Case:** Deploy when working with orthogonal arm configurations for improved computational efficiency

---

### 2. Jacobian Matrix Analysis

**Files:** `my_exo_jac.m`, `my_exo_orth_jac.m`

#### Overview
The Jacobian matrix maps velocities from joint space to Cartesian space, enabling analysis of manipulator dexterity and singularity detection.

#### Function: my_exo_jac (General Jacobian)
- **Input:** Joint position vector `q`
- **Output:** 6×6 Jacobian matrix (3 translational + 3 rotational DOF)
- **Purpose:** Computes the geometric Jacobian for general arm configurations
- **Mathematical Form:** $J = \frac{\partial x}{\partial q}$ where x is the end-effector pose

**Usage:**
```matlab
q = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6];
J = my_exo_jac(q);  % Compute Jacobian

% Velocity mapping: x_dot = J * q_dot
q_velocities = [q_dot_1; q_dot_2; q_dot_3; q_dot_4; q_dot_5; q_dot_6];
cartesian_velocity = J * q_velocities;
```

#### Function: my_exo_orth_jac (Orthogonal Jacobian)
- **Specialization:** Optimized for orthogonal shoulder configurations
- **Accuracy:** More precise singularity analysis for orthogonal kinematics
- **Integration:** Used in conjunction with `MGD_ortho` for consistent results

#### Applications
- **Singularity Detection:** Locate configurations where det(J) = 0
- **Manipulability Analysis:** Assess dexterity at different poses
- **Velocity Scaling:** Inverse kinematics with velocity constraints
- **Redundancy Resolution:** Handle kinematic redundancy in path planning

---

### 3. Inverse Kinematics (MGI)

**Files:** `MGI_NewtonRaphson_contraintes.m`, `MGI_NewtonRaphson_contraintes_poids.m`

#### Overview
These functions solve the inverse kinematics problem using iterative Newton-Raphson methods, enabling the robot to reach desired Cartesian positions.

#### Function: MGI_NewtonRaphson_contraintes (Newton-Raphson IK)
- **Inputs:**
  - `xd`: Desired Cartesian position [x, y, z] or full pose
  - `q0`: Initial joint configuration (critical for convergence)
- **Output:** Joint configuration `q` achieving the desired end-effector position
- **Algorithm:** Iterative numerical solver with Jacobian-based updates

**Usage:**
```matlab
desired_position = [x_target, y_target, z_target];
initial_config = [0, 0, 0, 0, 0, 0];  % Starting configuration
q_solution = MGI_NewtonRaphson_contraintes(desired_position, initial_config);
```

**Algorithm Principle:**
$$q_{k+1} = q_k + J^{-1}(q_k) \cdot (x_d - x(q_k))$$

- Iterates until convergence (position error below threshold)
- Requires non-singular Jacobian for stability
- Convergence depends on initial guess proximity to solution

#### Function: MGI_NewtonRaphson_contraintes_poids (Weighted IK)
- **Enhancement:** Weighted Jacobian pseudo-inverse for redundancy resolution
- **Inputs:** Same as standard MGI, with weight configuration internal to function
- **Functionality:** Prioritizes optimization objectives (e.g., joint limit avoidance, singularity avoidance)

**Weight Configuration:**
```matlab
% Modify weights directly in function body to prioritize:
% - Joint comfort (minimize excessive movement)
% - Singularity avoidance (maintain manipulability)
% - Obstacle avoidance (constraint weighting)
% - Energy efficiency (weighted joint velocities)
```

**Documentation of Internal Weights:**
The function contains weight matrices that can be tuned for application-specific constraints:
```matlab
W = diag([w1, w2, w3, w4, w5, w6]);  % Joint space weights
```

**Usage Example with Weight Optimization:**
```matlab
desired_position = [x_target, y_target, z_target];
initial_config = [0, 0, 0, 0, 0, 0];
q_weighted = MGI_NewtonRaphson_contraintes_poids(desired_position, initial_config);
% Prioritizes objectives defined by internal weights
```

#### Convergence Considerations
- **Convergence Criterion:** Euclidean distance between current and desired position < tolerance
- **Maximum Iterations:** Built-in safeguard to prevent infinite loops
- **Failure Cases:** Returns last best estimate if divergence detected
- **Initial Guess Importance:** Choose q0 close to expected solution for reliable convergence

---

### 4. Workspace & Manipulability Analysis

**File:** `espace_de_travail.m`

#### Overview
Comprehensive analysis of the robot's reachable workspace and dexterity characteristics across the entire workspace.

#### Functionality
The workspace analysis script generates four complementary visualizations:

1. **3D Workspace Envelope**
   - Visualization of all reachable points by the end-effector
   - Generated through dense sampling and Monte Carlo exploration
   - Shows accessibility across full workspace

2. **Smoothed Workspace Boundary**
   - Convex hull or alpha-shape representation of workspace
   - Reveals workspace topology and accessible regions
   - Identifies workspace boundaries and singularity surfaces

3. **Planar Cross-Section (z = constant)**
   - 2D slice of workspace at fixed height coordinate
   - Useful for path planning in horizontal planes
   - Shows reachable workspace area at specific heights

4. **Manipulability Heatmap**
   - Gradient color map representing dexterity at each point
   - Calculated from Jacobian determinant: $\mu = |\det(J)|$
   - Red regions: High manipulability (good dexterity)
   - Blue regions: Low manipulability (singularity proximity)

#### Usage
```matlab
espace_de_travail;  % Execute to generate all four figures
% Figures display:
% - Figure 1: Full 3D workspace
% - Figure 2: Workspace boundary
% - Figure 3: Z-slice cross-section
% - Figure 4: Manipulability gradient visualization
```

#### Applications
- **Path Planning:** Identify high-dexterity corridor regions
- **Task Reachability:** Verify if target positions lie within workspace
- **Singularity Avoidance:** Navigate around low-manipulability zones
- **Performance Characterization:** Compare different configurations or designs

---

### 5. Trajectory Generation

**File:** `trajectoire_en_cloche.m`

#### Overview
Generates smooth, bell-shaped (Gaussian-like) end-effector trajectories with kinematic compatibility.

#### Functionality
- **Trajectory Profile:** Smooth polynomial trajectory with continuous velocity and acceleration
- **Bell-Curve Shape:** Ensures gradual start/stop and smooth motion profile
- **Parameterization:** Time-based trajectory with configurable duration and control points
- **Inverse Kinematics Integration:** Converts Cartesian trajectory to joint space commands

#### Usage
```matlab
% Execute trajectory generation
trajectoire_en_cloche;  % Generates trajectory through joint space
% Output: Joint trajectory array (N×6) containing joint values over time
```

#### Trajectory Characteristics
- **Smoothness:** Continuously differentiable for precise control
- **Velocity Bounds:** Enforces maximum joint velocities
- **Acceleration Constraints:** Respects actuator acceleration limits
- **Collision Avoidance:** Can be extended with obstacle checking

#### Applications
- **Smooth Task Execution:** Reduces actuator strain and vibration
- **Real-time Control:** Pre-computed trajectory following
- **Performance Benchmarking:** Baseline trajectory for controller testing

---

### 6. 3D Robot Animation & Visualization

**File:** `check.m`

#### Overview
Real-time 3D visualization of the robot performing predefined trajectories with complete kinematic rendering.

#### Functionality
- **Real-time Animation:** 3D rendering of robot links and joints moving along trajectory
- **Trajectory Visualization:** Displays end-effector path as motion progresses
- **Frame Display:** Shows coordinate frames at each joint
- **Dynamic Updates:** Updates link positions based on forward kinematics

#### Usage
```matlab
check;  % Execute to launch 3D animation
% Display shows:
% - Robot structure in current configuration
% - Trajectory path traced by end-effector
% - Animation of robot motion along trajectory
% - Interactive 3D viewport
```

#### Configuration & Integration
The animation requires specification of:
1. **Inverse Kinematics Function:** Selects which MGI solver to use:
   - `MGI_NewtonRaphson_contraintes`: Standard IK
   - `MGI_NewtonRaphson_contraintes_poids`: Weighted optimization
   
2. **Trajectory Source:** Uses trajectory from `trajectoire_en_cloche.m` or custom trajectory

#### Important Notes
- **IK Function Consistency:** Ensure the script uses consistent IK function (with or without weights)
- **Synchronization:** For weighted IK trajectories, call `MGI_NewtonRaphson_contraintes_poids`
- **Performance:** Rendering speed depends on simulation step size and trajectory resolution

**Safety Consideration:** When switching between weighted and unweighted IK, verify the correct function is invoked in `check.m`:
```matlab
% For standard IK:
q = MGI_NewtonRaphson_contraintes(xd, q0);

% For weighted IK:
q = MGI_NewtonRaphson_contraintes_poids(xd, q0);
```

---

## Integration & Workflow

### Typical Analysis Workflow
1. **Define Goal:** Specify desired end-effector path
2. **Workspace Verification:** Run `espace_de_travail.m` to ensure reachability
3. **Trajectory Planning:** Use `trajectoire_en_cloche.m` for smooth Cartesian path
4. **Inverse Kinematics:** Compute joint trajectories via `MGI_NewtonRaphson_contraintes`
5. **Visualization:** Execute `check.m` to observe 3D robot motion
6. **Performance Analysis:** Use `my_exo_jac.m` to detect singularities along trajectory

### Function Dependencies
```
check.m (animation)
├── trajectoire_en_cloche.m (trajectory)
├── MGI_NewtonRaphson_contraintes.m (IK)
├── MGD_angle.m (forward kinematics)
└── my_exo_jac.m (Jacobian analysis)

espace_de_travail.m (workspace analysis)
├── MGD_angle.m
└── my_exo_jac.m
```

---

## System Requirements

- **MATLAB Version:** R2020b or later
- **Toolboxes:** Symbolic Math Toolbox, Robotics System Toolbox (optional for visualization)
- **Performance:** Modern multi-core processor recommended for dense workspace sampling

---

## Notes & Best Practices

- **Initial Guess Quality:** For MGI convergence, ensure initial joint configuration q0 is reasonably close to solution
- **Singularity Awareness:** Always check Jacobian determinant before relying on IK results near singularities
- **Weight Tuning:** When using weighted IK, document weight rationale for reproducibility
- **Workspace Planning:** Use manipulability heatmap to identify optimal working regions
- **Continuous Integration:** Combine multiple tools (MGD + Jacobian + MGI) for robust kinematic solutions


