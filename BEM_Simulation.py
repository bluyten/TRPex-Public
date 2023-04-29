import numpy as np
import BEM_Simulator

## SHOULD NOT BE TWEAKED -----------------------------------------------------------------------------------------

# Geometry
Lgrain_min = 104.9/1000  # m
Lgrain = 105/1000        # m
Lgrain_max = 105.1/1000  # m

Lhalf_min = 49.9/1000   # m
Lhalf = 50/1000         # m
Lhalf_max = 50.1/1000   # m

Lcase = 118.9/1000      # m

Dport_min = 24.95/1000  # m
Dport = 25/1000         # m
Dport_max = 25.05/1000  # m

Dout_min = 76.7/1000    # m
Dout = 76.8/1000        # m
Dout_max = 76.9/1000    # m

Dt = 8.37/1000          # m
De = 16.74/1000         # m
din = 1.5/1000          # m

# Chemistry
M = 39.86/1000          # kg/mol        NAKKA
gamma = 1.1362          # -             NAKKA

# Propellant
rho_f = 1.489*1000      # kg/m^3
rho_o = 2.109*1000      # kg/m^3
OF = 65/35              # -

## CAN BE SLIGHTLY TWEAKED ---------------------------------------------------------------------------------------

# Time
t_end_1 = 6             # s
t_end_2 = 5             # s
dt = 1/1000             # s

# Chemistry
Tc = 1600               # K             (Nakka: 1600)

# Regression
P0 = 1*10**6            # Pa            (Nakka: 8.4 mm/s @ 10 MPa)
r0_1 = 3.90/1000        # m/s           
r0_2 = (8.4/1000)*0.435 # m/s
n_1 = 0.30              # -             (Nakka: 0.225)
n_2 = 0.39              # -
sigma = 0.0033          # 1/K           (Olde)

# Geometry
Vadd_1 = 0.0007         # m^3
Vadd_2 = 0.0008         # m^3

## SHOULD BE TWEAKED ---------------------------------------------------------------------------------------------

# Ambient
Patm_1 = 101325         # Pa
Patm_2 = 102660         # Pa
Tatm_1= 273.15 + 5      # K
Tatm_2 = 273.15 + 3.9   # K

# Performance
eta_c_1 = 1             # -
eta_c_2 = 0.97          # -
eta_n_1 = 0.85          # -
eta_n_2 = 0.86          # -

# Angle
theta_0 = 0.001/180*np.pi       # rad
theta_end_1 = 3.5/180*np.pi     # rad
theta_end_2 = 0.001/180*np.pi   # rad
tb = 5                          # s

# Time
t_ign_S_1 = 0.3         # s
t_ign_S_2 = 0.5         # s
t_ign_T_1 = 0.8         # s
t_ign_T_2 = 1.2         # s

# ________________________________________________________________________________________________________________
# ________________________________________________________________________________________________________________

ambient_1 = (Patm_1, Tatm_1)
ambient_2 = (Patm_2, Tatm_2)
geometry = (Dout, Dport, Lgrain, Lhalf, Lcase, din, Dt, De)
regression_1 = (r0_1, P0, n_1, sigma)
regression_2 = (r0_2, P0, n_2, sigma)
chemistry = (M, Tc, gamma)
propellant = (rho_o, rho_f, OF)
performance_1 = (eta_c_1, eta_n_1)
performance_2 = (eta_c_2, eta_n_2)
time_1 = (t_end_1, dt, t_ign_S_1, t_ign_T_1)
time_2 = (t_end_2, dt, t_ign_S_2, t_ign_T_2)
angle_1 = (theta_0, theta_end_1, tb)
angle_2 = (theta_0, theta_end_2, tb)

# Setup Simulator With Data
# Configuration 1 Fit
BEM_Simulator.setup(ambient_1, geometry, regression_1, chemistry, propellant, performance_1, time_1, angle_1, Vadd_1)
# Configuration 2 Fit --> COMMENT IF LOOKING AT CONFIG 1
BEM_Simulator.setup(ambient_2, geometry, regression_2, chemistry, propellant, performance_2, time_2, angle_2, Vadd_2)

# ________________________________________________________________________________________________________________
# ________________________________________________________________________________________________________________


# Simulate Pressure and Thrust Profiles --> UNCOMMENT TO USE


#BEM_Simulator.pressurePlot(config1 = False, config2 = True, ref = True, calibration = False, image = False)
BEM_Simulator.thrustPlot(config1 = False, config2 = True, ref = True, calibration = False, image = False)
#BEM_Simulator.tdmsPlot(config1 = True, P = True)


# ----------------------------------------------------------------------------------------------------------------


# Calculate Relevant Performance Characteristics --> UNCOMMENT TO USE


# Threshold for Defining Burn Start and End
Pmin = 1*10**6

#BEM_Simulator.characteristics(Pmin, config1 = False, config2 = True)


# ----------------------------------------------------------------------------------------------------------------


# Animate Regression of Grain --> UNCOMMENT TO USE


# IMPORTANT: This Code Only Works if the User Has GhostScript Installed and Added to PATH
# ALSO: Make sure that the following file structure is created
#                              *This Folder*/Frames/eps/
#                                                  /png/
#                                           BEM_Simulation.py
#                                           BEM_Simulator.py
#                                           BEM_Animator.py
fps = 20
theta_end = 10/180*np.pi

#BEM_Simulator.animate(fps, theta_end)


# ----------------------------------------------------------------------------------------------------------------


# Perform Sensitivity Analysis --> UNCOMMENT TO USE


P_MIN = 920     # mbar --> All Time Low Europe
P_NOM = 1013.25 # mbar
P_MAX = 1070    # mbar --> All Time High Europe
P_factors = [P_MIN/P_NOM, 1, P_MAX/P_NOM]

T_MIN = 273.15 - 5  # K
T_NOM = 273.15 + 10 # K
T_MAX = 273.15 + 25 # K
T_factors = [T_MIN/T_NOM, 1, T_MAX/T_NOM]

factors = [0.8, 0.9, 1]
parameter = "Tatm"                            # Can be Vadd, r0, n, Tc, theta_end, eta_c, Patm or Tatm
success = True

#success = BEM_Simulator.sensitivity(parameter, T_factors, config1 = False)

if not success:
    print("Error: Parameter for Sensitivity Analysis Was Unvalid!")


# ----------------------------------------------------------------------------------------------------------------


# Perform Tolerance Sensitivity Analysis --> UNCOMMENT TO USE


geometry_min = (Dout_min, Dport_max, Lgrain_min, Lhalf_min)
geometry_max = (Dout_max, Dport_min, Lgrain_max, Lhalf_max)

#BEM_Simulator.toleranceSensitivity(geometry_min, geometry_max, config1 = False)


# ________________________________________________________________________________________________________________
# ________________________________________________________________________________________________________________