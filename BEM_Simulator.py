import numpy as np
import matplotlib.pyplot as plt
import BEM_Animator
import TDMS_Reader

## VARIABLE DEFINITION ________________________________________________________________________________________________________________

# Ambient
g0 = 9.80665            # m/s^2
Patm = 0                # Pa
Tatm = 0                # K
Tref = 273.15 + 10      # K

# Geometry
Dout = 0                # m
Dport = 0               # m
Lgrain = 0              # m
Lhalf = 0               # m
Lcase = 0               # m
din = 0                 # m
At = 0                  # m^2
Aratio = 0              # -
Pratio = 0              # -
Vadd = 0                # m^3
V01 = 0                 # m^3
V02 = 0                 # m^3

# Regression
r0 = 0                  # m/s
P0 = 0                  # Pa
n = 0                   # -
sigma = 0               # 1/K

# Chemistry
Ra = 8.314462618153     # J/Kmol
R = 0                   # J/kgK
Tc = 0                  # K
gamma = 0               # -
Gamma = 0               # -

# Propellant
rho = 0                 # kg/m^3

# Performance
cstar = 0               # m/s
eta_c = 0               # -
eta_n = 0               # -

# Time
t_end = 0               # s
dt = 0                  # s
t_ign_S = 0             # s
t_ign_T = 0             # s

# Angle
theta_0 = 0             # rad
theta_step = 0          # rad
theta_end = 0           # rad
tb_theta = 0            # s

# ALL DATA
data = ()

# TDMS READ DATA
# Config2: startPS - stopPS - dtPS
times = [5252*10**3, 5383*10**3, 4*10**(-5)]
# Config2: startLC - stopLC - dtLC
times += [10515*10**3, 10765*10**3, 2*10**(-5)]
# Config1: startPS - stopPS - dtPS
times += [1920*10**3 + 38*10**3, 2120*10**3, 4*10**(-5)]
# Config1: startLC - stopLC - dtLC
times += [2*1920*10**3 + 88*10**3, 2*2120*10**3, 2*10**(-5)]


## SETUP ______________________________________________________________________________________________________________________________

# Setup All Relevant Parameters
def setup(ambient, geometry, regression, chemistry, propellant, performance, time, angle, V00):
    """ambient = Patm, Tatm [°C ] \n
    geometry = Dout, Dport, Lgrain, Lhalf, Lcase, Dt, De \n
    regression = r0, P0, n, sigma \n
    chemistry = M, Tc, gamma \n
    propellant = rho_o, rho_f, OF \n
    performance = eta_c, eta_n \n
    time = t_end, dt, t_ign_S, t_ign_T \n
    angle = theta_0, theta_end, tb_theta \n"""

    global Patm, Tatm, Dout, Dport, Lgrain, Lhalf, Lcase, din, At, Aratio, Pratio, Vadd, r0, P0, n, sigma
    global R, Tc, gamma, Gamma, rho, cstar, eta_c, eta_n, t_end, dt, t_ign_S, t_ign_T, theta_0, theta_step, theta_end, tb_theta
    global data

    data = (ambient, geometry, regression, chemistry, propellant, performance, time, angle, V00)

    # Ambient
    Patm, Tatm = ambient

    # Geometry
    Dout, Dport, Lgrain, Lhalf, Lcase, din, Dt, De = geometry
    At = np.pi/4*Dt**2
    Ae = np.pi/4*De**2
    Aratio = Ae/At
    Vadd = V00

    # Regression
    r0, P0, n, sigma = regression

    # Chemistry
    R = Ra/chemistry[0]
    Tc, gamma = chemistry[1:]
    Gamma = np.sqrt(gamma)*(2/(gamma + 1))**((gamma + 1)/(2*gamma - 2))

    # Propellant
    rho_o, rho_f, OF = propellant
    rho = OF/(OF + 1)*rho_o + 1/(OF + 1)*rho_f

    # Performance
    cstar = np.sqrt(R*Tc/gamma*((gamma + 1)/2)**((gamma + 1)/(gamma - 1)))
    eta_c, eta_n = performance

    # Time
    t_end, dt, t_ign_S, t_ign_T = time

    # Angle
    theta_0, theta_end, tb_theta = angle
    theta_step = (theta_end - theta_0)/tb_theta*dt

    # Compute Initial Chamber Volume
    calculateV0()
    # Compute Pressure Ratio
    calculatePratio()

# Calculate Initial Chamber Volume
def calculateV0():
    global V01, V02

    # Port Volume
    V0 = np.pi/4*Dport**2*Lgrain
    # Additional Case Volume (Left and Right of Grain)
    V0 += (Lcase - Lgrain)*np.pi/4*Dout**2
    # CURRENT SETUP: Volume Lost by Adding Inhibitor Plates
    V01 = V0 - 2*din*np.pi/4*(Dout**2 - Dport**2)
    # NEW SETUP: Volume Gained by Adding Slit
    Lslit = Lgrain - 2*Lhalf
    V02 = V0 + 2*Lslit*np.pi/4*(Dout**2 - Dport**2)

    # Artificially Added Volume to Fit Model
    V01 += Vadd
    V02 += Vadd

# Calculate Pressure Ratio of the Nozzle
def calculatePratio():
    global Aratio, Pratio, Gamma, gamma

    eps = 10**(-8)
    # ! This Range is Based on the Known Characteristics of the Nozzle - for Different Nozzle: Restart With Large Range !
    pratios = np.arange(0.0492938, 0.0492939, 0.000000001)
    for prat in pratios:
        diff = abs(Aratio - Gamma/np.sqrt(2*gamma/(gamma - 1)*prat**(2/gamma)*(1 - prat**((gamma - 1)/gamma))))
        if diff < eps:
            Pratio = prat
            break


## SIMULATION _________________________________________________________________________________________________________________________

# Runge-Kutta Numerical Solver of the Fourth Order
def RK4(f, Y0, a, b, h):
    
    # Divide Time Domain [a, b] into Steps of Size h
    L = b - a
    N = int(L/h)
    t = np.linspace(a, b, N + 1)
    
    # Make RK4 Compatible with Both Normal (1D) and Vector (nD) Functions
    if type(Y0) == float or type(Y0) == int:
        Y0 = [Y0]
    M = len(Y0)
    
    # Y Will Store Solution
    Y = np.zeros((M, N + 1))
    Y[:, 0] = Y0
    
    for i in range(N):
        k1 = np.array(f(t[i], Y[:, i]))
        k2 = np.array(f(t[i] + h/2, Y[:, i] + h/2*k1))
        k3 = np.array(f(t[i] + h/2, Y[:, i] + h/2*k2))
        k4 = np.array(f(t[i] + h, Y[:, i] + h*k3))
        Y[:, i+1] = Y[:, i] + h/6*(k1 + 2*k2 + 2*k3 + k4)
        
    # Y is Vector Solution! --> Y[0] = All First Components as a Function of Time
    return Y, t

# Compute Input for RK: dY/dt = f(Y, t) - This Function Calculates that f
# ORIGINAL SETUP
def dY1(t, Y):

    # Prevent Chamber Pressure from Dropping Below Atmospheric
    Y[2] = max(Y[2], Patm)

    w, Vc, Pc, theta = Y

    # Define End Condition for Regression
    if w >= (Dout - Dport)/2:
        r = 0
    else:
        r = r0*np.exp(sigma*(Tatm - Tref))*(Pc/P0)**n

    # Calculate Burning Surface Area
    S = angleBurn1(w, theta)[0]
    
    T = 18/16*Tc
    cstar = np.sqrt(R*T/gamma*((gamma + 1)/2)**((gamma + 1)/(gamma - 1)))
    # Model Ignition: Only a Fraction of the Surface is Burning + Only at Fraction of Combustion Temperature
    if t < t_ign_S:
        frac = t/t_ign_S
        S *= frac
    if t < t_ign_T:
        frac = t/t_ign_T
        T *= frac

    # Compute Derivatives
    derivatives = [0, 0, 0, 0]
    derivatives[0] = r
    derivatives[1] = r*S
    derivatives[2] = R*T/Vc*(rho*S*r - Pc*At/(eta_c*cstar)) - Pc/Vc*derivatives[1]
    derivatives[3] = theta_step/dt

    return derivatives
# NEW SETUP
def dY2(t, Y):

    # Prevent Chamber Pressure from Dropping Below Atmospheric
    Y[2] = max(Y[2], Patm)

    w, Vc, Pc, theta = Y

    # Define End Condition for Regression
    if w >= Lhalf/2:
        r = 0
    else:
        r = r0**np.exp(sigma*(Tatm - Tref))*(Pc/P0)**n    

    # Calculate Burning Surface Area
    S = angleBurn2(w, theta)[0]

    T = Tc
    # Model Ignition: Only a Fraction of the Surface is Burning + Only at Fraction of Combustion Temperature
    if t < t_ign_S:
        frac = t/t_ign_S
        S *= frac
    if t < t_ign_T:
        frac = t/t_ign_T
        T *= frac

    # Compute Derivatives
    derivatives = [0, 0, 0, 0]
    derivatives[0] = r
    derivatives[1] = r*S
    derivatives[2] = R*T/Vc*(rho*S*r - Pc*At/(eta_c*cstar)) - Pc/Vc*derivatives[1]
    derivatives[3] = theta_step/dt

    return derivatives

# Calculate Burning Surface Area for a Given Web Thickness and Burn Angle
# ORIGINAL SETUP
def angleBurn1(w, theta):
    D = Dport + 2*w

    # Smallest Radius (at Igniter)
    R1 = np.minimum(D/2, Dout/2)
    # Largest Radius (at Nozzle, or at Wall Towards End of the Burn)
    R2 = np.minimum(R1 + Lgrain*np.tan(theta), Dout/2)

    # Length of the Grain that is Still Burning (Towards End of the Burn)
    # Handle 0-Division if theta is Just an int (Not an Array)
    if isinstance(theta, int):
        if theta != 0:
            Lx = (Dout/2 - R1)/np.tan(theta)
        else: 
            Lx = Lgrain
    else:
        Lx = (Dout/2 - R1)/np.tan(theta)

    # Length of the Truncated Cone Side
    s = np.minimum(Lx, Lgrain)/np.cos(theta)
    # Total Burning Surface Area
    S = np.pi*(R1 + R2)*s

    return S, R1, R2, np.minimum(Lx, Lgrain)
# NEW SETUP
def angleBurn2(w, theta):    
    D = Dport + 2*w
    L = Lhalf - 2*w

    # Radii and Lengths as Defined Below
    # --- = Empty Space
    # ___ = Unburnt Propellant

    # D/2 --- R11___R12 --- (Slit) --- R21___R22 ---
    #      w      L      w          w      L      w

    # Smallest Radius of First Half (at Igniter Side, But Some Part Burnt Away)
    R11 = np.minimum(D/2 + w*np.tan(theta), Dout/2)
    # Largest Radius of First Half (at Slit, or at Wall Towards End of the Burn)
    R12 = np.minimum(R11 + L*np.tan(theta), Dout/2)

    # Length of the First Half that is Still Burning (Towards End of the Burn)
    # Handle 0-Division if theta is Just an int (Not an Array)
    if isinstance(theta, int):
        if theta != 0:
            L1 = (Dout/2 - R11)/np.tan(theta)
        else: 
            L1 = L
    else:
        L1 = (Dout/2 - R11)/np.tan(theta)

    # Length of the First Truncated Cone Side
    s1 = np.minimum(L, L1)/np.cos(theta)

    # Smallest Radius of Second Half (at Slit, or at Wall Towards End of the Burn)
    R21 = np.minimum(D/2 + (Lgrain - w - L)*np.tan(theta), Dout/2)
    # Largest Radius of Second Half (at Nozzle, or at Wall Towards End of the Burn)
    R22 = np.minimum(D/2 + (Lgrain - w)*np.tan(theta), Dout/2)

    # Length of the Second Half that is Still Burning (Towards End of the Burn)
    # Handle 0-Division if theta is Just an int (Not an Array)
    if isinstance(theta, int):
        if theta != 0:
            L2 = (Dout/2 - R21)/np.tan(theta)
        else: 
            L2 = L
    else:
        L2 = (Dout/2 - R21)/np.tan(theta)

    # Length of the Second Truncated Cone Side
    s2 = np.minimum(L, L2)/np.cos(theta)

    # Largest End-Burning Surface (Igniter Side)
    S1 = np.pi*((Dout/2)**2 - R11**2)
    # Inner Surface of First Truncated Cone
    S2 = np.pi*(R11 + R21)*s1
    # Largest Mid-Burning Surface (Igniter Side)
    S3 = np.pi*((Dout/2)**2 - R12**2)

    # Smallest Mid-Burning Surface (Nozzle Side)
    S4 = np.pi*((Dout/2)**2 - R21**2)
    # Inner Surface of Second Truncated Cone
    S5 = np.pi*(R21 + R22)*s2
    # Smallest End-Burning Surface (Nozzle Side)
    S6 = np.pi*((Dout/2)**2 - R22**2)

    # Total Burning Surface Area
    S = S1 + S2 + S3 + S4 + S5 + S6

    return S, R11, R12, R21, R22, np.minimum(L, L1), np.minimum(L, L2)

# Calculate Thrust Based on Pressure And Cf
def thrust(Y):
    Pc = Y[2]
    Cf = Gamma*np.sqrt(2*gamma/(gamma - 1)*(1 - Pratio**((gamma - 1)/gamma))) + (Pratio - Patm/Pc)*Aratio
    F = eta_n*Cf*Pc*At
    F -= np.min(F)
    return F

# Simulate Both Burns
def simulate():

    # Numerically Solve CURRENT SETUP Using RK4
    Y01 = [0, V01, Patm, theta_0]
    Y1, t1 = RK4(dY1, Y01, 0, t_end, dt)
    w1, _, _, theta1 = Y1
    # Re-Compute Burning Surface Area at Every Point in Time
    surface1 = angleBurn1(w1, theta1)
    # Calculate Thrust at Every Point in Time
    F1 = thrust(Y1)

    # Numerically Solve NEW SETUP Using RK4
    Y02 = [0, V02, Patm, theta_0]
    Y2, t2 = RK4(dY2, Y02, 0, t_end, dt)
    w2, _, _, theta2 = Y2
    # Re-Compute Burning Surface Area at Every Point in Time
    surface2 = angleBurn2(w2, theta2)
    # Calculate Thrust at Every Point in Time
    F2 = thrust(Y2)

    return t1, Y1, surface1, F1, t2, Y2, surface2, F2

# Calculate Total Burn Time
def burnTime(P, t, Pmin):
    i = 0

    # Find Start Time (Moment When P Rises Above Pmin)
    while P[i] < Pmin:
        i += 1
    tb_start = t[i]

    # Find Stop Time (Moment When P Drops Below Pmin)
    while P[i] > Pmin or t[i] < tb_start + 1:
        i += 1
    tb_end = t[i]

    return tb_start, tb_end

# Calculate Total Impulse
def impulse(F, t, tb_start, tb_end):

    # Slice Thrust Profile to Include Only Burn Time
    i_start = np.where(t == tb_start)[0][0]
    i_stop = np.where(t == tb_end)[0][0]
    F = F[i_start:i_stop]

    # Calculate Total Area Using Trapezoidal Rule
    Itot = 0
    dt = t[1] - t[0]
    for i in range(len(F) - 1):
        Itot += (F[i] + F[i + 1])*dt/2
    
    return Itot


## POST-PROCESSING (USER INTERFACE) ___________________________________________________________________________________________________

# Simulate Chamber Pressure Profiles Throughout Burn
def pressurePlot(config1 = True, config2 = True, ref = False, calibration = False, image = False):

    # Simulate Burns
    t1, Y1, _, _, t2, Y2, _, _ = simulate()
    Pc1 = Y1[2]
    Pc2 = Y2[2]
    
    # Plot Both Pressure Profiles
    ax = referencePlot("P", config1, calibration, image = image) if ref else plt.axes()
    colour = 'r' if ref else ''
    label1 = "Simulation" if ref else "Original Setup"
    label2 = "Simulation" if ref else "New Setup"
    if config1:
        ax.plot(t1, Pc1/10**6, colour, label = label1)
    if config2:
        ax.plot(t2, Pc2/10**6, colour, label = label2)
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.grid(axis = 'y')
    if config1*config2 or ref:
        ax.legend()
    plt.show()

# Simulate Thrust Profiles Throughout Burn
def thrustPlot(config1 = True, config2 = True, ref = False, calibration = False, image = False):

    # Simulate Burns
    t1, _, _, F1, t2, _, _, F2 = simulate()

    # Plot Both Thrust Profiles
    # Plot Both Pressure Profiles
    ax = referencePlot("F", config1, calibration, image = image) if ref else plt.axes()
    colour = 'r' if ref else ''
    label1 = "Simulation" if ref else "Original Setup"
    label2 = "Simulation" if ref else "New Setup"
    if config1:
        ax.plot(t1, F1, colour, label = label1)
    if config2:
        ax.plot(t2, F2, colour, label = label2)
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$F$" + " " + r"$[N]$")
    ax.grid(axis = 'y')
    if config1*config2 or ref:
        ax.legend()
    plt.show()

def tdmsPlot(config1, P, calibration = False):
    ax = referencePlot("P", config1, calibration) if P else referencePlot("F", config1, calibration)
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$") if P else ax.set_ylabel(r"$F$" + " " + r"$[N]$")
    ax.grid(axis = 'y')
    plt.show()

# Plot Reference Data
def referencePlot(type, config1, calibration, image = False):
    
    name = "Data_Old.tdms" if config1 else "Data_New.tdms"
    _, ax = plt.subplots()

    # Pressure Plot
    if type == "P":

        # Background Image Data in Case of Old Configuration
        if config1 and image:
            img = plt.imread("Figures\Reference_Test_P.png")
            ax.imshow(img, aspect = "auto", extent = [0, t_end, 0, 8 + Patm/10**6])

        # TDMS Data
        #t, data = TDMS_Reader.readChannel(name, "PS0", 5220*10**3 + 32*10**3, 5420*10**3 - 37*10**3, 0.00004, 1000, calibration)
        t, data = TDMS_Reader.readChannel(name, "PS0", times[6*config1], times[6*config1 + 1], times[6*config1 + 2], 1000, calibration)
        ax.plot(t, data/10**6, label = "TDMS data")

    # Thrust Plot
    else:

        # Background Image Data in Case of Old Configuration
        if config1 and image:
            img = plt.imread("Figures\Reference_Test_F.png")
            ax.imshow(img, aspect = "auto", extent = [0, t_end, 0, 600])
        
        # TDMS Data
        #t, data = TDMS_Reader.readChannel(name, "LC1", 2*5220*10**3 + 75*10**3, 2*5420*10**3 - 75*10**3, 0.00002, 1000, calibration)
        t, data = TDMS_Reader.readChannel(name, "LC1", times[6*config1 + 3], times[6*config1 + 4], times[6*config1 + 5], 1000, calibration)
        ax.plot(t, data, label = "TDMS data")
    return ax

# Return Relevant Performance Characteristics
def characteristics(Pmin, config1 = True, config2 = True, calibration = False):

    # Simulate Burns
    t1, Y1, _, F1, t2, Y2, _, F2 = simulate()
    Pc1 = Y1[2]
    Pc2 = Y2[2]

    # Fetch Test Data
    if config1:
        t1a_test, Pc1_test = TDMS_Reader.readChannel("Data_Old.tdms", "PS0", times[6], times[7], times[8], 1000, calibration)
        t1b_test, F1_test = TDMS_Reader.readChannel("Data_Old.tdms", "LC1", times[9], times[10], times[11], 1000, calibration)
    if config2:
        t2a_test, Pc2_test = TDMS_Reader.readChannel("Data_New.tdms", "PS0", times[0], times[1], times[2], 1000, calibration)
        t2b_test, F2_test = TDMS_Reader.readChannel("Data_New.tdms", "LC1", times[3], times[4], times[5], 1000, calibration)


    if config1 or config2:
        print("________________________________________________________________________________ \n")

    if config1:
        print("ORIGINAL CONFIGURATION - Characteristics:\n")

        print("SIMULATED:\n")

        # Calculate Max Pressure and Thrust
        Pmax1 = np.max(Pc1)
        Fmax1 = np.max(F1)
        tPmax1 = t1[np.where(Pc1 == Pmax1)[0][0]]
        tFmax1 = t1[np.where(F1 == Fmax1)[0][0]]
        print("\t Maximum Chamber Pressure \t Pc = %0.2f MPa \t at t = %0.2f s" %(Pmax1/10**6, tPmax1))
        print("\t Maximum Thrust \t \t F = %0.2f N \t at t = %0.2f s" %(Fmax1, tFmax1))

        # Calculate Burn Time
        tb_start1, tb_end1 = burnTime(Pc1, t1, Pmin)
        print("\t Burn Time \t \t \t tb = %0.2f s \t from t = %0.2f to %0.2f s" %(tb_end1 - tb_start1, tb_start1, tb_end1))

        # Calculate Total Impulse
        Itot1 = impulse(F1, t1, tb_start1, tb_end1)
        print("\t Total Impulse \t \t \t Itot = %0.2f Ns" %(Itot1))

        # Calculate Specific Impulse
        Vc1_final = Y1[1][-1]
        Vp1 = Vc1_final - V01
        Mp1 = Vp1*rho
        Isp1 = Itot1/Mp1/g0
        print("\t Specific Impulse \t \t Isp = %0.2f s" %(Isp1))

        # Calculate Average Thrust
        Fav1 = Itot1/(tb_end1 - tb_start1)
        print("\t Average Thrust \t \t F = %0.2f N" %(Fav1))

        print("EXPERIMENT:\n")

        # Calculate Max Pressure and Thrust
        Pmax1_test = np.max(Pc1_test)
        Fmax1_test = np.max(F1_test)
        tPmax1_test = t1a_test[np.where(Pc1_test == Pmax1_test)[0][0]]
        tFmax1_test = t1b_test[np.where(F1_test == Fmax1_test)[0][0]]
        print("\t Maximum Chamber Pressure \t Pc = %0.2f MPa \t at t = %0.2f s" %(Pmax1_test/10**6, tPmax1_test))
        print("\t Maximum Thrust \t \t F = %0.2f N \t at t = %0.2f s" %(Fmax1_test, tFmax1_test))

        # Calculate Burn Time
        tb_start1_test, tb_end1_test = burnTime(Pc1_test, t1a_test, Pmin)
        print("\t Burn Time \t \t \t tb = %0.2f s \t from t = %0.2f to %0.2f s" %(tb_end1_test - tb_start1_test, tb_start1_test, tb_end1_test))

        # Calculate Total Impulse
        Itot1_test = impulse(F1_test, np.array(t1b_test), tb_start1_test, tb_end1_test)
        print("\t Total Impulse \t \t \t Itot = %0.2f Ns" %(Itot1_test))

        # Calculate Specific Impulse
        print("\t Specific Impulse \t \t Unknown (No Mass Data)")

        # Calculate Average Thrust
        Fav1_test = Itot1_test/(tb_end1_test - tb_start1_test)
        print("\t Average Thrust \t \t F = %0.2f N" %(Fav1_test))

        # Motor Quality
        eta = eta_c*eta_n
        print("\t Motor Quality \t \t \t eta = %0.2f %%" %(eta*100))


    if config1*config2:
        print("\n" + "-------------------------------------------------------------------------------- \n")
    
    if config2:

        print("NEW CONFIGURATION - Characteristics:\n")

        print("SIMULATED:\n")

        # Calculate Max Pressure and Thrust
        Pmax2 = np.max(Pc2)
        Fmax2 = np.max(F2)
        tPmax2 = t2[np.where(Pc2 == Pmax2)[0][0]]
        tFmax2 = t2[np.where(F2 == Fmax2)[0][0]]
        print("\t Maximum Chamber Pressure \t Pc = %0.2f MPa \t at t = %0.2f s" %(Pmax2/10**6, tPmax2))
        print("\t Maximum Thrust \t \t F = %0.2f N \t at t = %0.2f s" %(Fmax2, tFmax2))

        # Calculate Burn Time
        tb_start2, tb_end2 = burnTime(Pc2, t2, Pmin)
        print("\t Burn Time \t \t \t tb = %0.2f s \t from t = %0.2f to %0.2f s" %(tb_end2 - tb_start2, tb_start2, tb_end2))

        # Calculate Total Impulse
        Itot2 = impulse(F2, t2, tb_start2, tb_end2)
        print("\t Total Impulse \t \t \t Itot = %0.2f Ns" %(Itot2))

        # Calculate Specific Impulse
        Vc2_final = Y2[1][-1]
        Vp2 = Vc2_final - V02
        Mp2 = Vp2*rho
        Isp2 = Itot2/Mp2/g0
        print("\t Specific Impulse \t \t Isp = %0.2f s" %(Isp2))

        # Calculate Average Thrust
        Fav2 = Itot2/(tb_end2 - tb_start2)
        print("\t Average Thrust \t \t F = %0.2f N" %(Fav2))

        print("EXPERIMENT:\n")

        # Calculate Max Pressure and Thrust
        Pmax2_test = np.max(Pc2_test)
        Fmax2_test = np.max(F2_test)
        tPmax2_test = t2a_test[np.where(Pc2_test == Pmax2_test)[0][0]]
        tFmax2_test = t2b_test[np.where(F2_test == Fmax2_test)[0][0]]
        print("\t Maximum Chamber Pressure \t Pc = %0.2f MPa \t at t = %0.2f s" %(Pmax2_test/10**6, tPmax2_test))
        print("\t Maximum Thrust \t \t F = %0.2f N \t at t = %0.2f s" %(Fmax2_test, tFmax2_test))

        # Calculate Burn Time
        tb_start2_test, tb_end2_test = burnTime(Pc2_test, t2a_test, Pmin)
        print("\t Burn Time \t \t \t tb = %0.2f s \t from t = %0.2f to %0.2f s" %(tb_end2_test - tb_start2_test, tb_start2_test, tb_end2_test))

        # Calculate Total Impulse
        Itot2_test = impulse(F2_test, np.array(t2b_test), tb_start2_test, tb_end2_test)
        print("\t Total Impulse \t \t \t Itot = %0.2f Ns" %(Itot2_test))

        # Calculate Specific Impulse
        print("\t Specific Impulse \t \t Unknown (No Mass Data)")

        # Calculate Average Thrust
        Fav2_test = Itot2_test/(tb_end2_test - tb_start2_test)
        print("\t Average Thrust \t \t F = %0.2f N" %(Fav2_test))

        # Motor Quality
        eta = eta_c*eta_n
        print("\t Motor Quality \t \t \t eta = %0.2f %%" %(eta*100))



    if config1 or config2:
        print("________________________________________________________________________________ \n")

# Animate the Regression of the Grain Using the BEM_Animator.py File
def animate(fps, theta, config1 = True, config2 = True):

    # Calculate New dt Based on Desired FPS
    global dt
    dt = 1/fps

    # Use Different (Larger, More Visible) Burn Angles
    global theta_step, tb_theta
    theta_end = theta
    theta_step = (theta_end - theta_0)/tb_theta*dt

    # Simulate Burns
    _, _, surface1, _, _, _, surface2, _ = simulate()
    _, R1, R2, L = surface1
    _, R11, R12, R21, R22, L1, L2 = surface2

    # Animate Burns
    if config1:
        BEM_Animator.animate1(R1, R2, L)
        print("Animation 1 Succesfull")
    if config2:
        BEM_Animator.animate2(R11, R12, R21, R22, L1, L2)
        print("Animation 2 Succesfull")

# Performs a Sensitivity Analysis of the Pressure Profile on Different Parameters
def sensitivity(parameter, factors, config1 = True, config2 = True):
    """Parameter Options are: Vadd, r0, n, Tc, theta_end, eta_c, Patm, Tatm\n"""
    
    # Determine Which Parameter to Tweak
    if parameter == "Vadd":
        return sensV(factors, config1, config2) 
    elif parameter == "r0":
        return sensR(factors, config1, config2)
    elif parameter == "n":
        return sensN(factors, config1, config2)
    elif parameter == "Tc":
        return sensT(factors, config1, config2)
    elif parameter == "theta_end":
        return sensTheta(factors, config1, config2)
    elif parameter == "eta_c":
        return sensEta(factors, config1, config2)
    elif parameter == "Patm":
        return sensPatm(factors, config1, config2)
    elif parameter == "Tatm":
        return sensTatm(factors, config1, config2)
    else:
        return False

# Plot Reference Data (< Experiment)

## SENSITIVITY ANALYSIS________________________________________________________________________________________________________________

# Sensitivity to Added Chamber Volume
def sensV(factors, config1, config2):
    global Vadd, data

    # Determine Different Vadd's
    Vs = np.array(factors)*Vadd

    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for V in Vs:
        # Update Chamber Volume
        Vadd = V
        calculateV0()

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$V_0$ = " + str(round(V01*100**3, 2)) + r" $cm^3$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$V_0$ = " + str(round(V02*100**3, 2)) + r" $cm^3$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True

# Sensitivity to Reference Regression Rate
def sensR(factors, config1, config2):
    global r0

    # Determine Different r0's
    rs = np.array(factors)*r0
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for r in rs:
        # Update Reference Regression Rate
        r0 = r

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$r_0$ = " + str(round(r*1000, 2)) + r" $mm/s$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$r_0$ = " + str(round(r*1000, 2)) + r" $mm/s$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True

# Sensitivity to Combustion Index
def sensN(factors, config1, config2):
    global n

    # Determine Different n's
    ns = np.array(factors)*n
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for ni in ns:
        # Update Combustion Index
        n = ni

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$n$ = " + str(round(n, 2)))
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 +  r"$n$ = " + str(round(n, 2)))
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True

# Sensitivity to Chamber Temperature
def sensT(factors, config1, config2):
    global Tc, cstar, gamma, R

    # Determine Different Tc's
    Ts = np.array(factors)*Tc
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for T in Ts:
        # Update Chamber Temperature and Characteristic Velocity
        Tc = T
        cstar = np.sqrt(R*Tc/gamma*((gamma + 1)/2)**((gamma + 1)/(gamma - 1)))

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$T_c$ = " + str(Tc) + r" $K$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$T_c$ = " + str(Tc) + r" $K$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()  
    return True    

# Sensitivity to Burning Angle
def sensTheta(factors, config1, config2):
    global theta_step, theta_end

    # Determine Different thetas
    thetas = np.array(factors)*theta_end
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for theta in thetas:
        # Update Burning Angle
        theta_end = theta
        theta_step = (theta_end - theta_0)/tb_theta*dt

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$\theta_{end}$ = " + str(round(theta_end/np.pi*180, 2)) + r"$°$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$\theta_{end}$ = " + str(round(theta_end/np.pi*180, 2)) + r"$°$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True

# Sensitivity to Combustion Efficiency
def sensEta(factors, config1, config2):
    global eta_c

    # Determine Different etas
    etas = np.array(factors)*eta_c
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for eta in etas:
        # Update Combustion Effiency
        eta_c = eta

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$\eta_c$ = " + str(round(eta_c*100, 2)) + r" $\%$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$\eta_c$ = " + str(round(eta_c*100, 2)) + r" $\%$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True

# Sensitivity to Ambient Pressure
def sensPatm(factors, config1, config2):
    global Patm

    # Determine Different Patms
    Patms = np.array(factors)*Patm
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Thrust Profile
    plt.figure()
    ax = plt.axes()
    for P in Patms:
        # Update Atmospheric Pressure
        Patm = P

        # Simulate Burns
        t1, _, _, F1, t2, _, _, F2 = simulate()
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, F1, label = label1 + r"$P_a$ = " + str(round(Patm/101325, 2)) + r" $atm$")
        if config2:
            ax.plot(t2, F2, label = label2 + r"$P_a$ = " + str(round(Patm/101325, 2)) + r" $atm$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$F$" + " " + r"$[N]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True
    
# Sensitivity to Ambient Temperature
def sensTatm(factors, config1, config2):
    global Tatm

    # Determine Different Patms
    Tatms = np.array(factors)*Tatm
    
    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for T in Tatms:
        # Update Atmopsheric Temperature
        Tatm = T

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + r"$T_a$ = " + str(round(Tatm - 273.15, 2)) + r" $°C$")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + r"$T_a$ = " + str(round(Tatm - 273.15, 2)) + r" $°C$")
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True
    
# Sensitivity to Manufacturing Tolerances
def toleranceSensitivity(geometry_min, geometry_max, config1 = True, config2 = True):
    """Geometry = (Dout, Dport, Lgrain, Lhalf)\n"""
    global Dout, Dport, Lgrain, Lhalf

    # Nominal Geometry
    geometry_nom = data[1]
    geometries = [geometry_min, geometry_nom, geometry_max]
    labels = ["Minimal", "Nominal", "Maximal"]

    # Store Original Data!
    original = data

    # For Each, Simulate and Plot Pressure Profile
    plt.figure()
    ax = plt.axes()
    for i in range(len(geometries)):
        # Update Geometry
        Dout, Dport, Lgrain, Lhalf = geometries[i][:4]

        # Simulate Burns
        t1, Y1, _, _, t2, Y2, _, _ = simulate()
        Pc1 = Y1[2]
        Pc2 = Y2[2]
    
        # Plot Both Pressure Profiles
        label1 = label2 = ""
        if config1*config2:
            label1 = "Original, "
            label2 = "New, "
        if config1:
            ax.plot(t1, Pc1/10**6, label = label1 + labels[i] + " Propellant")
        if config2:
            ax.plot(t2, Pc2/10**6, label = label2 + labels[i] + " Propellant")
        
        # Print Out Performance Characteristics
        Pmin = 10**6
        print("\n \t" + labels[i] + " Propellant:")
        characteristics(Pmin, config1 = config1, config2 = config2)
    
    ax.set_xlabel(r"$t$" + " " + r"$[s]$")
    ax.set_ylabel(r"$P_c$" + " " + r"$[MPa]$")
    ax.legend()

    # Restore Original Data!
    setup(original[0], original[1], original[2], original[3], original[4], original[5], original[6], original[7], original[8])

    plt.show()
    return True
    
            






