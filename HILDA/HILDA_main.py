"""
Boxes:
    I:  stratified interior of depth D and width 1-delta
    LS: surface layer of I of depth hs
    HD: deep, well-mixed, polar layer of depth D and width delta
    HS: surface layer of HD of depth hs
    Atmosphere

    D:          Average depth of ocean                                  = 3800 m
    hs:         depth of surface layers                                 = 50 m
    delta:      width of total polar zone                               = 0.112-0.221 -> 0.16
    delta_s:    part of HS free of sea ice                              = 0.10

Absolute Areas and Area Fractions by Ocean Basin in the Low Latitude Ocean:
            North Pacific   South Pacific   North Atlantic  South Atlantic  North Indian    South Indian
    z,      Area    Frac.   Area    Frac.   Area    Frac.   Area    Frac.   Area    Frac.   Area    Frac.
    m       10**13 m**2     10**13 m**2     10**13 m**2     10**13 m**2     10**13 m**2     10**13 m**2
    0       7.79    0.259   8.05    0.268   4.43    0.148   3.36    0.112   1.17    0.039   5.21    0.174
    1000    7.29    0.263   7.80    0.282   3.69    0.133   3.15    0.114   0.99    0.036   4.75    0.172
    2000    7.03    0.268   7.40    0.282   3.35    0.128   3.07    0.117   0.87    0.033   4.49    0.171
    3000    6.60    0.280   6.71    0.284   2.83    0.120   2.84    0.120   0.67    0.029   3.97    0.168

ocean physics (deep, thermohaline circulation) are represented by four parameters:
    k:  an eddy diffusion coefficient                                   = 3.2*10**-5 m^2/s
    w:  a deep upwelling velocity in the stratified interior (I);       = 2.0*10**-8 m/s
    q:  a rate of lateral exchange between I and HD;                    = 7.5*10**-11 1/s
    u:  an exchange velocity between HS and HD                          = 1.9*10**-6 m/s

two-way exchange between the atmosphere and the surface layers of the ocean:
    g_LS: gas exchange velocity for low-latitude ocean surface layer    = 15.1 mol/m**2/yr  , at 280 ppm
    g_HS: gas exchange velocity for high-latitude ocean surface layer   = 15.1 mol/m**2/yr  , at 280 ppm

(Bidirectional exchange between the thin low- and high-latitude surface layers is neglected as is gas exchange through
sea ice)

Temperature:
    {T}(x,y,0): annual mean sea surface temperature                     = 2-8 °C -> 5°C
    T_LS:       low latitude surface layer temperature                  = 21.00 °C (winter: 19.54°C)
    T_HS:       high latitude surface layer temperature                 = 0.53 °C (winter: -0.34°C)

Carbon 14:
    lamda:      radioactive decay rate of carbon 14                     = 3.84*10**-12 1/s
    {14C}:      preindustrial, ocean mean value of carbon 14            = 0.84
    delta14C                                                            = -160 %%
    14C_LS:     preindustrial delta14C values in range -30 to -50%%     = 0.96
    delta14C_LS                                                         = -40 %%
    14C_HS:     preindustrial carbon 14 in HSv                          = 0.88
    delta14C_HS                                                         = -120 %%
    g_LS:       derived from (15) with 14C_LS and 14C_HS given          = 2.32*10**-7 m/s
    sumCO2                                                              = 2.33 mol/m**3
    fluxCO2                                                             = 5.41*10**-7 mol C/m**2/s; 17.1 mol C/m**2/yr

steady, conservation equations for tracer P:
    z:      vertical coordinate defined positive upward
    P:      tracer
    P_I(z): concentration in I at vertical coordinate z
    P_LS:   concentration in LS
    P_HD:   concentration in HD
    P_HS:   concentration in HS
    P_A:    concentration in the Atmosphere
    S_I(z): distribution of sources and sinks in I
    S_LS:   distribution of sources and sinks in LS
    S_HD:   distribution of sources and sinks in HD, includes sources and sinks on the bottom of HD
    S_HS:   distribution of sources and sinks in HS

    I:
    (1) k*d**2*P_I(z)/d*z**2 - w*d*P_I(z)/d*z - q*(P_I(z)-P_HD) + S_I(z) = 0

    HD:
    (2) (1-delta)*w*(P_HS-P_HD) + delta*u*(P_HS-P_HD) + (1-delta)*q*Integral[-D,0](P_I(z)-P_HD)*d*z + delta * S_HD = 0

        (3) q*Integral[-D,0](P_I(z)-P_HD)*d*z = q*D*([P_I](z)-P_HD)
                [P_I](z): average concentration in I
                q*D: exchange velocity

    LS:
    (4) g_LS(P_A-P_LS) - k*d*P_I(z)/d*z|_z=0 + [w*(P_I(z)|_z=0-P_LS)] + S_LS = 0
            equation may be considered as an upper boundary condition on the interior

        matching condition between I and LS:
        (5) P_I(z) - P_LS = 0, z=0
            -> the advective term in square brackets in (4) drops out

    HS:
    (6) delta_S*g_HS*(P_A-P_HS) + (1-delta)*w*(P_LS-P_HS) - delta*u*(P_HS-P_HD) + delta*S_HS = 0

    the net flux out of the bottom boundary balances the net influx plus net sources on the bottom:
    (7) k*d*P_I(z)/d*z-w*(P_I(Z)-P_HD)+ S_B = 0, z=-D
            S_B = distribution of sources and sinks on the ocean bottom

Temperature:
    s1,s2 = 0.5*[w*k**-1+-(w**2*k**-2+4*q*k**-1)**0.5]
    a4 =    (w-k*s2)*exp(-s2*D)*[(w-k*s1)*exp(-s1*D)]**-1-1
    a3 =    1 - s1*s2**-1 + a4*[1 + delta*(k*s1*(1-delta))**-1*u]
    a2 =    a3**-1*[a4*(k*s1)**-1*(w + delta*(1-delta)**-1*u)]
    a1 =    a3**-1*[(1-s2*s1**-1) + a4*(1-w*(k*s1)**-1)]


    general solution of (1) for the interior temperature with S_I(z)=0:
    (8) T_I(z) = T_HD + A1*exp(s1*z) + B1*exp(s2*z)

    Application of conditions (5) and (7) given (8) and S_B=0 yields:
    (A1)    T_LS-T_HD-A1-B1 = 0
    (A2)    A1*(k*s1-w)*exp(-s1*D) + B1*(k*s2-w)*exp(-s2*D) = 0

    condition of no net heat flux to the deep ocean:
    (9) w*(T_HS-T_LS) + delta*(1-delta)**-1*u*(T_HS-T_HD) + k*d*T_I(z)/d*z|_z = 0

    Obtained from (8) and (9):
    (A3)    w*(T_HS-T_LS) + delta*(1-delta)**-1*u*(T_HS-T_HD) + k*s1*A1 + k*s2*B1 = 0

    solution of (A1)-(A3) for the temperature in HD:
    (10)    T_HD = a1*T_LS + a2*T_HS

14C:
    14C stands henceforth for the 13C fractionation-corrected ratio of 14C/12C [cf. Siegenthaler, 1986]
    the error introduced in neglecting biological effects with this approach is less than 10% of the signal produced
    by circulation and radioactive decay [Fiadeiro, 1982; Bacastow and Maier-Reimer, 1990]

    s3,s4 = 0.5*j[w*k**-1+-(w**2*k**-2+4*(q+lamda)*k**-1)**0.5]

    conversion to standard d14C units:
    (11)    d14C = 1000*(model units -1)

    d14C in preindustrial atmosphere = 0 %%(permil)

    S_I(z) = -lamda*14C_I(z)
    S_HD = -lamda*D*14C_HD

    solution of (1) for 14C:
    (12)    14C_I(z) = q*(q+lamda)**-1*14C_HD + A2*exp(s3*z) + B2*exp(s4*z)

    (5) and (7) together with (12) yields, given S_B=0 for an abiotic approach:
    (A4)    14C_LS-q*(q+lamda)**-1*(14C_HD)-A2-B2 = 0
    (A5)    A2*(k*s3-w)*exp(-s3*D) + B2*(k*s4-w)*exp(-s4*D) + lamda*w*(q+lamda)**-1*(14C_HD) = 0

    neglecting radioactive decay in the thin surface layers (->S_LS, S_HS = 0), (4) and (6) yield:
    (A6)    g_LS*(1-14C_LS) - k*s3*A2 - k*s4*B2 = 0
    (A7)    delta_S*g_HS*(1-14C_HS) + (1-delta)*w*(14C_LS-14C_HS) - delta*u*(14C_HS-14C_HD) = 0

    combination of (2) and the net air-sea exchange - radioactive decay balance yields:
    (13)    (1-delta)*g_LS*(1-14C_LS) + delta_S*g_HS*(1-14C_HS) - (1-delta)*lamda*Integral[-D,0]14C_I(z)*d*z -
            delta*lamda*D*14C_HD = 0

    canceling terms involving the vertical integral of 14C_I(z) of (2) and (13) yields:
    (A8)    q*[(1-delta)*g_LS] - q*(1-delta)*g_LS*14C_LS + 14C_HS*[(1-delta)*lamda*w+delta*lamda*u-q*delta_S*g_HS] -
            lamda*14C_HS*[(1-delta)*w+delta*u+D*(q+delta*lamda)] = 0

    conservation equation for 14C:
    (14)    (1-delta)*g_LS*(1-14C_LS) + delta_s*g_HS*(1-14C_HS) - lamda*D*{14C} = 0

    better expression for exchange of 14C across the high-latitude sea surface:
        {delta_s(t)*g_HS(t)*[1-14C_HS(t)]} , the angle brackets represent yearly averaging

    g_LS for carbon 14
    (15)    g_LS = 1.23*10**-8*[(1-delta)*(1-14C_LS) + delta*(1-14C_HS)]**-t
    Oeschger & Siegenthaler modell -> implementieren von 14C Produktionsraten
    Eduart Bart
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Function representing the system of equations
def model(y, t):
    pass
    # Implement your system of equations here
    # For example: dydt = ... (Equations for temperature and concentration)

    # return dydt

def main():
    # Constants and parameters:
    # Box parameters:
    D = 3800  # Average depth of ocean (m)
    hs = 50  # Depth of surface layers (m)
    delta = 0.16  # Width of total polar zone
    delta_s = 0.10 # Part of HS free of sea ice
    # Ocean physics parameters:
    k = 3.2 * 10 ** -5 # Eddy diffusion coefficient (m**2/s)
    w = 2.0 * 10 ** -8 # Deep upwelling velocity in I (m/s)
    q = 7.5 * 10 ** -11 # Rate of lateral exchange between I and HD (1/s)
    u = 1.9 * 10 ** -6 # Exchange velocity between HS and HD (m/s)
    # Air-sea gas exchange parameters:
    g_LS = 15.1 # Gas exchange velocity for LS (mol/m**2/yr), at 280 ppm
    g_HS = 15.1 # Gas exchange velocity for HS (mol/m**2/yr), at 280 ppm
    # Temperature:
    T_ini_x_y_0 = 5 # Annual mean sea surface temperature: 2-8°C -> 5°C
    T_LS = 21.00 # low latitude surface layer temperature (°C) (winter: 19.54°C)
    T_HS = 0.53 # high latitude surface layer temperature (°C) (winter: -0.34°C)
    # Carbon 14:
    lamda = 3.84*10**-12 # radioactive decay rate of carbon 14 (1/s)



    # Initial conditions
    # ... (set up initial conditions for concentrations and temperatures)

    # Time settings
    num_time_steps = 100
    time = np.linspace(0, 10, num_time_steps)  # Adjust time range as needed

    # Solve the system of equations over time
    solution = odeint(model, initial_conditions, time)

    # Plot concentrations and temperatures
    plt.plot(time, solution[:, 0], label='P_I')
    plt.plot(time, solution[:, 1], label='T_I')
    # ... (plot other variables)

    plt.xlabel('Time')
    plt.ylabel('Concentration/Temperature')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()