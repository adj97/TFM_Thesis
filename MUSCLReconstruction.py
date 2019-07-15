import json         # read json format parameter file

# Monotonic Upwind Scheme for Conservation Laws
# As formulated in
#    Linear and Parabolic Reconstruction,
#    van Leer, B.
#    1979,
#    Towards the Ultimate Conservative Difference Scheme,
#    V. A Second Order Sequel to Godunov's Method,
#    J. Com. Phys,
#    32,
#    101-136
#
# Available on https://en.wikipedia.org/wiki/MUSCL_scheme


# 2nd Order
def muscl2(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_muscl[i] is equivalent to rho_{cell+i} for i={-1, 0, 1}
    rho_muscl = [variable_array[cell],       # center cell
                 variable_array[cell+1],     # right cell
                 variable_array[cell-1]]    # left cell

    # Avoid division by zero
    epsilon = 1e-6

    # STEP 2
    # Compute r_i
    r_i = (rho_muscl[0]-rho_muscl[-1])/(rho_muscl[1]-rho_muscl[0]+epsilon)

    # STEP 3
    # Calculate reconstructed values
    left_reco = rho_muscl[0]-0.5*limiter(r_i)*(rho_muscl[1]-rho_muscl[0])
    right_reco = rho_muscl[0] + 0.5 * limiter(r_i) * (rho_muscl[1] - rho_muscl[0])

    # End function
    return [left_reco, right_reco]


# 3rd Order
def muscl3(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_muscl[i] is equivalent to rho_{cell+i} for i={-1, 0, 1}
    rho_muscl = [variable_array[cell],       # center cell
                 variable_array[cell+1],     # right cell
                 variable_array[cell-1]]    # left cell

    # Avoid division by zero
    epsilon = 1e-6

    # STEP 2
    # Compute r_i
    r_i = (rho_muscl[0]-rho_muscl[-1])/(rho_muscl[1]-rho_muscl[0]+epsilon)

    # STEP 3
    # Compute delta values
    delta_r = rho_muscl[1]-rho_muscl[0]
    delta_l = rho_muscl[0]-rho_muscl[-1]

    # Constant
    kappa = 1/3

    # STEP 4a
    # Partial reconstructed values
    left_reco = (1-kappa)*delta_r+(1+kappa)*delta_l
    right_reco = (1-kappa)*delta_l+(1+kappa)*delta_r

    # STEP 4b
    # Final reconstructed values
    left_reco = rho_muscl[0] - 0.25*limiter(r_i)*left_reco
    right_reco = rho_muscl[0] + 0.25*limiter(r_i)*right_reco

    # End function
    return [left_reco, right_reco]


# Read parameter file
with open('params.txt') as file:
    params = json.load(file)

# Find required slope limiter
chosen_limiter = params['limiter']


# Slope Limiter
def limiter(r):
    # These limiters are given in https://en.wikipedia.org/wiki/Flux_limiter

    if chosen_limiter == "Charm":
        # Charm
        # Zhou, 1995
        # not 2nd order TDV accurate
        limited_slope = r*(3*r+1)/((r+1)**2) if r > 0 else 0

    elif chosen_limiter == "HCUS":
        # HCUS
        # Waterson & Deconinck, 1995
        # not 2nd order TDV accurate
        limited_slope = 1.5*(r+abs(r))/(r+2)

    elif chosen_limiter == "HQUICK":
        # HQUICK
        # Waterson & Deconinck, 1995
        # not 2nd order TDV accurate
        limited_slope = 2*(r+abs(r))/(r+3)

    elif chosen_limiter == "Koren":
        # Koren
        # Koren, 1993
        # 3rd order TDV accurate for sufficiently smooth data
        limited_slope = max(0, min(2*r, min((1+2*r)/3, 2)))

    elif chosen_limiter == "MinMod":
        # MinMod
        # Roe, 1986
        # 2nd order TDV accurate
        limited_slope = max(0, min(1, r))

    elif chosen_limiter == "MonotonizedCentral":
        # Monotonized Central
        # van Leer, 1977
        # 2nd order TDV accurate
        limited_slope = max(0, min(2*r, (1+r)/2, 2))

    elif chosen_limiter == "Osher":
        # Osher
        # Chakravarthy & Osher, 1983
        # 2nd order TDV accurate
        b = 1.5
        limited_slope = max(0, min(r, b))

    elif chosen_limiter == "Ospre":
        # Ospre
        # Waterson & Deconinck, 1995
        # 2nd order TDV accurate, symmetric
        limited_slope = 1.5*(r**2+r)/(r**2+r+1)

    elif chosen_limiter == "Smart":
        # Smart
        # Gaskell & Lau, 1988
        # not 2nd order TDV accurate
        limited_slope = max(0, min(2*r, (1+3*r)/4, 4))

    elif chosen_limiter == "Superbee":
        # Superbee
        # Roe, 1986
        # 2nd order TDV accurate, symmetric
        limited_slope = max(0, min(2*r, 1), min(r, 2))

    elif chosen_limiter == "Sweby":
        # Sweby
        # Sweby, 1984
        # not 2nd order TDV accurate, symmetric
        b = 1.5
        limited_slope = max(0, min(b*r, 1), min(r, b))

    elif chosen_limiter == "UMIST":
        # UMIST
        # Lien & Leschziner, 1994
        # 2nd order TDV accurate
        limited_slope = max(0, min(2*r, (1+3*r)/4, (3+r)/4, 2))

    elif chosen_limiter == "vanAlbada1":
        # van Albada 1
        # van Albada, et al., 1982
        # 2nd order TDV accurate, symmetric
        limited_slope = (r**2+r)/(r**2+1)

    elif chosen_limiter == "vanAlbada2":
        # van Albada 2
        # Kermani, 2003
        # not 2nd order TDV accurate
        # Alternate form for high spatial order schemes
        limited_slope = (2*r)/(r**2+1)

    elif chosen_limiter == "vanLeer":
        # van Leer
        # van Leer, 1974
        # 2nd order TDV accurate, symmetric
        limited_slope = (r+abs(r))/(1+abs(r))

    else:
        # Wrong slope limiter specifier

        print('ERROR: Wrong slope limiter specifier')
        exit(1)

    return limited_slope
