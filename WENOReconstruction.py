import numpy as np  # numerical programming
import math         # mathematical functions

# Weighted Essentially Non-Oscillatory Reconstruction Schemes
# As formulated in
#    Procedure 2.2,
#    ENO and WENO Schemes for Hyperbolic Conservation Laws,
#    Chi-Wang Shu,
#    NASA/CR-97-206253,
#    ICASE Report No. 97-65,
#    November 1997
#
# Available on
#   https://www3.nd.edu/~zxu2/acms60790S13/Shu-WENO-notes.pdf


# 3rd Order
def weno3(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_weno[i] is equivalent to rho_{cell+i} for i={-2, -1, 0, 1, 2}
    rho_weno = [variable_array[cell],                  # center cell
                variable_array[cell+1],      # right
                variable_array[cell-1]]     # left

    # STEP 2
    # Define known ENO coefficients
    # For ~c = c[r-1][j]
    c = [[1/2, 1/2],      # r=0
         [-1/2, 3/2],      # r=1
         [3/2, -1/2]]     # r=-1  (for the ~c transformation)

    # Generate 3 right and 3 left polynomial ENO reconstructed values
    rhoplusENO = [0, 0]  # right
    rhominENO = [0, 0]   # left

    # Each ENO value
    for r in range(0, 2):

        # Initialise two sums
        poly_reco_val_plus = 0  # right
        poly_reco_val_min = 0   # left

        # Each local rho value and respective weight
        for j in range(0, 2):
            poly_reco_val_plus += c[r][j]*rho_weno[j-r]           # right
            poly_reco_val_min += c[r-1][j] * rho_weno[j - r]      # left

        # Store each left or right ENO value
        rhoplusENO[r] = poly_reco_val_plus  # right
        rhominENO[r] = poly_reco_val_min  # left

    # STEP 3
    # Define the known ENO polynomial convex weights
    # \tilde d[r] = d[1-r] with r={0, 1, 2}
    d = [2/3, 1/3]

    # STEP 4
    # Calculate the smoothness indicators
    beta = [0, 0]
    beta[0] = (rho_weno[1] - rho_weno[0]) ** 2
    beta[1] = (rho_weno[0] - rho_weno[-1]) ** 2

    # STEP 5a
    # Calculate 3 alpha weights
    # ~alpha = alpha * (d_{1-r}/d_r)
    alpha = [0, 0]      # right
    alpha_t = [0, 0]    # left

    # Avoid division by zero
    epsilon = 1e-6

    # Each alpha value
    for r in range(0, 2):

        # Right alpha
        alpha[r] = d[r]/((epsilon+beta[r])**2)

        # Left alpha (transformation)
        alpha_t[r] = alpha[r]*(d[1-r]/d[r])

    # STEP 5b
    # Calculate omega weights
    # Omega and ~omega separate arrays as no simple transformation
    omega = [0, 0]       # right
    omega_t = [0, 0]     # left

    # For each weight
    for r in range(0, 2):

        # right weights
        omega[r] = alpha[r]/sum(alpha)

        # left weights
        omega_t[r] = alpha_t[r]/sum(alpha_t)

    # STEP 6
    # Evaluate the two 5th order reconstructions
    # Initialise sums
    rhoplusWENO = 0   # right
    rhominWENO = 0    # left

    # Each weight and ENO reconstructed value in the sum
    for r in range(0, 2):

        # Right
        rhoplusWENO += omega[r]*rhoplusENO[r]

        # Left
        rhominWENO += omega_t[r]*rhominENO[r]

    # STEP 7
    # Assign reconstructed values to output objects
    left_reco = rhominWENO
    right_reco = rhoplusWENO

    # End function
    return [left_reco, right_reco]


# 5th Order
def weno5(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_weno[i] is equivalent to rho_{cell+i} for i={-2, -1, 0, 1, 2}
    rho_weno = np.concatenate((variable_array[cell:cell+3],           # center and right
                               variable_array[cell-1:cell-3:-1]))     # left

    # STEP 2
    # Define known ENO coefficients
    # For ~c = c[r-1][j]
    c = [[1/3, 5/6, -1/6],      # r=0
         [-1/6, 5/6, 1/3],      # r=1
         [1/3, -7/6, 11/6],     # r=2
         [11/6, -7/6, 1/3]]     # r=-1  (for the ~c transformation)

    # Generate 3 right and 3 left polynomial ENO reconstructed values
    rhoplusENO = [0, 0, 0]  # right
    rhominENO = [0, 0, 0]   # left

    # Each ENO value
    for r in range(0, 3):

        # Initialise two sums
        poly_reco_val_plus = 0  # right
        poly_reco_val_min = 0   # left

        # Each local rho value and respective weight
        for j in range(0, 3):
            poly_reco_val_plus += c[r][j]*rho_weno[j-r]           # right
            poly_reco_val_min += c[r-1][j] * rho_weno[j - r]      # left

        # Store each left or right ENO value
        rhoplusENO[r] = poly_reco_val_plus  # right
        rhominENO[r] = poly_reco_val_min  # left

    # STEP 3
    # Define the known ENO polynomial convex weights
    # \tilde d[r] = d[2-r] with r={0, 1, 2}
    d = [3/10, 3/5, 1/10]

    # STEP 4
    # Define smoothness sum coefficients
    coeff1 = [1, -2, 1]     # r-invariant
    coeff2 = [[3, -4, 1],   # r=0
              [1, 0, -1],   # r=1
              [1, -4, 3]]   # r=2

    # Calculate the smoothness indicators
    beta = [0, 0, 0]

    # Each smoothness indicator
    for r in range(0, 3):

        # Initialise sums
        sum1 = 0
        sum2 = 0

        # Each new term
        for j in range(0, 3):
            sum1 += coeff1[j]*rho_weno[j-r]
            sum2 += coeff2[r][j]*rho_weno[j-r]

        # Smoothness indicator form
        beta[r] = (13/12)*(sum1**2)+(1/4)*(sum2**2)

    # STEP 5a
    # Calculate 3 alpha weights
    # ~alpha = alpha * (d_{2-r}/d_r)
    alpha = [0, 0, 0]      # right
    alpha_t = [0, 0, 0]    # left

    # Avoid division by zero
    epsilon = 1e-6

    # Each alpha value
    for r in range(0, 3):

        # Right alpha
        alpha[r] = d[r]/((epsilon+beta[r])**2)

        # Left alpha (transformation)
        alpha_t[r] = alpha[r]*(d[2-r]/d[r])

    # STEP 5b
    # Calculate omega weights
    # Omega and ~omega separate arrays as no simple transformation
    omega = [0, 0, 0]       # right
    omega_t = [0, 0, 0]     # left

    # For each weight
    for r in range(0, 3):

        # right weights
        omega[r] = alpha[r]/sum(alpha)

        # left weights
        omega_t[r] = alpha_t[r]/sum(alpha_t)


    # STEP 6
    # Evaluate the two 5th order reconstructions
    # Initialise sums
    rhoplusWENO = 0   # right
    rhominWENO = 0    # left

    # Each weight and ENO reconstructed value in the sum
    for r in range(0, 3):

        # Right
        rhoplusWENO += omega[r]*rhoplusENO[r]

        # Left
        rhominWENO += omega_t[r]*rhominENO[r]

    # STEP 7
    # Assign reconstructed values to output objects
    left_reco = rhominWENO
    right_reco = rhoplusWENO

    # End function
    return [left_reco, right_reco]


# 7th Order
def weno7(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_weno[i] is equivalent to rho_{cell+i} for i={-3, -2, -1, 0, 1, 2, 3}
    rho_weno = np.concatenate((variable_array[cell:cell+4],           # right stencil cells
                               variable_array[cell-1:cell-4:-1]))     # left stencil cells

    # STEP 2
    # Define known ENO coefficients
    # For ~c = c[r-1][j]
    c = [[1/4, 13/12, -5/12, 1/12],       # r=0
         [-1/12, 7/12, 7/12, -1/12],      # r=1
         [1/12, -5/12, 13/12, 1/4],       # r=2
         [-1/4, 13/12, -23/12, 25/12],    # r=3
         [25/12, -23/12, 13/12, -1/4]]    # r=-1  (for the ~c transformation)

    # Generate 3 right and 3 left polynomial ENO reconstructed values
    rhoplusENO = [0, 0, 0, 0]  # right
    rhominENO = [0, 0, 0, 0]   # left

    # Each ENO value
    for r in range(0, 4):

        # Initialise two sums
        poly_reco_val_plus = 0  # right
        poly_reco_val_min = 0   # left

        # Each local rho value and respective weight
        for j in range(0, 4):
            poly_reco_val_plus += c[r][j] * rho_weno[j-r]           # right
            poly_reco_val_min += c[r-1][j] * rho_weno[j - r]      # left

        # Store each left or right ENO value
        rhoplusENO[r] = poly_reco_val_plus  # right
        rhominENO[r] = poly_reco_val_min  # left

    # STEP 3
    # Define the known ENO polynomial convex weights
    # \tilde d[r] = d[2-r] with r={0, 1, 2}
    d = [4 / 35, 18 / 35, 12 / 35, 1 / 35]

    # STEP 4
    # Define smoothness sum coefficients
    # b[r,j,l] three indices
    # r = 0, 1, 2, 3
    # j = 0, 1, 2, 3
    # l = 0, ..., 3-j
    b = [[[2107, -9402, 7042, -1854], [11003, -17246, 4642], [7043, -3882], [547]],  # r=0
         [[547, -2522, 1922, -494], [3443, -5966, 1602], [2843, -1642], [267]],      # r=1
         [[267, -1642, 1602, -494], [2843, -5966, 1922], [3443, -2522], [547]],      # r=2
         [[547, -3882, 4642, -1854], [7043, -17246, 7042], [11003, -9402], [2107]]]  # r=3

    # Calculate the smoothness indicators
    beta = [0, 0, 0, 0]

    # Each smoothness indicator
    for r in range(0, 4):

        # Initialise outer sum
        outer_sum = 0

        # J loop
        for j in range(0, 4):

            # Initialise inner sum
            inner_sum = 0

            # L loop
            for l in range(0, 3-j):
                inner_sum += b[r][j][l]*rho_weno[j+l-r]

            # Multiply inner sum
            inner_sum *= rho_weno[j-r]

            # Add to outer sum
            outer_sum += inner_sum

        # Allocate smoothness indicator
        beta[r] = outer_sum

    # STEP 5a
    # Calculate 3 alpha weights
    # ~alpha = alpha * (d_{2-r}/d_r)
    alpha = [0, 0, 0, 0]      # right
    alpha_t = [0, 0, 0, 0]    # left

    # Avoid division by zero
    epsilon = 1e-6

    # Each alpha value
    for r in range(0, 4):

        # Right alpha
        alpha[r] = d[r]/((epsilon+beta[r])**2)

        # Left alpha (transformation)
        alpha_t[r] = alpha[r]*(d[3-r]/d[r])

    # STEP 5b
    # Calculate omega weights
    # Omega and ~omega separate arrays as no simple transformation
    omega = [0, 0, 0, 0]       # right
    omega_t = [0, 0, 0, 0]     # left

    # For each weight
    for r in range(0, 4):

        # right weights
        omega[r] = alpha[r]/sum(alpha)

        # left weights
        omega_t[r] = alpha_t[r]/sum(alpha_t)

    # STEP 6
    # Evaluate the two 5th order reconstructions
    # Initialise sums
    rhoplusWENO = 0   # right
    rhominWENO = 0    # left

    # Each weight and ENO reconstructed value in the sum
    for r in range(0, 4):

        # Right
        rhoplusWENO += omega[r]*rhoplusENO[r]

        # Left
        rhominWENO += omega_t[r]*rhominENO[r]

    # STEP 7
    # Monotonicity preserving bounds
    # A.  Suresh  and  H.  T.  Huynh,
    # Accurate monotonicity preserving scheme with Runge-Kutta time-stepping,
    # J. Comput. Phys.136, 83 (1997).

    # STEP 7a
    # Zone center curvature measures
    d_j = rho_weno[1]-2*rho_weno[0]+rho_weno[-1]
    d_jp1 = rho_weno[2]-2*rho_weno[1]+rho_weno[0]
    d_jm1 = rho_weno[0]-2*rho_weno[-1]+rho_weno[-2]

    # STEP 7b
    # Minmod of local curvature
    right_dmd = minmod(d_j, d_jp1)
    right_dlc = minmod(d_j, d_jp1)
    left_dmd = minmod(d_j, d_jm1)
    left_dlc = minmod(d_j, d_jm1)

    # STEP 7c
    # Define curvature constants
    alpha_curv = 2
    beta_curv = 4

    # STEP 7d
    # Median bounds
    right_md = 0.5*(rho_weno[0]+rho_weno[1])-0.5*right_dmd
    left_md = 0.5*(rho_weno[0]+rho_weno[-1])-0.5*left_dmd
    # Upper limit bound
    right_ul = rho_weno[0] + alpha_curv*(rho_weno[0]-rho_weno[-1])
    left_ul = rho_weno[0] + alpha_curv*(rho_weno[0]-rho_weno[1])
    # Large curvature
    right_lc = rho_weno[0] + 0.5*(rho_weno[0]-rho_weno[-1])+beta_curv*right_dlc/3
    left_lc = rho_weno[0] + 0.5*(rho_weno[0]-rho_weno[1])+beta_curv*left_dlc/3

    # STEP 7e
    # Define minimum and maximum left and right bounds
    right_min = max(min(rho_weno[0], rho_weno[1], right_md), min(rho_weno[0], right_ul, right_lc))
    right_max = min(max(rho_weno[0], rho_weno[1], right_md), max(rho_weno[0], right_ul, right_lc))
    left_min = max(min(rho_weno[0], rho_weno[-1], left_md), min(rho_weno[0], left_ul, left_lc))
    left_max = min(max(rho_weno[0], rho_weno[-1], left_md), max(rho_weno[0], left_ul, left_lc))

    # STEP 8
    # Allocate the recosntructed monotonic bounded values
    right_reco = median(rhoplusWENO, right_min, right_max)
    left_reco = median(rhominWENO, left_min, left_max)

    # End func
    return [left_reco, right_reco]


# Minmod function
def minmod(arg1,arg2):
    return 0.5*(math.copysign(1.0, arg1)+math.copysign(1.0, arg2))*min(abs(arg1), abs(arg2))


# Median function
def median(x,y,z):
    return x+minmod(y-x, z-x)
