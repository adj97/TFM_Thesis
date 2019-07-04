import numpy as np  # numerical programming

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


# 5th Order
def weno5(cell, variable_array):

    # STEP 1
    # Construct required density values
    # rho_weno[i] is equivalent to rho_{cell+i} for i={-2, -1, 0, 1, 2}
    rho_weno = np.concatenate((variable_array[cell],                  # center cell
                               variable_array[cell+1:cell+3],         # right stencil cells
                               variable_array[cell-1:cell-3:-1]))     # left stencil cells

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
    for r in range(0,3):

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
    for r in range(0,3):

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
    # rho_weno[i] is equivalent to rho_{cell+i} for i={-2, -1, 0, 1, 2}
    rho_weno = np.concatenate((variable_array[cell:cell+3],           # right stencil cells
                               variable_array[cell-1:cell-3:-1]))     # left stencil cells

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
    for r in range(0,3):

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
    for r in range(0,3):

        # Right
        rhoplusWENO += omega[r]*rhoplusENO[r]

        # Left
        rhominWENO += omega_t[r]*rhominENO[r]

    # STEP 7
    # Assign reconstructed values to output objects
    left_reco = rhominWENO
    right_reco = rhoplusWENO


    # End func
    return [left_reco, right_reco]