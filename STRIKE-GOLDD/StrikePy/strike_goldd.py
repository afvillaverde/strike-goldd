"""
=================================================================================================================
 -> STRIKE-GOLDD (Structural Identifiability taken as Extended-Generalized Observability with Lie Derivatives) <-
-----------------------------------------------------------------------------------------------------------------
  * A Matlab toolbox implementation in Python for structural identifiability and observability analysis 
    of nonlinear models by David Rey Rostro (davidreyrostro@gmail.com)
  * Based on the Matlab STRIKE-GOLDD version by Alejandro Fernandez Villaverde (afvillaverde@uvigo.gal)
-----------------------------------------------------------------------------------------------------------------
  * User input can be specified in file 'options.py' or by creating an options file in the folder custom_options
=================================================================================================================
"""

import os
from sympy.matrices import Matrix, zeros
from time import time
from math import ceil
from functions.rationalize import rationalize_all_numbers
from datetime import datetime
from pathlib import Path

def strike_goldd(*args):

    ########################################################################
    # StrikePy directory...
    dir = os.path.dirname(os.path.abspath(__file__))
    ########################################################################
    # Create 'results' folder in StrikePy directory in case it does not exist
    path = Path('{}/results'.format(dir))
    path.mkdir(parents=True, exist_ok=True)
    ########################################################################
    print('\n\n ---------------------------------------')
    print(' >>> STRIKE-GOLDD toolbox for Python <<<')
    print(' --------------------------------------- \n')

    tStart = time()
    ########################################################################
    # Read options and load model
    import numpy as np
    import sympy as sym
    import importlib
    from functions.elim_and_recalc import elim_and_recalc
    import symbtools as st

    # If no parameters are passed when calling the function the default options file is used
    if len(args) == 0:
        import options
        model = importlib.import_module('models.{}'.format(options.modelname))
        print(' Analyzing the {} model... \n'.format(options.modelname))

    # If the name of a file located in the custom_options folder is passed as a parameter to the function, this file is used as the options for the analysis.
    if len(args) == 1:
        options = importlib.import_module('custom_options.{}'.format(args[0]))
        model = importlib.import_module('models.{}'.format(options.modelname))
        print(' Analyzing the {} model... \n'.format(options.modelname))

    tic = time()
    ########################################################################
    # Initialize variables:
    identifiables   = []    # identifiable parameters.
    nonidentif      = []    # unidentifiable parameters.
    obs_states      = []    # observable states.
    unobs_states    = []    # unobservable states.
    obs_inputs      = []    # observable inputs.
    unobs_inputs    = []    # unobservable inputs.
    lastrank        = None
    unidflag        = 0
    skip_elim       = 0
    isFISPO         = 0
    ########################################################################
    # Remove unknown parameters that have already been classified as identifiable:
    if len(model.p) > 1:
        try:
            len(model.p[0])
        except:
            auxiliar = []
            for i in range(len(model.p)):
                auxiliar.append([model.p[i]])
            model.p = auxiliar

    try:
        options.prev_ident_pars
    except:
        options.prev_ident_pars = []

    if len(options.prev_ident_pars) != 0:
        if len(model.p) == 1:
            for npi in range(len(options.prev_ident_pars)):
                if options.prev_ident_pars[npi] in model.p:
                    model.p = []
        if len(model.p) > 1:
            for npi in range(len(options.prev_ident_pars)):
                if any(options.prev_ident_pars[npi] in arr for arr in model.p):
                    known = [options.prev_ident_pars[npi]]
                    model.p.remove(known)
    ########################################################################
    # Dimensions of the problem:

    m = len(model.h)  # number of outputs
    n = len(model.x)  # number of states
    q = len(model.p)  # number of unknown parameters

    if len(model.h) > 1:
        try:
            len(model.h[0])
        except:
            auxiliar = []
            for i in range(m):
                auxiliar.append([model.h[i]])
            model.h = auxiliar

    if len(model.x) > 1:
        try:
            len(model.x[0])
        except:
            auxiliar = []
            for i in range(n):
                auxiliar.append([model.x[i]])
            model.x = auxiliar

    try:
        model.w
    except:
        model.w = []

    # number of unknown inputs
    nw = len(model.w)

    try:
        model.u
    except:
        model.u = []

    # number of known inputs
    nu = len(model.u)

    print(" >>> The model contains:\n{} states:\n {}".format(n,np.array(model.x)))
    print("{} outputs:\n {}".format(m,np.array(model.h)))
    print("{} known inputs:\n {}".format(nu,np.array(model.u)))
    print("{} unknown inputs:\n {}".format(nw,np.array(model.w)))
    print("{} parameters:\n {}".format(q,np.array(model.p)))
    ########################################################################
    # Check which states are directly measured, if any.
    # Basically it is checked if any state is directly on the output, then that state is directly measurable.
    if m == 1:
        saidas = model.h
    else:
        saidas = []
        for i in range(m):
            saidas.append(model.h[i])

    if n == 1:
        estados = model.x
    else:
        estados = []
        for i in range(n):
            estados.append(model.x[i][0])

    ismeasured = []
    for i in range(n):
        ismeasured.append(0) # creo lista de n ceros

    if len(saidas) == 1:
        for i in range(n):
            if estados[i] in saidas:
                ismeasured[i] = 1
    else:
        for i in range(n):
            if any(estados[i] in arr for arr in saidas):
                ismeasured[i] = 1

    meas_x_indices = [] # indices of the measured states
    for i in range(n):
        if ismeasured[i] == 1:
            meas_x_indices.append(i)

    unmeas_x_indices = [] # indices of the unmeasured states
    for i in range(n):
        if ismeasured[i] == 0:
            unmeas_x_indices.append(i)

    meas_x = [] # names of the measured states
    if len(meas_x_indices) == 1 and n == 1:
        meas_x = estados
    if len(meas_x_indices) == 1 and n != 1:
        meas_x.append(estados[meas_x_indices[0]])
    if len(meas_x_indices) > 1:
        for i in range(len(meas_x_indices)):
            meas_x.append([estados[meas_x_indices[i]]])
    ########################################################################
    r = n + q + nw  # number of unknown variables to observe / identify
    nd = ceil((r - m)/m)  # minimum number of Lie derivatives for Oi to have full rank
    print('\n\n >>> Building the observability-identifiability matrix requires at least {} Lie derivatives'.format(nd))
    print('     Calculating derivatives: ', end="")

    tic = time()
    ########################################################################
    # Check if the size of nnzDerU and nnzDerW are appropriate
    if len(model.u) > len(options.nnzDerU):
        raise Exception(' The number of known inputs is higher than the size of nnzDerU and must have the same size. \n            Go to the options file and modify it. \n            For more information about the error see point 7 of the StrikePy instruction manual. ')
    if len(model.w) > len(options.nnzDerW):
        raise Exception(' The number of unknown inputs is higher than the size of nnzDerW and must have the same size. \n            Go to the options file and modify it. \n            For more information about the error see point 7 of the StrikePy instruction manual. ')
    ########################################################################
    # Input derivates:

    # Create array of known inputs and set certain derivatives to zero:
    input_der = []
    if len(model.u) > 0:
        for ind_u in range(len(model.u)):  # create array of derivatives of the inputs
            if len(model.u) == 1:
                locals()['{}'.format(model.u[ind_u])] = sym.Symbol('{}'.format(model.u[ind_u]))  # the first element is the underived input
                auxiliar = [locals()['{}'.format(model.u[ind_u])]]
            else:
                locals()['{}'.format(model.u[ind_u][0])] = sym.Symbol('{}'.format(model.u[ind_u][0]))  # the first element is the underived input
                auxiliar = [locals()['{}'.format(model.u[ind_u][0])]]
            for k in range(nd):
                if len(model.u) == 1:
                    locals()['{}_d{}'.format(model.u[ind_u], k + 1)] = sym.Symbol('{}_d{}'.format(model.u[ind_u], k + 1))
                    auxiliar.append(locals()['{}_d{}'.format(model.u[ind_u], k + 1)])
                else:
                    locals()['{}_d{}'.format(model.u[ind_u][0], k + 1)] = sym.Symbol('{}_d{}'.format(model.u[ind_u][0], k + 1))
                    auxiliar.append(locals()['{}_d{}'.format(model.u[ind_u][0], k + 1)])
            if len(model.u) == 1:
                input_der = auxiliar
                if len(input_der) >= options.nnzDerU[0] + 1:
                    for i in range(len(input_der[(options.nnzDerU[0] + 1):])):
                        input_der[(options.nnzDerU[0] + 1) + i] = 0
            else:
                input_der.append(auxiliar)
                if len(input_der[0]) >= options.nnzDerU[ind_u] + 1:
                    for i in range(len(input_der[0][(options.nnzDerU[ind_u] + 1):])):
                        input_der[ind_u][(options.nnzDerU[ind_u] + 1) + i] = 0
    zero_input_der_dummy_name = sym.Symbol('zero_input_der_dummy_name')

    # Create array of unknown inputs and set certain derivatives to zero:
    w_der = []
    if len(model.w) > 0:
        for ind_w in range(len(model.w)):  # create array of derivatives of the inputs
            if len(model.w) == 1:
                locals()['{}'.format(model.w[ind_w])] = sym.Symbol('{}'.format(model.w[ind_w]))  # the first element is the underived input
                auxiliar = [locals()['{}'.format(model.w[ind_w])]]
            else:
                locals()['{}'.format(model.w[ind_w][0])] = sym.Symbol('{}'.format(model.w[ind_w][0]))  # the first element is the underived input
                auxiliar = [locals()['{}'.format(model.w[ind_w][0])]]
            for k in range(nd + 1):
                if len(model.w) == 1:
                    locals()['{}_d{}'.format(model.w[ind_w], k + 1)] = sym.Symbol('{}_d{}'.format(model.w[ind_w], k + 1))
                    auxiliar.append(locals()['{}_d{}'.format(model.w[ind_w], k + 1)])
                else:
                    locals()['{}_d{}'.format(model.w[ind_w][0], k + 1)] = sym.Symbol('{}_d{}'.format(model.w[ind_w][0], k + 1))
                    auxiliar.append(locals()['{}_d{}'.format(model.w[ind_w][0], k + 1)])
            if len(model.w) == 1:
                w_der = auxiliar
                if len(w_der) >= options.nnzDerW[0] + 1:
                    for i in range(len(w_der[(options.nnzDerW[0] + 1):])):
                        w_der[(options.nnzDerW[0] + 1) + i] = 0
            else:
                w_der.append(auxiliar)
                if len(w_der[0]) >= options.nnzDerW[ind_w] + 1:
                    for i in range(len(w_der[0][(options.nnzDerW[ind_w] + 1):])):
                        w_der[ind_w][(options.nnzDerW[ind_w] + 1) + i] = 0

        if sym.shape(Matrix(w_der).T)[0] == 1:
            w1vector = []
            for i in range(len(w_der) -1):
                w1vector.append([w_der[i]])
            w1vector_dot = []
            for i in range(len(w_der)):
                if i !=0:
                    w1vector_dot.append([w_der[i]])

        else:
            w1vector = []
            for k in range(sym.shape(Matrix(w_der))[1] - 1):
                for i in w_der:
                    w1vector.append([i[k]])
            w1vector_dot = []
            for k in range(sym.shape(Matrix(w_der))[1]):
                for i in w_der:
                    if k != 0:
                        w1vector_dot.append([i[k]])

        # -- Include as states only nonzero inputs / derivatives:
        nzi = []
        for fila in range(len(w1vector)):
            if w1vector[fila][0] != 0:
                nzi.append([fila])
        nzj = []
        for fila in range(len(w1vector)):
            if w1vector[fila][0] != 0:
                nzj.append([1])
        nz_w1vec = []
        for fila in range(len(w1vector)):
            if w1vector[fila][0] != 0:
                nz_w1vec.append(w1vector[fila])
        w1vector = nz_w1vec
        w1vector_dot = w1vector_dot[0:len(nzi)]

    else:
        w1vector = []
        w1vector_dot = []

    ########################################################################
    # Augment state vector, dynamics:
    if len(model.x) == 1:
        xaug = []
        xaug.append(model.x)
        xaug = np.append(xaug, model.p, axis=0)
        if len(w1vector) != 0:
            xaug = np.append(xaug, w1vector, axis=0)

        faug = []
        faug.append(model.f)
        faug = np.append(faug, zeros(len(model.p), 1), axis=0)
        if len(w1vector) != 0:
            faug = np.append(faug, w1vector_dot, axis=0)

    else:
        xaug = model.x
        xaug = np.append(xaug, model.p, axis=0)
        if len(w1vector) != 0:
            xaug = np.append(xaug, w1vector, axis=0)

        faug = model.f
        faug = np.append(faug, zeros(len(model.p), 1), axis=0)
        if len(w1vector) != 0:
            faug = np.append(faug, w1vector_dot, axis=0)
    ########################################################################
    # Build Oi:
    onx = np.array(zeros(m * (1 + nd), n + q + len(w1vector)))
    jacobiano = sym.Matrix(model.h).jacobian(xaug)
    onx[0:len(model.h)] = np.array(jacobiano)  # first row(s) of onx (derivative of the output with respect to the vector states+unknown parameters).
    toc = time()
    totaltime = toc - tic  # execution time
    ind = 0  # Lie derivative index (sometimes called 'k')
    lasttime = 0
    ########################################################################
    past_Lie = model.h
    extra_term = np.array(0)

    # loop as long as I don't complete the preset Lie derivatives or go over the maximum time set for each derivative
    while ind < nd and lasttime < options.maxLietime:
        tic = time()
        Lieh = Matrix((onx[(ind * m):(ind + 1) * m][:]).dot(faug))
        if ind > 0:
            if len(model.u) > 0:
                for i in range(ind):
                    if len(model.u) == 1:
                        column = len(input_der) - 1
                        if i < column:
                            lo_u_der = input_der[i]
                            if lo_u_der == 0:
                                lo_u_der = zero_input_der_dummy_name
                            lo_u_der = np.array([lo_u_der])
                            hi_u_der = input_der[i + 1]
                            hi_u_der = Matrix([hi_u_der])

                            intermedio = sym.Matrix([past_Lie]).jacobian(lo_u_der) * hi_u_der
                            if extra_term:
                                extra_term = extra_term + intermedio
                            else:
                                extra_term = intermedio
                    else:
                        column = len(input_der[0]) - 1
                        if i < column:
                            lo_u_der = []
                            hi_u_der = []
                            for fila in input_der:
                                lo_u_der.append(fila[i])
                                hi_u_der.append(fila[i + 1])
                            for i in range(len(lo_u_der)):
                                if lo_u_der[i] == 0:
                                    lo_u_der[i] = zero_input_der_dummy_name
                            lo_u_der = np.array(lo_u_der)
                            hi_u_der = Matrix(hi_u_der)
                            intermedio = sym.Matrix([past_Lie]).jacobian(lo_u_der) * hi_u_der
                            if extra_term:
                                extra_term = extra_term + intermedio
                            else:
                                extra_term = intermedio
        if extra_term:
            ext_Lie = Lieh + extra_term
        else:
            ext_Lie = Lieh
        past_Lie = ext_Lie
        onx[((ind + 1) * m):(ind + 2) * m] = sym.Matrix(ext_Lie).jacobian(xaug)
        lasttime = time() - tic
        totaltime = totaltime + lasttime
        ind = ind + 1
        print(end=' {}'.format(ind))

    if ind == nd: #If I have done all the minimum derivatives to build onx (I have not exceeded the time limit)....
        increaseLie = 1
        while increaseLie == 1:  # while increaseLie is 1 I will increase the size of onx
            print('\n >>> Observability-Identifiability matrix built with {} Lie derivatives'.format(nd))
            print('     (calculated in {} seconds)'.format(totaltime))
            # =============================================================================================
            # The observability/identifiability matrix is saved in a .txt file
            file = open("results/obs_ident_matrix_{}_{}_Lie_deriv.txt".format(options.modelname, nd), "w")
            file.write('onx = {}'.format(str(onx.tolist())))
            file.close()
            # =============================================================================================
            # Check identifiability by calculating rank:
            print(' >>> Calculating rank of matrix with size {}x{}...'.format(sym.shape(Matrix(onx))[0], sym.shape(Matrix(onx))[1]))
            tic = time()
            rational_onx = rationalize_all_numbers(Matrix(onx))
            rango = st.generic_rank(Matrix(rational_onx))
            toc = time() - tic
            print('     Rank = {} (calculated in {} seconds)'.format(rango, toc))
            if rango == len(xaug): # If the onx matrix already has full rank... all is observable and identifiable
                obs_states = model.x
                obs_inputs = model.w
                identifiables = model.p
                increaseLie = 0 # stop increasing the number of onx rows with derivatives

            else: # With that number of Lie derivatives the array is not full rank.
                #----------------------------------------------------------
                # If there are unknown inputs, we may want to check id/obs of (x,p,w) and not of dw/dt:
                if len(model.w) > 0:
                    [identifiables, nonidentif, obs_states, unobs_states, obs_inputs, unobs_inputs] = elim_and_recalc(unmeas_x_indices, rango, onx, model.p, model.x, unidflag, w1vector, identifiables, obs_states, obs_inputs)

                    # Check which unknown inputs are observable:
                    obs_in_no_der = []
                    if len(model.w) == 1 and len(obs_inputs) > 0:
                        if model.w == obs_inputs:
                            obs_in_no_der = model.w
                    if len(model.w) > 1 and len(obs_inputs) > 0:
                        for elemento in model.w:
                            if len(obs_inputs) == 1:
                                if elemento == obs_inputs:
                                    obs_in_no_der = elemento
                            else:
                                for input in obs_inputs:
                                    if elemento == input:
                                        obs_in_no_der.append(elemento[0])
                    if len(identifiables) == len(model.p) and len(obs_states)+len(meas_x) == len(model.x) and len(obs_in_no_der) == len(model.w):
                        obs_states = model.x
                        obs_inputs = obs_in_no_der
                        identifiables = model.p
                        increaseLie = 0  # -> with this we skip the next 'if' block and jump to the end of the algorithm
                        isFISPO = 1
                #----------------------------------------------------------
                # If possible (& necessary), calculate one more Lie derivative and retry:
                if nd < len(xaug) and lasttime < options.maxLietime and rango != lastrank and increaseLie == 1:
                    tic = time()
                    ind = nd
                    nd = nd + 1 # One is added to the number of derivatives already made
                    extra_term = np.array(0) # reset for each new Lie derivative
                    # - Known input derivatives: ----------------------------------
                    if len(model.u) > 0: # Extra terms of extended Lie derivatives
                        # may have to add extra input derivatives (note that 'nd' has grown):
                        input_der = []
                        for ind_u in range(len(model.u)):  # create array of derivatives of the inputs
                            if len(model.u) == 1:
                                locals()['{}'.format(model.u[ind_u])] = sym.Symbol('{}'.format(model.u[ind_u]))  # the first element is the underived input
                                auxiliar = [locals()['{}'.format(model.u[ind_u])]]
                            else:
                                locals()['{}'.format(model.u[ind_u][0])] = sym.Symbol('{}'.format(model.u[ind_u][0]))  # the first element is the underived input
                                auxiliar = [locals()['{}'.format(model.u[ind_u][0])]]
                            for k in range(nd):
                                if len(model.u) == 1:
                                    locals()['{}_d{}'.format(model.u[ind_u], k + 1)] = sym.Symbol('{}_d{}'.format(model.u[ind_u], k + 1))
                                    auxiliar.append(locals()['{}_d{}'.format(model.u[ind_u], k + 1)])
                                else:
                                    locals()['{}_d{}'.format(model.u[ind_u][0], k + 1)] = sym.Symbol('{}_d{}'.format(model.u[ind_u][0], k + 1))
                                    auxiliar.append(locals()['{}_d{}'.format(model.u[ind_u][0], k + 1)])
                            if len(model.u) == 1:
                                input_der = auxiliar
                                if len(input_der) >= options.nnzDerU[0] + 1:
                                    for i in range(len(input_der[(options.nnzDerU[0] + 1):])):
                                        input_der[(options.nnzDerU[0] + 1) + i] = 0
                            else:
                                input_der.append(auxiliar)
                                if len(input_der[0]) >= options.nnzDerU[ind_u] + 1:
                                    for i in range(len(input_der[0][(options.nnzDerU[ind_u] + 1):])):
                                        input_der[ind_u][(options.nnzDerU[ind_u] + 1) + i] = 0

                        for i in range(ind):
                            if len(model.u) == 1:
                                column = len(input_der) - 1
                                if i < column:
                                    lo_u_der = input_der[i]
                                    if lo_u_der == 0:
                                        lo_u_der = zero_input_der_dummy_name
                                    lo_u_der = np.array([lo_u_der])
                                    hi_u_der = input_der[i + 1]
                                    hi_u_der = Matrix([hi_u_der])

                                    intermedio = sym.Matrix([past_Lie]).jacobian(lo_u_der) * hi_u_der
                                    if extra_term:
                                        extra_term = extra_term + intermedio
                                    else:
                                        extra_term = intermedio
                            else:
                                column = len(input_der[0]) - 1
                                if i < column:
                                    lo_u_der = []
                                    hi_u_der = []
                                    for fila in input_der:
                                        lo_u_der.append(fila[i])
                                        hi_u_der.append(fila[i + 1])
                                    for i in range(len(lo_u_der)):
                                        if lo_u_der[i] == 0:
                                            lo_u_der[i] = zero_input_der_dummy_name
                                    lo_u_der = np.array(lo_u_der)
                                    hi_u_der = Matrix(hi_u_der)
                                    intermedio = sym.Matrix([past_Lie]).jacobian(lo_u_der) * hi_u_der
                                    if extra_term:
                                        extra_term = extra_term + intermedio
                                    else:
                                        extra_term = intermedio

                    #- Unknown input derivatives:----------------
                    # add new derivatives, if they are not zero:
                    if len(model.w) > 0:
                        prev_size = len(w1vector)
                        w_der = []
                        for ind_w in range(len(model.w)):  # create array of derivatives of the inputs
                            if len(model.w) == 1:
                                locals()['{}'.format(model.w[ind_w])] = sym.Symbol(
                                    '{}'.format(model.w[ind_w]))  # the first element is the underived input
                                auxiliar = [locals()['{}'.format(
                                    model.w[ind_w])]]
                            else:
                                locals()['{}'.format(model.w[ind_w][0])] = sym.Symbol('{}'.format(
                                    model.w[ind_w][0]))  # the first element is the underived input
                                auxiliar = [locals()['{}'.format(
                                    model.w[ind_w][0])]]
                            for k in range(nd + 1):
                                if len(model.w) == 1:
                                    locals()['{}_d{}'.format(model.w[ind_w], k + 1)] = sym.Symbol(
                                        '{}_d{}'.format(model.w[ind_w], k + 1))
                                    auxiliar.append(locals()['{}_d{}'.format(model.w[ind_w], k + 1)])
                                else:
                                    locals()['{}_d{}'.format(model.w[ind_w][0], k + 1)] = sym.Symbol(
                                        '{}_d{}'.format(model.w[ind_w][0], k + 1))
                                    auxiliar.append(locals()['{}_d{}'.format(model.w[ind_w][0], k + 1)])
                            if len(model.w) == 1:
                                w_der = auxiliar
                                if len(w_der) >= options.nnzDerW[0] + 1:
                                    for i in range(len(w_der[(options.nnzDerW[0] + 1):])):
                                        w_der[(options.nnzDerW[0] + 1) + i] = 0
                            else:
                                w_der.append(auxiliar)
                                if len(w_der[0]) >= options.nnzDerW[ind_w] + 1:
                                    for i in range(len(w_der[0][(options.nnzDerW[ind_w] + 1):])):
                                        w_der[ind_w][(options.nnzDerW[ind_w] + 1) + i] = 0

                        if sym.shape(Matrix(w_der).T)[0] == 1:
                            w1vector = []
                            for i in range(len(w_der) - 1):
                                w1vector.append([w_der[i]])
                            w1vector_dot = []
                            for i in range(len(w_der)):
                                if i != 0:
                                    w1vector_dot.append([w_der[i]])

                        else:
                            w1vector = []
                            for k in range(sym.shape(Matrix(w_der))[1] - 1):
                                for i in w_der:
                                    w1vector.append([i[k]])
                            w1vector_dot = []
                            for k in range(sym.shape(Matrix(w_der))[1]):
                                for i in w_der:
                                    if k != 0:
                                        w1vector_dot.append([i[k]])

                        # -- Include as states only nonzero inputs / derivatives:
                        nzi = []
                        for fila in range(len(w1vector)):
                            if w1vector[fila][0] != 0:
                                nzi.append([fila])
                        nzj = []
                        for fila in range(len(w1vector)):
                            if w1vector[fila][0] != 0:
                                nzj.append([1])
                        nz_w1vec = []
                        for fila in range(len(w1vector)):
                            if w1vector[fila][0] != 0:
                                nz_w1vec.append(w1vector[fila])
                        w1vector = nz_w1vec
                        w1vector_dot = w1vector_dot[0:len(nzi)]

                        ########################################################################
                        # Augment state vector, dynamics:
                        if len(model.x) == 1:
                            xaug = []
                            xaug.append(model.x)
                            xaug = np.append(xaug, model.p, axis=0)
                            if len(w1vector) != 0:
                                xaug = np.append(xaug, w1vector, axis=0)

                            faug = []
                            faug.append(model.f)
                            faug = np.append(faug, zeros(len(model.p), 1), axis=0)
                            if len(w1vector) != 0:
                                faug = np.append(faug, w1vector_dot, axis=0)

                        else:
                            xaug = model.x
                            xaug = np.append(xaug, model.p, axis=0)
                            if len(w1vector) != 0:
                                xaug = np.append(xaug, w1vector, axis=0)

                            faug = model.f
                            faug = np.append(faug, zeros(len(model.p), 1), axis=0)
                            if len(w1vector) != 0:
                                faug = np.append(faug, w1vector_dot, axis=0)
                        ########################################################################
                        # -- Augment size of the Obs-Id matrix if needed:
                        new_size = len(w1vector)
                        onx = np.append(onx, zeros((ind + 1)*m, new_size-prev_size), axis=1)
                    ########################################################################
                    newLie = Matrix((onx[(ind * m):(ind + 1) * m][:]).dot(faug))
                    if extra_term:
                        past_Lie = newLie + extra_term
                    else:
                        past_Lie = newLie
                    newOnx = sym.Matrix(past_Lie).jacobian(xaug)
                    onx = np.append(onx, newOnx, axis=0)
                    lasttime = time() - tic
                    totaltime = totaltime + lasttime
                    lastrank = rango

                # If that is not possible, there are several possible causes:
                # This is the case when you have onx with all possible derivatives done and it is not full rank, the maximum time for the next derivative has passed
                # or the matrix no longer increases in rank as derivatives are increased.
                else:
                    if nd >= len(xaug):  # The maximum number of Lie derivatives has been reached
                        unidflag = 1
                        print('\n >>> The model is structurally unidentifiable as a whole')
                    else:
                        if rango == lastrank:
                            onx = onx[0:(-1 -(m - 1))]
                            nd = nd - 1 # It is indicated that the number of derivatives needed was one less than the number of derivatives made
                            unidflag = 1
                        else:
                            if lasttime >= options.maxLietime:
                                print('\n => More Lie derivatives would be needed to see if the model is structurally unidentifiable as a whole.')
                                print('    However, the maximum computation time allowed for calculating each of them has been reached.')
                                print('    You can increase it by changing <<maxLietime>> in options (currently maxLietime = {})'.format(options.maxLietime))
                                unidflag = 0
                    if skip_elim == 0 and isFISPO == 0:
                        # Eliminate columns one by one to check identifiability of the associated parameters:
                        [identifiables, nonidentif, obs_states, unobs_states, obs_inputs, unobs_inputs] = elim_and_recalc(unmeas_x_indices, rango, onx, model.p, model.x, unidflag, w1vector, identifiables, obs_states, obs_inputs)

                        # Check which unknown inputs are observable:
                        obs_in_no_der = []
                        if len(model.w) == 1 and len(obs_inputs) > 0:
                            if model.w == obs_inputs:
                                obs_in_no_der = model.w
                        if len(model.w) > 1 and len(obs_inputs) > 0:
                            for elemento in model.w:  # for each unknown input
                                if len(obs_inputs) == 1:
                                    if elemento == obs_inputs:
                                        obs_in_no_der = elemento
                                else:
                                    for input in obs_inputs:
                                        if elemento == input:
                                            obs_in_no_der.append(elemento[0])

                        if len(identifiables) == len(model.p) and (len(obs_states) + len(meas_x)) == len(model.x) and len(obs_in_no_der) == len(model.w):
                            obs_states = model.x
                            obs_inputs = obs_in_no_der
                            identifiables = model.p
                            increaseLie = 0 # -> with this we skip the next 'if' block and jump to the end of the algorithm
                            isFISPO = 1
                        increaseLie = 0

    else: # If the maxLietime has been reached, but the minimum of Lie derivatives has not been calculated:
        print('\n => More Lie derivatives would be needed to analyse the model.')
        print('    However, the maximum computation time allowed for calculating each of them has been reached.')
        print('    You can increase it by changing <<maxLietime>> in options (currently maxLietime = {})'.format(options.maxLietime))
        tic = time()
        print('\n >>> Calculating rank of matrix with size {}x{}...'.format(sym.shape(Matrix(onx))[0], sym.shape(Matrix(onx))[1]))
        # =============================================================================================
        # The observability/identifiability matrix is saved in a .txt file
        file = open("results/obs_ident_matrix_{}_{}_Lie_deriv.txt".format(options.modelname, nd), "w")
        file.write('onx = {}'.format(str(onx.tolist())))
        file.close()
        # =============================================================================================
        rational_onx = rationalize_all_numbers(Matrix(onx))
        rango = st.generic_rank(Matrix(rational_onx))
        toc = time() - tic
        print('\n     Rank = {} (calculated in {} seconds)'.format(rango, toc))
        [identifiables, nonidentif, obs_states, unobs_states, obs_inputs, unobs_inputs] = elim_and_recalc(unmeas_x_indices, rango, onx, identifiables, obs_states, obs_inputs)
    #======================================================================================
    # Build the vectors of identifiable / unidentifiable parameters, and of observable / unobservable states and inputs:
    if len(identifiables) != 0:
        p_id = Matrix(identifiables).T
        p_id = np.array(p_id).tolist()[0]
    else:
        p_id = []

    if len(nonidentif) != 0:
        p_un = Matrix(nonidentif).T
        p_un = np.array(p_un).tolist()[0]
    else:
        p_un = []

    if len(obs_states) != 0:
        obs_states = Matrix(obs_states).T
        obs_states = np.array(obs_states).tolist()[0]

    if len(unobs_states) != 0:
        unobs_states = Matrix(unobs_states).T
        unobs_states = np.array(unobs_states).tolist()[0]

    if len(obs_inputs) != 0:
        obs_inputs = Matrix(obs_inputs).T
        obs_inputs = np.array(obs_inputs).tolist()[0]

    if len(unobs_inputs) != 0:
        unobs_inputs = Matrix(unobs_inputs).T
        unobs_inputs = np.array(unobs_inputs).tolist()[0]
    #========================================================================================
    # The observability/identifiability matrix is saved in a .txt file
    file = open("results/obs_ident_matrix_{}_{}_Lie_deriv.txt".format(options.modelname, nd), "w")
    file.write('onx = {}'.format(str(onx.tolist())))
    file.close()

    # The summary of the results is saved in a .txt file
    file = open("results/id_results_{}_{}.txt".format(options.modelname, datetime.today().strftime('%d-%m-%Y')), "w")
    file.write('\n RESULTS SUMMARY:')

    # Report results:
    print('\n ------------------------ ')
    print('     RESULTS SUMMARY:')
    print(' ------------------------ ')
    if len(p_id) == len(model.p) and len(obs_states) == len(model.x) and len(obs_inputs) == len(model.w):
        print('\n >>> The model is Fully Input-State-Parameter Observable (FISPO):')
        file.write('\n >>> The model is Fully Input-State-Parameter Observable (FISPO):')
        if len(model.w) > 0:
            print('\n     All its unknown inputs are observable.')
            file.write('\n     All its unknown inputs are observable.')
        print('\n     All its states are observable.')
        print('\n     All its parameters are locally structurally identifiable.')
        file.write('\n     All its states are observable.\n     All its parameters are locally structurally identifiable.')
    else:
        if len(p_id) == len(model.p):
            print('\n >>> The model is structurally identifiable:')
            print('\n     All its parameters are structurally identifiable.')
            file.write('\n >>> The model is structurally identifiable:\n     All its parameters are structurally identifiable.')
        else:
            if unidflag:
                print('\n >>> The model is structurally unidentifiable.')
                print('\n >>> These parameters are identifiable:\n      {} '.format(p_id))
                print('\n >>> These parameters are unidentifiable:\n      {}'.format(p_un))
                file.write('\n >>> The model is structurally unidentifiable.\n >>> These parameters are identifiable:\n      {}\n >>> These parameters are unidentifiable:\n      {}'.format(p_id,p_un))
            else:
                print('\n >>> These parameters are identifiable:\n      {}'.format(p_id))
                file.write('\n >>> These parameters are identifiable:\n      {}'.format(p_id))

        if len(obs_states) > 0:
            print('\n >>> These states are observable (and their initial conditions, if unknown, are identifiable):\n      {}'.format(obs_states))
            file.write('\n >>> These states are observable (and their initial conditions, if unknown, are identifiable):\n      {}'.format(obs_states))
        if len(unobs_states) > 0:
            print('\n >>> These states are unobservable (and their initial conditions, if unknown, are unidentifiable):\n      {}'.format(unobs_states))
            file.write('\n >>> These states are unobservable (and their initial conditions, if unknown, are unidentifiable):\n      {}'.format(unobs_states))

        if len(meas_x) != 0: # para mostrarlo en una fila, como el resto
            meas_x = Matrix(meas_x).T
            meas_x = np.array(meas_x).tolist()[0]
        else:
            meas_x = []

        if len(meas_x) > 0:
            print('\n >>> These states are directly measured:\n      {}'.format(meas_x))
            file.write('\n >>> These states are directly measured:\n      {}'.format(meas_x))
        if len(obs_inputs) > 0:
            print('\n >>> These unmeasured inputs are observable:\n      {}'.format(obs_inputs))
            file.write('\n >>> These unmeasured inputs are observable:\n      {}'.format(obs_inputs))
        if len(unobs_inputs) > 0:
            print('\n >>> These unmeasured inputs are unobservable:\n      {}'.format(unobs_inputs))
            file.write('\n >>> These unmeasured inputs are unobservable:\n      {}'.format(unobs_inputs))
        if len(model.u) > 0:
            print('\n >>> These inputs are known:\n      {}'.format(model.u))
            file.write('\n >>> These inputs are known:\n      {}'.format(model.u))
    #========================================================================================
    totaltime = time() - tStart # time is measured in seconds
    print('\n Total execution time: {} {} \n\n'.format(totaltime, 'seconds'))
    file.write('\n Total execution time: {} {} \n\n'.format(totaltime, 'seconds'))
    file.close()
