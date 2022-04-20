"""
    Function to check the observability / identifiability of each variable
"""

def elim_and_recalc(unmeas_xred_indices, rangoinicial, numonx, p, x, unidflag, w1vector, *args):
    ################################################################################
    # Imports:
    import options
    import sympy as sym
    from sympy import Matrix
    import symbtools as st
    from functions.rationalize import rationalize_all_numbers
    ################################################################################
    numonx = rationalize_all_numbers(Matrix(numonx))
    # Depending on the number of arguments you pass to the function, there are two cases:

    # called when there is no 'w'
    if len(args) == 0:
        pred = p
        xred = x
        wred = w1vector
        identifiables = []
        obs_states = []
        obs_inputs = []
        q = len(pred)
        n = len(xred)
        nw = len(wred)

    # called when there are 'w'
    if len(args) == 3:
        pred = p
        xred = x
        wred = w1vector
        identifiables = args[0]
        obs_states = args[1]
        obs_inputs = args[2]
        q = len(pred)
        n = len(xred)
        nw = len(wred)

    r = sym.shape(Matrix(numonx))[1] # before: q+n+nw; but with unknown inputs there may also be derivatives
    new_ident_pars = identifiables
    new_nonid_pars = []
    new_obs_states = obs_states
    new_unobs_states = []
    new_obs_in = obs_inputs
    new_unobs_in = []

    #========================================================================
    # ELIMINATE A PARAMETER:
    #========================================================================
    # At each iteration we remove a different column (= parameter) from onx:
    for ind in range(q): # for each parameter of p...
        if q <= 1: # check if the parameter has already been marked as identifiable
            isidentifiable = pred[ind] in identifiables
        else:
            isidentifiable = any(pred[ind] in arr for arr in identifiables)
        if isidentifiable:
            print('\n Parameter {} has already been classified as identifiable.'.format(pred[ind]))
        else:
            indices = []
            for i in range(r):
                indices.append(i)
            indices.pop(n + ind)
            column_del_numonx = Matrix(numonx).col(indices) # one column is removed
            num_rank = st.generic_rank(Matrix(column_del_numonx)) # the range is calculated without that column
            if num_rank == rangoinicial:
                if unidflag == 1:
                    print('\n    => Parameter {} is structurally unidentifiable'.format(pred[ind]))
                    new_nonid_pars.append(pred[ind])
                else:
                    print('\n    => We cannot decide about parameter {} at the moment'.format(pred[ind]))
            else:
                print('\n    => Parameter {} is structurally identifiable'.format(pred[ind]))
                new_ident_pars.append(pred[ind])

    # ========================================================================
    # ELIMINATE A STATE:
    # ========================================================================
    # At each iteration we try removing a different state from 'xred':
    if options.checkObser == 1:
        for ind in range(len(unmeas_xred_indices)): # for each unmeasured state
            original_index = unmeas_xred_indices[ind]
            if len(obs_states) <= 1:
                isobservable = xred[original_index] in obs_states
            else:
                isobservable = any(xred[original_index] in arr for arr in obs_states)
            if isobservable:
                print('\n State %s has already been classified as observable.'.format(xred[ind]))
            else:
                indices = []
                for i in range(r):
                    indices.append(i)
                indices.pop(original_index) # remove the column that we want to check
                column_del_numonx = Matrix(numonx).col(indices)
                num_rank = st.generic_rank(Matrix(column_del_numonx))
                if num_rank == rangoinicial:
                    if unidflag == 1:
                        print('\n    => State {} is unobservable'.format(xred[original_index]))
                        new_unobs_states.append(xred[original_index])
                    else: # if this function was called because the necessary number of derivatives was not calculated...
                        print('\n    => We cannot decide about state {} at the moment'.format(xred[original_index]))
                else:
                    print('\n    => State {} is observable'.format(xred[original_index]))
                    new_obs_states.append(xred[original_index])

    # ========================================================================
    # ELIMINATE AN UNKNOWN INPUT:
    # ========================================================================
    # At each iteration we try removing a different column from onx:
    for ind in range(nw): # for each unknown input...
        if len(obs_inputs) <= 1:  # check if the unknown input has already been marked as observable
            isobservable = wred[ind] in obs_inputs
        else:
            isobservable = any(wred[ind] in arr for arr in obs_inputs)
        if isobservable:
            print('\n Input %s has already been classified as observable.'.format(wred[ind]))
        else:
            indices = []
            for i in range(r):
                indices.append(i)
            indices.pop(n+q+ind) # remove the column that we want to check
            column_del_numonx = Matrix(numonx).col(indices)
            num_rank = st.generic_rank(Matrix(column_del_numonx))
            if num_rank == rangoinicial:
                if unidflag == 1:
                    print('\n    => Input {} is unobservable'.format(wred[ind]))
                    new_unobs_in.append(wred[ind])
                else:
                    print('\n    => We cannot decide about input {} at the moment'.format(wred[ind]))
            else:
                print('\n    => Input {} is observable'.format(wred[ind]))
                new_obs_in.append(wred[ind])
    return new_ident_pars, new_nonid_pars, new_obs_states, new_unobs_states, new_obs_in, new_unobs_in
