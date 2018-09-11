# Copyright (c) 2013, SciPy Developers, Robert McGibbon
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

import numpy as np
from clear_para_file import *

def fmin_bfgs(f,fprime, x0, gtol=1e-5, callback=None, maxiter=None):
    """Minimize a function, via only information about its gradient, using BFGS
    
    The difference between this and the "standard" BFGS algorithm is that the
    line search component uses a weaker criterion, because it can't check
    for sure that the function value actually decreased.
    
    Parameters
    ----------
    fprime : callable f(x, *args)
        gradient of the objective function to be minimized
    x0 : ndarray
        Initial guess
    args : tuple, optional
        Extra arguments to be passed to `fprime`
    gtol : float, optional
        gradient norm must be less than `gtol` before succesful termination
    callback : callable, optional
        An optional user-supplied function to call after each iteration.
        Called as `callback(xk)`, where `xk` is the current parameter vector.
    maxiter : int, optional
        Maximum number of iterations to perform.
    
    Returns
    -------
    xopt : ndarray
        Parameters which minimize `f`, and which are a root of the gradient,
        `fprime`
    gopt : ndrarray
        value of the gradient at `xopt`, which should be near zero
    Hopt : ndarray
        final estimate of the hessian matrix at `xopt`
    n_grad_calls : int
        number of gradient calls made
    """
    
    x0 = np.asarray(x0).flatten()
    if maxiter is None:
            maxiter = len(x0)*200
    
    gfk, old_drag = fprime(x0,1)  # initial gradient
    n_grad_calls = 1  # number of calls to fprime()
    
    k = 0  # iteration counter
    N = len(x0)  # degreees of freedom
    I = np.eye(N, dtype=int)
    
    Hk = I  # initial guess of the Hessian
    xk = x0
    sk = [2*gtol]
    
    gnorm = np.linalg.norm(gfk)
    iterations = 0
    alpha_guess = 0.02 # Initial guess for the step size
    ls_grad_calls = 1
    while (gnorm > gtol) and (k < maxiter) and ls_grad_calls < 20: # 20 is number of max iters
        iterations += 1
	# search direction
        pk = -np.dot(Hk, gfk)
        
        alpha_k, gfkp1, ls_grad_calls, old_drag = _line_search(f,fprime, xk, gfk, pk, alpha_guess, old_drag)
	alpha_guess = alpha_k * 1.5 # Update alpha guess size as gradient decreases
        n_grad_calls += ls_grad_calls
        
        # advance in the direction of the step
        xkp1 = xk + alpha_k * pk
        sk = xkp1 - xk
        xk = xkp1

        if gfkp1 is None:
            gfkp1 = fprime(xkp1,1)
            n_grad_calls += 1
        
        yk = gfkp1 - gfk
        gfk = gfkp1
        
        if callback is not None:
            callback(xk)
        
        k += 1
        gnorm = np.linalg.norm(gfk)
        if gnorm < gtol:
            break
            
        try:  #this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (np.dot(yk, sk))
        except ZeroDivisionError:
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"
        if np.isinf(rhok):  #this is patch for numpy
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"

        # main bfgs update here. this is copied straight from
        # scipy.optimize._minimize_bfgs
        A1 = I - sk[:, np.newaxis] * yk[np.newaxis, :] * rhok
        A2 = I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok
        Hk = np.dot(A1, np.dot(Hk, A2)) + rhok * sk[:, np.newaxis] \
                * sk[np.newaxis, :]
                
    
    if k >= maxiter:
        print "Warning: %d iterations exceeded" % maxiter
        print "         Current gnorm: %f" % gnorm
        print "         grad calls: %d" % n_grad_calls
        print "         iterations: %d" % k
        
    
    elif gnorm < gtol:
        print "Optimization terminated successfully."
        print "         Current gnorm: %f" % gnorm
        print "         grad calls: %d" % n_grad_calls
        print "         iterations: %d" % k
    
    xopt = xk
    gopt = gfk
    Hopt = Hk
    n_grad_calls = n_grad_calls
    iterations = iterations    
    return xopt, gopt, Hopt, n_grad_calls, iterations


def test_fmin_bfgs_1():
    import scipy.optimize
    #def callback(xk):
    #    pass
    print fmin_bfgs(scipy.optimize.rosen_der, [-50, 20]).xopt

    #scipy.optimize.fmin_bfgs(scipy.optimize.rosen, [10, 20], scipy.optimize.rosen_der)


def _line_search(f,fprime, xk, gk, pk, alpha_guess, old_drag, curvature_condition=0.9,
    update_rate=0.5, maxiters=20):
    
    alpha = alpha_guess
    drag_old = old_drag
    with open('drag.txt','a') as myfile:
	myfile.write('Line search')
	myfile.write('Old drag ')
	myfile.write(str(drag_old))
	myfile.write('\n')

    with open('lcst_iterations.txt','a') as myfile:
	myfile.write('Line search')

    for i in range(maxiters):
	clear_para_file()
    	drag = f(xk + alpha * pk)
    	if drag == 'error' or drag > drag_old:
		alpha *= update_rate	
        else:
		break
    
    new_drag = drag
    gk,not_used = fprime(xk + alpha * pk,0)
    with open('alpha.txt','a') as myfile:
                myfile.write('Line search chosen alpha')
                myfile.write('\n')
                myfile.write(str(alpha))

    # i+1 is the final number of calls to f()
    return alpha, gk, i+1, new_drag

if __name__ == '__main__':
    test_fmin_bfgs_1()
