
import numpy as np
import matplotlib.pyplot as plt
import GPy
from time import time
from scipy.interpolate import Rbf

def predict(mu, t, params, dataset, dyas=True):

    if dyas: params = params[:, [1, 2, 3, 6, 7, 8]]

    dmd = DMD(svd_ran=10)
    dmd.fit(dataset)
    dmd.original_time['t0'] = 12000
    dmd.original_time['tend'] = 20000
    dmd.original_time['dt'] = 1
    dmd.dmd_time['tend'] = t

    future_state = dmd.reconstructed_data.real[:, -1] # output at t = 30 s

    kernel = GPy.kern.RBF(input_dim=params.shape[1], variance=1., lengthscale=1.)
    model = GPy.models.GPRegression(params, future_state.reshape(-1, 1), kernel, normalizer=True)
    model.optimize_restarts(num_restarts=10, verbose=False)

    return model.predict(np.atleast_2d(mu))[0][0][0]
    
if __name__ == '__main__':
    n_params_train = 70
    params_top = np.load('parTop.npy')
    params_bottom = np.load('parBot.npy')

    params = np.concatenate((params_top, params_bottom), axis=1)

    params_train = params[:n_params_train, :]
    lift_train = np.load('lift_170par_20000t.npy')[:n_params_train, :]
    print(predict(mu, t, params_train, lift_train))
