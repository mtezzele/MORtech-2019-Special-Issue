
import numpy as np
import matplotlib.pyplot as plt
from athena import ActiveSubspaces, Normalizer

def plot_as_at_given_time(params, lift, time):
    output = lift[:, time]

    ss = ActiveSubspaces()
    ss.compute(inputs=params, outputs=output, method='local', nboot=1000)
    ss.partition(1)

    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    axs[0].set(xlim=(-2., 2.), ylim=(-0.1, 0.5))
    axs[1].set(xlim=(0, 11), ylim=(-0.82, 0.82))
    
    axs[0].scatter(params.dot(-ss.W1), output, c='blue', s=40, alpha=0.9, edgecolors='k', label='True lift')
    axs[0].set_xlabel('Active variable true ' + r'$W_1^T \mathbf{\mu}}$', fontsize=18)
    axs[0].set_ylabel(r'$f \, (\mathbf{\mu}, t=$' + str(time/1000.0) + r'$)$', fontsize=18)
    
    axs[1].set_xlabel('Eigenvector component', fontsize=18)
    axs[1].set_xticks(range(1, 11))
    axs[1].set_xticklabels([r'$c_1$', r'$c_2$', r'$c_3$', r'$c_4$', r'$c_5$', r'$d_1$', r'$d_2$', r'$d_3$', r'$d_4$', r'$d_5$'], fontsize=14)
    axs[1].set_ylabel(r'$W_1$' + ' magnitude', fontsize=18)
    axs[1].axhline(linewidth=1, color='black')
    axs[1].scatter(range(1, 11), -ss.W1, c='blue', s=60, alpha=0.9, edgecolors='k')
    
    axs[0].grid(linestyle='dotted')
    axs[1].grid(linestyle='dotted')

    axs[0].legend()
    plt.show()
    
if __name__ == '__main__':
    n_params_train = 70
    params_top = np.load('parameters/parTop.npy')
    params_bottom = np.load('parameters/parBot.npy')
    params = np.concatenate((params_top, params_bottom), axis=1)

    x_low = np.array([0, 0, 0, 0, 0, -0.03, -0.03, -0.03, -0.03, -0.03])
    x_up = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0, 0, 0, 0, 0])
    normalizer = Normalizer(x_low, x_up)
    params = normalizer.normalize(params)
    params_train = params[:n_params_train, :]
    
    lift_train = np.load('outputs/lift_170par_20000t.npy')[:n_params_train, :]
    
    desired_time = 12000
    plot_as_at_given_time(params_train, lift_train, time=desired_time)

    