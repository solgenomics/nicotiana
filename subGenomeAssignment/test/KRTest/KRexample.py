import matplotlib.pyplot as plt 
import numpy as np 
from sklearn.neighbors.kde import KernelDensity
from scipy.stats import norm

if __name__ == "__main__":
    np.random.seed(1)

    x = np.concatenate((np.random.normal(0, 1, int(0.3*100)), np.random.normal(5, 1, int(0.7*100))))[:, np.newaxis]
    # [:,np.newaxis] converts the 1D array into a column vector
    plot_x = np.linspace(-5, 10, 1000)[:, np.newaxis]
    true_dens = 0.3*norm(0, 1).pdf(plot_x) + 0.7*norm(5, 1).pdf(plot_x)

    log_dens = KernelDensity(bandwidth=1).fit(x).score_samples(plot_x)

    plt.figure(),
    plt.fill(plot_x, true_dens, fc='#AAAAFF', label='true_density')
    plt.plot(plot_x, np.exp(log_dens), 'r', label='estimated_density')
    for _ in range(x.shape[0]):
        plt.plot(x[:, 0], np.zeros(x.shape[0])-0.01, 'g*') 
    plt.legend()

    plt.show()
