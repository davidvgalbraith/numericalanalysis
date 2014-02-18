import matplotlib.pyplot as plt
import numpy as np


y1 =  [2.0000e+00,   2.9995e+00,   6.6663e-01,   2.1565e-01,   7.8729e-02,   2.9744e-02,   1.1334e-02,   4.3271e-03,   1.6526e-03,   6.3124e-04]

y2 = [   2.0000e+00,   1.0000e+00,   5.5556e-01,   3.3333e-01,   2.0988e-01,   1.3992e-01,   9.1907e-02,   6.0357e-02,   3.9781e-02,   2.6368e-02,   1.7562e-02,   1.1674e-02,   7.7600e-03,   5.1614e-03,   3.4382e-03,   2.2912e-03,   1.5263e-03,   1.0168e-03, 6.7746e-04]

y3 = [   2.0000e+00,   8.3333e-01,   3.6111e-01,   1.6204e-01,   7.4846e-02,   3.5365e-02,    1.6997e-02,   8.2697e-03,   4.0587e-03,   2.0124e-03,   1.0020e-03,   4.9863e-04]
x1 = np.arange(0, len(y1), 1)
x2 = np.arange(0, len(y2), 1)
x3 = np.arange(0, len(y3), 1)

plt.plot(x1, y1, label="Conjugate Gradient")
plt.plot(x2, y2, label="Jacobi Method")
plt.plot(x3, y3, label="Gauss-Seidel Method")
plt.yscale("log")
#plt.axis([.001, .01, .99, 1.01])
plt.title("Convergence Rates for Various Algorithms")
plt.xlabel("Iteration")
plt.ylabel("Error (Log scale)")
plt.legend()
plt.show()
