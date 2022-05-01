import numpy as np
import matplotlib.pyplot as plt
from tools import *
from funcs import *


func = [
    lambda x: 1 / 6 * x ** 6 - 52 / 25 * x ** 5 + 39 / 80 * x ** 4 + 71 / 10 * x ** 3 - 79 / 20 * x ** 2 - x + 1 / 10,
    lambda x: np.sin(x) + np.sin(10 / 3 * x),
    lambda x: -sum([k * np.sin((k + 1) * x + k) for k in range(1, 6)]),
    lambda x: -(16 * x ** 2 - 24 * x + 5) * np.exp(-x),
    lambda x: (3 * x - 1.4) * np.sin(18 * x),
    lambda x: -(x + np.sin(x)) * np.exp(-x ** 2),
    lambda x: np.sin(x) + np.sin(10 / 3 * x) + np.log(x) - 0.84 * x + 3,
    lambda x: -sum([k * np.cos((k + 1) * x + k) for k in range(1, 6)]),
    lambda x: np.sin(x) + np.sin(2 / 3 * x),
    lambda x: -x * np.sin(x)]
omega = [(-1.5, 11),
         (2.7, 7.5),
         (-10, 10),
         (1.9, 3.9),
         (0, 1.2),
         (-10, 10),
         (2.7, 7.5),
         (-10, 10),
         (3.1, 20.4),
         (0, 10)]


def search_min(f, a, b, r=8, N=150, mode1='Global', mode2='notSmooth',  xi=1e-9):
    eps = (b-a) * 1e-4
    q = [(a, f(a), deriv(f, a)), (b, f(b), deriv(f, b))]    # точки на кривой
    k = 0
    for i in range(N):
        # посчитаем константу Липшица
        K = K_Lipschitz(q, mode1, r=r)
        R, which_is_min_in_R = find_new_R(q, K, mode2)
        t = np.argmax(R) + 1
        if q[t][0] - q[t-1][0] < eps:
            break
        if mode2 == 'notSmooth':
            x = x_hat(*q[t-1], *q[t], K[t-1])    # точка очередного испытания
        elif mode2 == 'Smooth':
            x = which_is_min_in_R[t-1]
        # новое измерение
        z = f(x)
        z_ = deriv(f, x)
        qq = merge(q, x, z, z_)
        k += 1
        q = qq.copy()
    return q, k, *plot_minoranas(q, K, mode2)


def main():

    N = 200
    i = 4
    mode1 = 'Local' # Global, Local
    mode2 = 'notSmooth' # notSmooth, Smooth
    r = 1.1   # 8

    print('i =', i)
    x = np.arange(*omega[i], 0.005)
    y = func[i](x)
    x_real, y_real = x[np.argmin(y)], min(y)
    print(f'Настоящий минимум: {x_real}, {y_real}', '\n')
    q, k, xx, yy, node_x, node_y = search_min(func[i], *omega[i], N=N, r=r, mode1=mode1, mode2=mode2)
    x_min, y_min, z_ = min(q, key=lambda x: x[1])
    print(f'Количество шагов: {k}')

    print(f'Найденный минимум: {x_min}, {y_min}')

    plt.plot(x, y)
    plt.plot(xx, yy)    # миноранты
    plt.scatter(node_x, node_y)    # точки измерений на кривой
    plt.show()

if __name__ == '__main__':
    main()
