import numpy as np
import matplotlib.pyplot as plt
from funcs import *

def deriv(f, x, eps=1e-3):
    # вычисление численных производных
    return (f(x+eps) - f(x-eps))/(2*eps)

def m_1(q):
    # бичарский вариант
    n = len(q)
    m = []
    for i in range(1, n):
        x1, _, z_1 = q[i - 1]
        x2, _, z_2 = q[i]
        m.append(np.abs(z_2-z_1)/(x2-x1))
    return m

def m_3(q):
    # подсчет характеристик отрезков для оценки константы Липшица
    def tau(zz, delta, t):
        return np.abs(2 * (zz + delta * t) / (1 + 2 * t ** 2))
    n = len(q)
    m = [None] * (n-1)
    for i in range(1, n):
        x1, z1, z_1 = q[i - 1]
        x2, z2, z_2 = q[i]
        h, delta = (x2-x1)/2, z_2-z_1
        zz = (z2 - z1)/h - z_2 - z_1
        if delta != 0:
            t1 = -zz/delta + ((zz/delta)**2 + 0.5)**0.5
            t2 = -zz/delta - ((zz/delta)**2 + 0.5)**0.5
        else:
            t1 = t2 = 0
        m[i-1] = max(tau(zz, delta, -1), tau(zz, delta, 1))
        if -1 < t1 < 1 and -1 < t2 < 1:
            m[i-1] = max(m[i-1], tau(zz, delta, t1), tau(zz, delta, t2))
    return m

def m_3_new(q, r=8):
    # кусок кода Славы, он просто бог проги
    qq = np.array(q)
    x = qq[:, 0]
    dx = x[1:] - x[:-1]
    y = qq[:, 1]
    y_dot = qq[:, 2]

    def calculate_tau(i, pt):    #
        dx_f = pt - x[i]
        dx_b = pt - x[i + 1]
        return 2 * abs(y[i + 1] - y[i] + y_dot[i + 1] * dx_b
                       - y_dot[i] * dx_f) / (dx_b**2 + dx_f**2)

    def calculate_nu(i):
        if abs(y_dot[i] - y_dot[i + 1]) < 1e-6:
            return calculate_tau(i, (x[i] + x[i + 1]) / 2)
        d = 0.5 * np.sqrt(
            (2 * (y[i] - y[i + 1]) + (y_dot[i + 1] + y_dot[i]) * dx[i])**2
            + (y_dot[i + 1] - y_dot[i])**2 * dx[i]**2
            )
        d_minus = x[i] + ((y[i] - y[i + 1] + y_dot[i + 1] * dx[i] - d)
                          / (y_dot[i + 1] - y_dot[i]))
        if x[i] <= d_minus <= x[i + 1]:
            tau_minus = calculate_tau(i, d_minus)
        else:
            tau_minus = -np.inf
        d_plus = x[i] + ((y[i] - y[i + 1] + y_dot[i + 1] * dx[i] + d)
                         / (y_dot[i + 1] - y_dot[i]))
        if x[i] <= d_plus <= x[i + 1]:
            tau_plus = calculate_tau(i, d_plus)
        else:
            tau_plus = -np.inf
        return max(
            calculate_tau(i, x[i]),
            calculate_tau(i, x[i + 1]),
            tau_minus,
            tau_plus
        )

    return [calculate_nu(i) for i in range(len(q)-1)]


def K_Lipschitz(q, mode1, xi=1e-9, r=1.1):
    # оценка константы Липшица для разных случаев
    n = len(q)
    # m = m_1(q)
    # m = m_3(q)    # не работает
    m = m_3_new(q, r=r)
    G = np.max(m)
    if mode1 == 'Global':
        if G == 0:
            K = [1] * (n-1)
        else:
            K = [r*G] * (n-1)
    elif mode1 == 'Local':
        mm = m.copy()  # m_i в формуле для K_i
        d = []
        for i in range(1, n):
            x1, *_ = q[i - 1]
            x2, *_ = q[i]
            d.append(x2 - x1)
            if i > 1:
                mm[i - 1] = max(m[i - 1], m[i - 2])
            if i < n - 1:
                mm[i - 1] = max(m[i - 1], m[i])
        d_max = np.max(d)
        K = []
        # print(mm)
        for i in range(n-1):
            K.append(r * max([mm[i], G*d[i]/d_max, xi]))
    return K

def merge(q, *point):
    res = []
    for i in range(len(q)):
        if q[i][0] < point[0]:
            res.append(q[i])
        else:
            res.append(point)
            res.extend(q[i:])
            break
    return res

def which_is_min_in(arr, *which):
    # функция с говорящим названием
    for x in which:
        if x == arr:
            return x
    raise KeyError

def find_new_R(q, K, mode2):
    # подсчет характеристики отрезков для определения следующего измерения
    # кроме самого R возвращает еще и список, в котором для каждого отрезка лежит точка,
    # в которой достигается этот минимум
    R = []
    which_is_min = []
    n = len(q)
    if mode2 == 'notSmooth':    # which_is_min is not used
        for i in range(1, n):
            x1, z1, z_1 = q[i - 1]
            x2, z2, z_2 = q[i]
            K1 = K[i-1]
            x = x_hat(x1, z1, z_1, x2, z2, z_2, K1)
            R.append(-min([z1, z2, f2(x, x2, z2, z_2, K1)]))

    if mode2 == 'Smooth':
        for i in range(1, n):
            x1, z1, z_1 = q[i - 1]
            x2, z2, z_2 = q[i]
            K1 = K[i-1]
            X_1, X_2, X_2hat, Z_2hat, cons_z_2hat = count_all(*q[i - 1], *q[i], K1)
            arr_z = [z1, z2]
            arr_x = [X_1, X_2]
            # print(X_2)
            if cons_z_2hat:
                arr_z.append(Z_2hat)
                arr_x.append(X_2hat)
            R.append(-min(arr_z))
            # print('z:', arr_z, 'x:', arr_x)
            which_is_min.append(arr_x[np.argmin(arr_z)])

    return R, which_is_min

def plot_minoranas(q, K, mode2, step=0.001):
    n = len(q)
    x, y = [], []
    curve_x, curve_y = [], []
    if mode2 == 'notSmooth':
        for i in range(1, n):
            x1, z1, z_1 = q[i - 1]
            x2, z2, z_2 = q[i]
            curve_x.extend([x1, x2])
            curve_y.extend([z1, z2])
            K1 = K[i-1]
            xx = np.arange(x1, x2, step)
            yy = np.maximum(f1(xx, x1, z1, z_1, K1), f2(xx, x2, z2, z_2, K1))
            x.extend(xx)
            y.extend(yy)

    if mode2 == 'Smooth':
        for i in range(1, n):
            x1, z1, z_1 = q[i - 1]
            x2, z2, z_2 = q[i]
            K1 = K[i-1]
            curve_x.extend([x1, x2])
            curve_y.extend([z1, z2])
            X_1, X_2, C0, C1, X_2hat, Z_2hat = count_curve_params(*q[i - 1], *q[i], K1)

        # кривая на промежутке x_{i-1},..., x_i'
            xx = np.arange(x1, X_1, step)
            yy = f1(xx, x1, z1, z_1, K1)
            x.extend(xx)
            y.extend(yy)

        # кривая на промежутке x_i',..., x_{i+1}'
            xx = np.arange(X_1, X_2, step)
            yy = f12(xx, C0, C1, K1)
            x.extend(xx)
            y.extend(yy)
        # кривая на промежутке x_{i+1},..., x_i'
            xx = np.arange(X_2, x2, step)
            yy = f2(xx, x2, z2, z_2, K1)
            x.extend(xx)
            y.extend(yy)
    x = np.array(x)
    y = np.array(y)
    return x, y, curve_x, curve_y
