# ФУНКЦИИ

def f1(x, x1, z1, z_1, K1):
    return z1 + z_1 * (x - x1) - K1 * (x1 - x) ** 2 / 2

def f2(x, x2, z2, z_2, K1):
    return z2 - z_2 * (x2 - x) - K1 * (x2 - x) ** 2 / 2

def f12(x, C0, C1, K1):
    return C0 + C1*x + K1*x**2/2

def c0(x2, z2, z_2, X_2, K1):
    return z2 - z_2*x2 - K1*x2**2/2 + K1*X_2**2

def c1(x2, z_2, X_2, K1):
    return z_2 - 2*K1*X_2 + K1*x2

# ТОЧКИ

def x_hat(x1, z1, z_1, x2, z2, z_2, K1):
    # print('знаменатель =', (z_1 - z_2 - K1 * (x2 - x1))) #
    # print('K1 =', K1)
    return (z2 - z1 - z_2 * x2 + z_1 * x1 - K1 * (x2 ** 2 - x1 ** 2) / 2) / (z_1 - z_2 - K1 * (x2 - x1))

def x_1(X_hat, x1, x2, z_1, z_2, K1):    # Это x_i'
    return X_hat - (x2-x1)/4 - (z_2-z_1) / (4*K1)

def x_2(X_hat, x1, x2, z_1, z_2, K1):    # Это x_i''
    return X_hat + (x2-x1)/4 + (z_2-z_1) / (4*K1)

def x_2hat(x2, X_2, z_2, K1):
    return 2*X_2 - x2 - z_2/K1

def z_2hat(C0, X_2hat, K1):
    return C0 - K1*X_2hat**2/2

# Условие dF/dx

def concider_z_2hat(c1, X_1, X_2, K1):
    return c1 + K1*X_1 < 0 and c1 + K1*X_2 > 0


def count_all(x1, z1, z_1, x2, z2, z_2, K1):    # посчитаем все сразу
    # x_1, x_2, x_2hat, z_2hat, cons_z_2hat
    X_hat = x_hat(x1, z1, z_1, x2, z2, z_2, K1)
    X_1 = x_1(X_hat, x1, x2, z_1, z_2, K1)
    X_2 = x_2(X_hat, x1, x2, z_1, z_2, K1)
    X_2hat = x_2hat(x2, X_2, z_2, K1)

    C0 = c0(x2, z2, z_2, X_2, K1)
    C1 = c1(x2, z_2, X_2, K1)
    Z_2hat = z_2hat(C0, X_2hat, K1)

    cons_z_2hat = concider_z_2hat(C1, X_1, X_2, K1)
    return X_1, X_2, X_2hat, Z_2hat, cons_z_2hat

def count_curve_params(x1, z1, z_1, x2, z2, z_2, K1):
    X_hat = x_hat(x1, z1, z_1, x2, z2, z_2, K1)
    # print('X_hat =', X_hat)
    X_1 = x_1(X_hat, x1, x2, z_1, z_2, K1)
    X_2 = x_2(X_hat, x1, x2, z_1, z_2, K1)
    X_2hat = x_2hat(x2, X_2, z_2, K1)  #

    C0 = c0(x2, z2, z_2, X_2, K1)
    C1 = c1(x2, z_2, X_2, K1)
    Z_2hat = z_2hat(C0, X_2hat, K1)  #
    # Z_2hat = f12(X_2hat, C0, C1, K1)
    return X_1, X_2, C0, C1, X_2hat, Z_2hat
