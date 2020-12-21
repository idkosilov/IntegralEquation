from math import *
import matplotlib.pyplot as plt
import sys


def kernel(x, s):
    return s + sin(x)


def test(x):
    return [*map(lambda x: x, x)]


def f(x):
    return [*map(lambda x: x + pi ** 2 * (1 / 2 * sin(x) + 1 / 3 * pi), x)]


def solve(matrix, fn):
    n = len(fn) - 1
    for i in range(n):
        lead, lead_index = abs(matrix[i][i]), i
        for j in range(i + 1, n + 1):
            if abs(matrix[j][i]) > lead:
                lead, lead_index = abs(matrix[i][j]), j
        matrix[i], matrix[lead_index] = matrix[lead_index], matrix[i]
        for j in range(i + 1, n + 1):
            if abs(matrix[j][i]) > sys.float_info.epsilon:
                fn[j] = -fn[j] * matrix[i][i] / matrix[j][i] + fn[i]
                for m in range(n, i-1, -1):
                    matrix[j][m] = -matrix[j][m] * matrix[i][i] / matrix[j][i] + matrix[i][m]
    for i in range(n, -1, -1):
        for j in range(i + 1, n + 1):
            fn[i] -= fn[j] * matrix[i][j]
        fn[i] /= matrix[i][i]
    return fn


def left_rectangle(xn, kernel, h):
    n = len(xn) - 1
    matrix = [[1.0 if i == j else 0.0 for i in range(n+1)] for j in range(n+1)]
    for i in range(n+1):
        for j in range(n):
            matrix[i][j] += h*kernel(xn[i], xn[j])
    return matrix


def simpson(xn, kernel, h):
    n = len(xn)-1
    matrix = [[1.0 if i == j else 0.0 for i in range(n + 1)] for j in range(n + 1)]
    for i in range(n + 1):
        for j in range(1, n, 2):
            matrix[i][j - 1] += kernel(xn[i], xn[j - 1]) * 2 * h / 6
            matrix[i][j] += kernel(xn[i], xn[j]) * 8 * h / 6
            matrix[i][j + 1] += kernel(xn[i], xn[j + 1]) * 2 * h / 6
    return matrix


def show_plots(x, u_lrq, u_simpson, u):
    xu = [x[0] + j * (x[-1] - x[0]) / 10000 for j in range(10001)]
    fig, axs = plt.subplots(1, 2, figsize=(7, 4), constrained_layout=True)
    fig.suptitle('Численные решения интегрального уравнение', fontsize=12)
    axs[0].plot(x, u_lrq, 'o', xu, u(xu), '-')
    axs[0].set_title('Квадратурная формула левых прямоугольников', fontsize=8)
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('U(x)')
    axs[1].plot(x, u_simpson, 'o', xu, u(xu), '-')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('U(x)')
    axs[1].set_title('Квадратурная формула Симпсона', fontsize=8)
    n = len(x) - 1
    fig.savefig(f'integration{n}.png', dpi=800)
    plt.show()


def main():
    a, b, n = 0, pi, 32
    h = (b-a)/n
    print(h)
    xn = [a + i*h for i in range(n+1)]
    U_left_rectangle = solve(left_rectangle(xn, kernel, h), f(xn))
    U_simpson = solve(simpson(xn, kernel, h), f(xn))
    print(U_left_rectangle)
    show_plots(xn, U_left_rectangle, U_simpson, test)


if __name__ == '__main__':
    main()

