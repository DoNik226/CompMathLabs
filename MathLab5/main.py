import numpy as np
from numpy import genfromtxt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import sys
import math
import pandas as pd
from math import factorial as fc


def in_float(s='Введите число', integer=False, check=[False, 0, 0]):
    flag = True
    while flag:
        flag = False
        try:
            if integer:
                val = int(input(s + ': '))
            else:
                val = float(input(s + ': '))
            if check[0] and (val < check[1] or val > check[2]):
                raise ValueError
        except ValueError:
            flag = True
            if check[0]:
                print(f'Попробуйте снова! Введенное число должно принадлежать интервалу [{check[1]}; {check[2]}]\n')
            else:
                print(f'Попробуйте снова!\n')
    return val


def parse():
    flag = True
    while flag:
        path = "input"
        try:
            a = genfromtxt(path, delimiter=',')
            if True in np.isnan(a) or a.shape[0] != 2:
                raise ValueError
            return a
        except ValueError:
            print('В файле должно быть 2 строчки, в каждой одинаковое количество чисел\n')
        except OSError:
            print('Такого файла нет.\n')
        print('Попробуйте снова!\n')


def input_vals():
    n = in_float(s='Введите количество точек', integer=True)
    print()
    a = []
    for i in range(int(n)):
        a.append([in_float('x'), in_float('y')])
        print()
    return np.array(a).transpose()

def f1(x):
    return math.sin(x)+5

def f2(x):
    return x**2+3*x-2

def fun():
    flag=True
    while flag:
        flag=False
        fun=int(input("Выберите функцию:\n" +
            "1. sin(x)+5\n" +
            "2. x^2+3*x-2\n"))
        if fun==1:
            f=f1
        elif fun==2:
            f=f2
        else:
            print("Введен неверный номер, попробуйте еще раз")
            flag=True
    a=int(input("Введите левую границу интервалла. a = "))
    b = int(input("Введите правую границу интервалла. b = "))
    n = int(input("Введите количество точек. n = "))
    a=min(a, b)
    b=max(a, b)
    x=[]
    y=[]
    h=(b-a)/(n-1)
    for i in range(n):
        x.append(a+i*h)
        y.append(f(a+i*h))
    return [x,y]

def to_df(array, f, eps, f_type='f'):
    np_arr = np.concatenate((array, [f], [eps]), axis=0)
    return pd.DataFrame(data=np_arr, index=["X", "Y", f_type, 'eps'])


def newline(p1, p2, color='black'):
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if (p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1] + (p2[1] - p1[1]) / (p2[0] - p1[0]) * (xmax - p1[0])
        ymin = p1[1] + (p2[1] - p1[1]) / (p2[0] - p1[0]) * (xmin - p1[0])

    l = mlines.Line2D([xmin, xmax], [ymin, ymax], color=color)
    ax.add_line(l)
    return l


def check_and_draw(x, y, approximate_function, z,  point):
    fig, ax = plt.subplots()
    xnew = np.linspace(np.min(x), np.max(x), 100)
    ynew = [approximate_function(x, y, i) for i in xnew]
    plt.plot(x, y, 'o', label='Входные точки')
    plt.plot(xnew, ynew, label='Функция аппр.')
    plt.plot(point[0], point[1], '.', markersize=12, label='Аппроксимация')
    plt.title(z)
    ax.legend()
    plt.grid(True)
    plt.show()


def interpolate_lagrange(x, y, x_cur):
    res = 0.0
    for i in range(0, len(x)):
        p = 1.0
        for j in range(0, len(x)):
            if (i != j):
                p *= (x_cur - x[j]) / (x[i] - x[j])
        res += p * y[i]
    return res


def interpolate_newton(x, y, x_cur):
    coef = [[0] * len(x) for i in range(len(x))]
    for i in range(len(x)):
        coef[i][0]=y[i]
    for i in range(1, len(x)):
        for j in range(0, len(x)-i):
            coef[j][i]=(coef[j+1][i-1]-coef[j][i-1])/(x[j+i]-x[j])

    return y[0]+sum(np.prod([x_cur-x[i] for i in range(k)])*coef[0][k] for k in range(1, len(x)))


def gauss_forward_interpolation(x, y, x_cur):
    coef = [[0] * len(x) for i in range(len(x))]
    t=(x_cur-x[(len(x)-1)//2])/((x[-1]-x[0])/(len(x)-1))
    for i in range(len(x)):
        coef[i][0] = y[i]
    for i in range(1, len(x)):
        for j in range(0, len(x) - i):
            coef[j][i] = (coef[j + 1][i - 1] - coef[j][i - 1])

    i = (len(x) - 1) // 2
    res = y[i]
    i-=1
    #print(res)
    p = 1
    f = 1
    j = 1
    for k in range(1, (len(x) - 1) // 2 + 1):
        p *= (t - (k - 1))
        res += p * coef[i][j] / math.factorial(f)
        j += 1
        f += 1
        p *= (t + k)
        res += p * coef[i][j] / math.factorial(f)
        j += 1
        i-=1
        f+=1
    return (res)

def gauss_backward_interpolation(x, y, x_cur):
    coef = [[0] * len(x) for i in range(len(x))]
    t=(x_cur-x[(len(x)-1)//2])/((x[-1]-x[0])/(len(x)-1))
    for i in range(len(x)):
        coef[i][0]=y[i]
    for i in range(1, len(x)):
        for j in range(0, len(x)-i):
            coef[j][i]=(coef[j+1][i-1]-coef[j][i-1])

    i=(len(x)-1)//2
    res=y[i]
    p=1
    f=1
    j=1
    for k in range(1, (len(x)-1)//2+1):
        p*=(t+(k-1))
        res+=p*coef[i][j]/math.factorial(f)
        i-=1
        j+=1
        f+=1
        p *= (t - k)
        res += p * coef[i][j] / math.factorial(f)
        j+= 1
        f+=1
    return(res)

def interpolate_gauss(x, y, x_cur):
    if x_cur > x[(len(x)-1)//2]:
        return gauss_backward_interpolation(x, y, x_cur)
    else:
        return gauss_forward_interpolation(x, y, x_cur)




def run():
    again = True
    while again:
        again = False
        in_type = input('Введите:\n\t* k - если вводить с клавиатуры\n\t* f - если вводить из файла\n\t* fun - если вводить функцию\n')
        if in_type.strip() == 'k':
            data = input_vals()
        elif in_type.strip() == 'f':
            data = parse()
        elif in_type.strip() == 'fun':
            data = fun()
        else:
            print('Введено неверно, попробуйте снова.')
            again = True
    print(f'Итоговый сет данных:\n{data[0]}\n{data[1]}')

    cur_x = in_float('Число, для которого интерполировать значение: ', check=[True, min(data[0]), max(data[0])])


    lagrange_result = interpolate_lagrange(data[0], data[1], cur_x)
    print(f'Lagrange\nОтвет методом Лагранжа: {lagrange_result}')
    check_and_draw(data[0], data[1], interpolate_lagrange, "Лаграндж", [cur_x, lagrange_result])
    print()

    flag=False
    for i in range(2, len(data[0])):
        if round(data[0][i]-data[0][i-1], 3)!=round(data[0][i-1]-data[0][i-2], 3):
            flag=True
            break
    if (flag):
        print("\nТаблица конечных разностей:")
        coef = [[0] * len(data[0]) for i in range(len(data[0]))]
        for i in range(len(data[0])):
            coef[i][0] = data[1][i]
        for i in range(1, len(data[0])):
            for j in range(0, len(data[0]) - i):
                coef[j][i] = (coef[j + 1][i - 1] - coef[j][i - 1]) / (data[0][j + i] - data[0][j])
        print("X:     ", end=' | ')
        for i in data[0]:
            print("{:10.3f}".format(i), end=' | ')
        print()
        print("Y:     ", end=' | ')
        for i in data[1]:
            print("{:10.3f}".format(i), end=' | ')
        print()
        for j in range(1, len(data[0])):
            print("delta", j, end=' | ')
            for i in range(0, len(data[0]) - j):
                print("{:10.3f}".format(coef[i][j]), end=' | ')
            print()
        print()
        newton_result = interpolate_newton(data[0], data[1], cur_x)
        print(f'Newton\nОтвет методом Ньютона: {newton_result}')
        check_and_draw(data[0], data[1], interpolate_newton, "Ньютон", [cur_x, newton_result])
        print()

    if (not flag):
        print("\nТаблица конечных разностей:")
        coef = [[0] * len(data[0]) for i in range(len(data[0]))]
        for i in range(len(data[0])):
            coef[i][0] = data[1][i]
        for i in range(1, len(data[0])):
            for j in range(0, len(data[0]) - i):
                coef[j][i] = (coef[j + 1][i - 1] - coef[j][i - 1])
        print("X:     ", end=' | ')
        for i in data[0]:
            print("{:10.3f}".format(i), end=' | ')
        print()
        print("Y:     ", end=' | ')
        for i in data[1]:
            print("{:10.3f}".format(i), end=' | ')
        print()
        for j in range(1, len(data[0])):
            print("delta", j, end=' | ')
            for i in range(0, len(data[0]) - j):
                print("{:10.3f}".format(coef[i][j]), end=' | ')
            print()
        print()
        if len(data[0])%2==1:
            gauss_result = interpolate_gauss(data[0], data[1], cur_x)
            print(f'Gauss\nОтвет методом Гаусса: {gauss_result}')
            check_and_draw(data[0], data[1], interpolate_gauss,  "Гаусс", [cur_x, gauss_result])
        else:
            print("Чтобы использовать формулу Гаусса, должно быть нечетное количество точек")

run()

