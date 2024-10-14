import numpy as np
from matplotlib import pyplot as plt
import math
import sys
import sympy
X, Y = sympy.symbols("x,y")

def f1(x, y):
    return y/x

def f2(x, y):
    return y+(1+x)*y**2

def f3(x, y):
    return x

def p1(x):
    return x

def p2(x):
    return -math.e**x/(x*math.e**x)

def p3(x):
    return x**2/2

def show(x, y, desc=""):
  plt.plot(x,y,'o')
  plt.xlabel("Value of x")
  plt.ylabel("Value of y")
  plt.title(desc)
  plt.show()

# Методом мод. Эйлера
def euler_modified(f, a, b, y0, h):
    n = int(math.ceil((b - a) / h))
    resultx = [a]
    resulty = [y0]
    resultf = []

    for i in range(n):
        xi = float(resultx[i])
        yi = float(resulty[i])
        xi1 = float(xi + h)
        y1 = f(xi, yi)
        y2 = f(xi1, yi + h * y1)
        resultf.append(y2)
        yi1 = float(yi + h / 2 * (y1 + y2))
        resultx.append(xi1)
        resulty.append(yi1)
    return resultx, resulty, resultf

#Метод Рунге-Кутта
def runge_kutt(f, a, b, y0, h):
    n = int(math.ceil((b - a) / h))
    resultx = [a]
    resulty = [y0]
    resultf = []
    for i in range(n):
        xi = float(resultx[i])
        yi = float(resulty[i])
        resultf.append(f(xi, yi))
        xi1 = float(xi + h)
        k1 = h*f(xi, yi)
        k2 = h*f(xi+h/2, yi + k1/2)
        k3 = h*f(xi+h/2, yi + k2/2)
        k4 = h*f(xi+h, yi + k3)
        yi1 = float(yi + (k1+2*k2+2*k3+k4)/6)
        resultx.append(xi1)
        resulty.append(yi1)
    return resultx, resulty, resultf


def runge_rule(f, a, b, y0, h, eps, m, p):
    eps_cur = eps
    y1=0
    i=2
    resx, resy, resf = m(f, a, b, y0, h)
    y=resy[1]
    r = 10000
    while r >= eps_cur:
        prev_resx, prev_resy, prev_resf = resx, resy, resf
        h /= 2
        resx, resy, resf = m(f, a, b, y0, h)
        y1=y
        y=resy[i]
        i*=2
        r = abs(resy[-1] - prev_resy[-1])/(2**p-1)
    return resx, resy, resf


def adams(f, p, a, b, y0, h, eps):
    while True:
        n = int(math.ceil((b - a) / h))
        resultx1, resulty1, resultf1 = runge_kutt(f, a, a+3*h, y0, h)
        resultx = resultx1[:4]
        resulty = resulty1[:4]
        resultf = resultf1[:3]
        for i in range(3, n):
            xi = float(resultx[i])
            yi = float(resulty[i])
            resultf.append(f(xi, yi))
            xi1 = float(xi+h)
            y1 = yi + h/24*(55*resultf[i] -59*resultf[i-1] + 37*resultf[i-2] -9*resultf[i-3])
            resultf.append(f(xi1, y1))
            yi1 = yi + h/24*(9*resultf[i+1] +19*resultf[i] - 5*resultf[i-1] + resultf[i-2])
            resultf[-1] = f(xi1, yi1)
            resultx.append(xi1)
            resulty.append(yi1)
        toch_resulty = []
        e=[]
        for i in range(len(resultx)):
            toch_resulty.append(p(resultx[i]))
            e.append(abs(resulty[i]-toch_resulty[i]))
        if max(e)<=eps:
            break
        h/=2

    return resultx, resulty, resultf, toch_resulty


def run():
    print("Выберите функцию:\n" +
        "1. y' = y/x\n" +
        "2. y' = y + (1+x) * y^2\n " +
        "3. y' = x\n")
    flag = True
    while flag:
        flag=False
        fun = int(input())
        if fun==1:
            f=f1
            p=p1
        elif fun==2:
            f=f2
            p=p2
        elif fun==3:
            f=f3
            p=p3
        else:
            print("Введены неверные данные. Попробуйте снова")
            flag=True
    a = float(input("Введите начало интервала а: "))
    b = float(input("Введите конец интервала b: "))
    h = float(input("Введите шаг h: "))
    y0 = float(input("Введите y0: "))
    eps = float(input("Введите точность eps: "))
    if a > b:
        a, b = b, a
    print("\nМодифицированный метод Эйлера\n")
    euler_resultx, euler_resulty, euler_resultf = runge_rule(f, a, b, y0, h, eps, euler_modified, 1)
    print("i       ", end=' | ')
    print("Xi      ", end=' | ')
    print("Yi      ", end=' | ')
    print("Phi_i   ", end=' | ')
    print()
    for i in range(len(euler_resultx)):
        print("{:8.3f}".format(i), end=' | ')
        print("{:8.3f}".format(euler_resultx[i]), end=' | ')
        print("{:8.3f}".format(euler_resulty[i]), end=' | ')
        if i!=len(euler_resultx)-1:
            print("{:8.3f}".format(euler_resultf[i]), end=' | ')
        print()

    print("\nМетод Рунге\n")
    runge_resultx, runge_resulty, runge_resultf = runge_rule(f, a, b, y0, h, eps, runge_kutt, 4)
    print("i       ", end=' | ')
    print("Xi      ", end=' | ')
    print("Yi      ", end=' | ')
    print("Phi_i   ", end=' | ')
    print()
    for i in range(len(runge_resultx)):
        print("{:8.3f}".format(i), end=' | ')
        print("{:8.3f}".format(runge_resultx[i]), end=' | ')
        print("{:8.3f}".format(runge_resulty[i]), end=' | ')
        if i != len(runge_resultx) - 1:
            print("{:8.3f}".format(runge_resultf[i]), end=' | ')
        print()

    print("\nМетод Адамса\n")
    adams_resultx, adams_resulty, adams_resultf, toch_resulty = adams(f, p, a, b, y0, h, eps)
    print("i       ", end=' | ')
    print("Xi      ", end=' | ')
    print("Yi      ", end=' | ')
    print("Phi_i   ", end=' | ')
    print("Y_toch  ", end=' | ')
    print()
    for i in range(len(adams_resultx)):
        print("{:8.3f}".format(i), end=' | ')
        print("{:8.3f}".format(adams_resultx[i]), end=' | ')
        print("{:8.3f}".format(adams_resulty[i]), end=' | ')
        if i != len(adams_resultx) - 1:
            print("{:8.3f}".format(adams_resultf[i]), end=' | ')
        else:
            print("        ", end=' | ')
        print("{:8.3f}".format(toch_resulty[i]), end=' | ')
        print()


    h=(b-a)/1000
    X=[]
    Y=[]
    for i in range(1001):
        X.append(a+i*h)
        Y.append(p(X[i]))


    plt.plot(euler_resultx, euler_resulty, label="Мод. метод Эйлера")
    plt.plot(X, Y, label="Точное значение функции")
    plt.grid(True)
    plt.legend()
    plt.show()
    plt.plot(runge_resultx, runge_resulty, label="Метод Рунге-Кутта")
    plt.plot(X, Y, label="Точное значение функции")
    plt.grid(True)
    plt.legend()
    plt.show()
    plt.plot(adams_resultx, adams_resulty, label="Метод Адамса")
    plt.plot(X, Y, label="Точное значение функции")
    plt.grid(True)
    plt.legend()
    plt.show()

run()