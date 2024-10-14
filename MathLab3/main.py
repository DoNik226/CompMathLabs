import math
import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
from scipy.integrate import quad
n=0
def f1(x):
    return x**2+x**3-x**4+x
def f2(x):
    return 3*x**3-2*x**2-7*x-8
def f3(x):
    return np.cos(x)*np.sin(x)

def left_rectangle_method(a, b, eps, f):
    global n
    s=0
    n=4
    h = (b - a) / n
    for i in range(n):
        x = a + i * h
        s += f(x)
    s *= h
    n*=2
    while(True):
        sum=0
        h=(b-a)/n
        for i in range(n):
            x=a+i*h
            sum+=f(x)
        sum*=h
        if (abs(s-sum)<eps):
            return sum
        else:
            n*=2
            s=sum
def right_rectangle_method(a, b, eps, f):
    global n
    s=0
    n=4
    h = (b - a) / n
    for i in range(1, n + 1):
        x = a + i * h
        s += f(x)
    s *= h
    n*=2
    while(True):
        sum=0
        h=(b-a)/n
        for i in range(1, n+1):
            x=a+i*h
            sum+=f(x)
        sum*=h
        if (abs(s-sum)<eps):
            print(s)
            return sum
        else:
            n*=2
            s=sum
def medium_rectangle_method(a, b, eps, f):
    global n
    s=0
    n=4
    h = (b - a) / n
    for i in range(n):
        x = a + i * h + h / 2
        s += f(x)
    s*= h
    n*=2
    while(True):
        sum=0
        h=(b-a)/n
        for i in range(n):
            x=a+i*h+h/2
            sum+=f(x)
        sum*=h
        if (abs(s-sum)/3<eps):
            return sum
        else:
            n*=2
            s=sum
def trapezoid_method(a, b, eps, f):
    global n
    s=0
    n=4
    h = (b - a) / n
    for i in range(n):
        x1 = a + i * h
        x2 = a + (i + 1) * h
        s += f(x2) + f(x1)
    s *= h / 2
    n*=2
    while(True):
        sum=0
        h=(b-a)/n
        for i in range(n):
            x1=a+i*h
            x2 = a + (i+1) * h
            sum+=f(x2)+f(x1)
        sum*=h/2
        if (abs(s-sum)/3<eps):
            return sum
        else:
            n*=2
            s=sum
def simpson_method(a, b, eps, f):
    global n
    s=0
    n=4
    h = (b - a) / n
    for i in range(1, n):
        x = a + i * h
        if (i % 2 == 1):
            s += 4 * f(x)
        else:
            s += 2 * f(x)
    s += f(a) + f(b)
    s *= h / 3
    n*=2
    while(True):
        sum=0
        h=(b-a)/n
        for i in range(1, n):
            x=a+i*h
            if (i%2==1):
                sum+=4*f(x)
            else:
                sum+=2*f(x)
        sum+=f(a)+f(b)
        sum*=h/3
        if (abs(s-sum)/15<eps):
            return sum
        else:
            n*=2
            s=sum


print("\nВведите номер функции для интегрирования:\n"+
        "1. y = x^2+x^3-x^4+x\n"+
        "2. y = 3x^3-2x^2-7x-8\n"+
        "3. y = cos(x)*sin(x)")
while(True):
    fun=float(input())
    if (fun == 1):
        f = f1
        break
    elif (fun == 2):
        f = f2
        break
    elif (fun == 3):
        f = f3
        break
    else:
        print("Введен неверный номер функции.\n" + "Попробуй еще раз")

print("Введите левый и правый пределы интегрирования")
while (True):
    a = float(input('a = '))
    b = float(input('b = '))
    if a > b:
        print("Введен неверный предел интегрирования. Попробуйте еще раз")
        continue
    break

print("Введите точность:")
while (True):
    eps = float(input('eps = '))
    if (eps > 0):
        break
    print("Введенное значение точности должно быть < 0. Попробуйте еще раз:")


print("\nВведите метод для интегрирования:\n" +
      "1. Метод прямоугольников\n" +
      "2. Метод трапеций\n" +
      "3. Метод Симпсона\n"+
      "4. Все")
while (True):
    method = float(input())
    if (method == 1 or method == 2 or method == 3 or method ==4):
        break
    print("Введен неверный номер метода. Попробуйте еще раз:")

if method == 1:
    print("\nВведите модификацию для метода прямоугольников:\n" +
          "1. Метод левых прямоугольников\n" +
          "2. Метод правых прямоугольников\n" +
          "3. Метод средних прямоугольников\n" +
          "4. Все 3 модификации")
    while (True):
        modification = float(input())
        if (modification == 1 or modification == 2 or modification == 3 or modification == 4):
            break
        print("Введен неверный номер модификации. Попробуйте еще раз:")
    if modification == 1:
        print("Методом левый прямоугольников:", left_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    elif modification == 2:
        print("Методом правых прямоугольников:", right_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    elif modification == 3:
        print("Методом средних прямоугольников:", medium_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    elif modification == 4:
        print("Методом левый прямоугольников:", left_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
        print("Методом правых прямоугольников:", right_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
        print("Методом средних прямоугольников:", medium_rectangle_method(a, b, eps, f), "Количество интервалов:", n)

elif method == 2:
    print("Методом трапеций:", trapezoid_method(a, b, eps, f), "Количество интервалов:", n)
elif method == 3:
    print("Методом Симпсона:", simpson_method(a, b, eps, f), "Количество интервалов:", n)
elif method == 4:
    print("Методом левый прямоугольников:", left_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    print("Методом правых прямоугольников:", right_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    print("Методом средних прямоугольников:", medium_rectangle_method(a, b, eps, f), "Количество интервалов:", n)
    print("Методом трапеций:", trapezoid_method(a, b, eps, f), "Количество интервалов:", n)
    print("Методом Симпсона:", simpson_method(a, b, eps, f), "Количество интервалов:", n)

print("Точное значение интеграла =", quad(f, a, b))