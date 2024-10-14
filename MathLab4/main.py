import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

x: int
y: int
content = []
sko = []


def draw(f, a, z):
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, 'ko')
    plt.plot(a, f)
    plt.title(z)
    plt.grid()
    plt.show()
def cof_cor(Xs, Ys):
    X_mean = np.mean(Xs)
    Y_mean = np.mean(Ys)
    numerator = 0
    denominator_x = 0
    denominator_y = 0
    for x, y in zip(Xs, Ys):
        numerator += (x-X_mean)*(y-Y_mean)
        denominator_x += (x-X_mean) ** 2
        denominator_y += (y-Y_mean) ** 2
    return numerator / sqrt(denominator_x * denominator_y)

def lin_approx(Xs, Ys):
    A = np.array([[sum(Xs ** 2), sum(Xs)],[sum(Xs), len(Xs)]])
    B = np.array([sum(Xs*Ys), sum(Ys)])
    a, b = np.linalg.solve(A, B)
    return a, b

def linear(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    #Подсчет аппроксимации
    a, b = lin_approx(Xs, Ys)
    print(f'Коэффициент корреляции = {cof_cor(Xs, Ys)}')
    f = a * Xs + b
    #СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    #Коэф_детер
    sumf=sum(f)/len(f)
    deter = 1 - sum(eps**2)/sum((Ys-sumf)**2)
    #мера отклонения
    S = sum(eps ** 2)
    print(f'Линейная аппроксимация\nf = {a} * x + {b}')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(a * i + b), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * Xs[i] + b - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F=[]
    X=[]
    h=(Xs[-1]+1-Xs[0]+1)/1000
    i=0
    while i<=1000:
        X.append(Xs[0]-1+i*h)
        F.append(a * X[i] + b)
        i+=1
    draw(F, X, "Линейная функция")

def squared(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    #Подсчет аппроксимации
    A = np.array([[len(Xs), sum(Xs), sum(Xs ** 2)],
                  [sum(Xs), sum(Xs ** 2), sum(Xs ** 3)],
                  [sum(Xs ** 2), sum(Xs ** 3), sum(Xs ** 4)]])
    B = np.array([sum(Ys), sum(Xs*Ys), sum((Xs**2)*Ys)])
    c, b, a = np.linalg.solve(A, B)
    f = a * Xs ** 2 + b * Xs + c
    # СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    # Коэф_детер
    sumf = sum(f) / len(f)
    deter = 1 - sum(eps ** 2) / sum((Ys - sumf) ** 2)
    #Мера отклонения
    S = sum(eps ** 2)
    print(f'\n\nКвадратичная аппроксимация\nf = {a} * x^2 + {b} * x +{c}')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(a * i ** 2 + b * i + c), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * Xs[i] ** 2 + b * Xs[i] + c - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F = []
    X = []
    h = (Xs[-1] + 1 - Xs[0] + 1) / 1000
    i = 0
    while i <= 1000:
        X.append(Xs[0]-1 + i * h)
        F.append(a * X[i] ** 2 + b * X[i] + c)
        i += 1
    draw(F, X, "Полиномиальная второй степени")

def triple(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    #Подсчет аппроксимации
    A = np.array([[len(Xs), sum(Xs), sum(Xs ** 2), sum(Xs ** 3)],
                  [sum(Xs), sum(Xs ** 2), sum(Xs ** 3), sum(Xs ** 4)],
                  [sum(Xs ** 2), sum(Xs ** 3), sum(Xs ** 4), sum(Xs ** 5)],
                  [sum(Xs ** 3), sum(Xs ** 4), sum(Xs ** 5), sum(Xs ** 6)]])
    B = np.array([sum(Ys), sum(Xs*Ys), sum((Xs**2)*Ys), sum((Xs **3 )*Ys)])
    d, c, b, a = np.linalg.solve(A, B)
    f = a * Xs ** 3 + b * Xs ** 2 + c * Xs + d
    # СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    # Коэф_детер
    sumf = sum(f) / len(f)
    deter = 1 - sum(eps ** 2) / sum((Ys - sumf) ** 2)
    #Мера отклонения
    S = sum(eps ** 2)
    print(f'\n\nТретичная аппроксимация\nf = {a} * x^3 + {b} * x^2 + {c} * x + {d}')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(a * i ** 3 + b * i ** 2 + c * i + d), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * Xs[i] ** 3 + b * Xs[i] ** 2 + c * Xs[i] + d - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F = []
    X = []
    h = (Xs[-1] + 1 - Xs[0] + 1) / 1000
    i = 0
    while i <= 1000:
        X.append(Xs[0]-1 + i * h)
        F.append(a * X[i] ** 3 + b * X[i] ** 2 + c * X[i] + d)
        i += 1
    draw(F, X, "Полиномиальная третьей степени")
    return a, b, c, S, delta

def power(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    log_Xs = np.log(Xs)
    log_Ys = np.log(Ys)
    if True in np.isnan(log_Xs) or True in np.isnan(log_Ys):
        raise ValueError
    #Подсчет аппроксимации
    b, a = lin_approx(log_Xs, log_Ys)
    a = np.exp(a)
    f = a * (Xs ** b)
    # СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    # Коэф_детер
    sumf = sum(f) / len(f)
    deter = 1 - sum(eps ** 2) / sum((Ys - sumf) ** 2)
    #Мера отклонения
    S = sum(eps ** 2)
    print(f'\n\nСтепенная аппроксимация\nf = {a} * x ** {b}')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(a * (i ** b)), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * (Xs[i] ** b) - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F = []
    X = []
    h = (Xs[-1] + 1 - Xs[0] + 1) / 1000
    i = 0
    while i <= 1000:
        X.append(Xs[0]-1 + i * h)
        F.append(a * (X[i] ** b))
        i += 1
    draw(F, X, "Степенная функция")

def exponential(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    log_Ys = np.log(Ys)
    if True in np.isnan(log_Ys):
        raise ValueError
    #Подсчет аппроксимации
    b, a = lin_approx(Xs, log_Ys)
    a = np.exp(a)
    f = a * (np.exp(Xs * b))
    # СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    # Коэф_детер
    sumf = sum(f) / len(f)
    deter = 1 - sum(eps ** 2) / sum((Ys - sumf) ** 2)
    #Мера отклонения
    S = sum(eps ** 2)
    print(f'\n\nЭкспоненциальная аппроксимация\nf = {a} * e**(x * {b})')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(a * (np.exp(i * b))), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * (np.exp(Xs[i] * b)) - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F = []
    X = []
    h = (Xs[-1] + 1 - Xs[0] + 1) / 1000
    i = 0
    while i <= 1000:
        X.append(Xs[0]-1 + i * h)
        F.append(a * (np.exp(X[i] * b)))
        i += 1
    draw(F, X, "Экспоненциальная функция")
    return a, b, '-', S, delta

def logarithm(array):
    Xs = np.array(array[0])
    Ys = np.array(array[1])
    log_Xs = np.log(Xs)
    if True in np.isnan(log_Xs):
        raise ValueError
    #Подсчет аппроксимации
    a, b = lin_approx(log_Xs, Ys)
    f = a * (log_Xs) + b
    # СКО
    eps = Ys - f
    delta = (sum(eps ** 2) / len(Xs)) ** 0.5
    sko.append(delta)
    # Коэф_детер
    sumf = sum(f) / len(f)
    deter = 1 - sum(eps ** 2) / sum((Ys - sumf) ** 2)
    #Мера отклонения
    S = sum(eps ** 2)
    print(f'\n\nЛогарифмическая аппроксимация\nf = {a} * ln(x) + {b}')
    print("X:     ", end=' | ')
    for i in Xs:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Y:     ", end=' | ')
    for i in Ys:
        print("{:10.3f}".format(i), end=' | ')
    print()
    print("Phi(X):", end=' | ')
    for i in log_Xs:
        print("{:10.3f}".format(a * (i) + b), end=' | ')
    print()
    print("Eps:   ", end=' | ')
    for i in range(len(Xs)):
        print("{:10.3e}".format(a * (log_Xs[i]) + b - Ys[i]), end=' | ')
    print()
    print(f'Мера отклонения = {S}')
    print(f'СКО = {delta}\n\n')
    print(f'Коэффициент детерминации = {deter}\n\n')
    F = []
    X = []
    h = (Xs[-1] + 1 - Xs[0] + 1) / 1000
    i = 0
    s=1
    if Xs[0]<=1:
        s=Xs[0]*0.999
    while i <= 1000:
        X.append(Xs[0]-s + i * h)
        F.append(a * (np.log(X[i])) + b)
        i += 1
    draw(F, X, "Логарифмическая функция")
    return a, b, '-', S, delta


def run():
    again = True
    while again:
        again = False
        in_type = input('Введите:\n\t* k - если вводить с клавиатуры\n\t* f - если хотите вводить из файла\n')
        if in_type.strip() == 'k':
            line_x = input("Введите значиния x (от 8 до 12 значений: ")
            content.append([float(x) for x in line_x.split(" ")])
            line_y = input("Введите значения y (от 8 до 12 значений: ")
            content.append([float(x) for x in line_y.split(" ")])
        elif in_type.strip() == 'f':
            with open("input1.txt") as f:
                for line in f:
                    content.append([float(x) for x in line.split(" ")])
            f.close()
        else:
            print('Введено неверно, попробуйте снова.')
            again = True
        if not again:
            if len(content[0])!=len(set(content[0])):
                print("Среди введенных х есть повторяющиеся, введите данные заново")
                again=True
            elif len(content[0])!=len(content[1]) or len(content[0])<8 or len(content[0])>12:
                print("Введено неправильное количество данных. Количество х и у должно совпадать и содержать от 8 до 12 значений")
                again = True
        if again:
            content.clear()

    global x
    x = np.array(content[0])
    global y
    y = np.array(content[1])
    for i in range(len(content[0])-1, 0, -1):
        for j in range(i):
            if content[0][j]>content[0][j+1]:
                content[0][j], content[0][j+1]=content[0][j+1],content[0][j]
                content[1][j], content[1][j+1]=content[1][j+1], content[1][j]
    linear(content)
    squared(content)
    triple(content)
    if min(content[0]) <= 0:
        print("Среди введенных x есть те, которые недопустимы для степенной функции, введите x>0\n\n")
    elif min(content[1]) <= 0:
        print("Среди введенных y есть те, которые недопустимы для степенной функции, введите y>0\n\n")
    else:
        power(content)
    if min(content[1]) <= 0:
        print("Среди введенных y есть те, которые недопустимы для экспоненциальной функции, введите y>0\n\n")
    else:
        exponential(content)
    if min(content[0]) <= 0:
        print("Среди введенных x есть те, которые недопустимы для логарифмической функции, введите x>0\n\n")
    else:
        logarithm(content)

    minsko=0
    for i in range(6):
        if min(sko)==sko[i]:
            minsko=i
            break
    if (minsko==0):
        print("Наилучшая аппроксимирующая функция: линейная функция")
    elif (minsko==1):
        print("Наилучшая аппроксимирующая функция: полиномиальная функция 2-й степени")
    if (minsko==2):
        print("Наилучшая аппроксимирующая функция: полиномиальная функция 3-й степени")
    if (minsko==3):
        print("Наилучшая аппроксимирующая функция: степенная функция")
    if (minsko==4):
        print("Наилучшая аппроксимирующая функция: экспоненциальная функция")
    if (minsko==5):
        print("Наилучшая аппроксимирующая функция: логарифмическая функция")
    print("СКО =", sko[minsko])
run()

