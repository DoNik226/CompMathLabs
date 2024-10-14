import numpy as np
from scipy import linalg
def get_triangular(matrix):
    n = len(matrix)
    per=0

    for i in range(n):
        for j in range(n):
            index_mass[i][j] = j + 1

    for i in range(n-1):
        point=0
        max_val = 0
        for j in range(i, n):
            if abs(matrix[j][i]) > max_val:
                max_val = abs(matrix[j][i])
                point=j

        if max_val == 0:
            return None, 0

        if matrix[i][i]==0: #!!!!!!
            matrix[i], matrix[point] = matrix[point], matrix[i]
            index_mass[i], index_mass[point] = index_mass[point], index_mass[i]
            print("Переставляем строки " , i + 1 , " и " , point + 1)
            per+=1
            print("Матрица после перестановки:")
            print_matrix(matrix)
        else:
            print("Перестановка не требуется");



        for k in range(i + 1, n):
            for j in range(n, i-1, -1):
                matrix[k][j] = matrix[k][j] - (matrix[k][i]/matrix[i][i])*matrix[i][j]

        print("Матрица после " , (i + 1) , "го преобразования:");
        print_matrix(matrix);
        print("-------");
    return matrix, per


def get_determinant(a, per):
    det = 1


    size = len(a)
    for i in range(size):
        det *= a[i][i]
    if per%2!=0:
        det*=(-1)

    return det


def get_roots(matrix):
    n = len(matrix)

    roots = [0]*(size)

    write_out_cnt = n - 1
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            roots[index_mass[n - 1][write_out_cnt] - 1] = matrix[i][n] / matrix[i][n - 1]
        else:
            root = matrix[i][n]
            point = 0
            for j in range(n - 1, -1, -1):
                if matrix[i][j]!=0:
                    if roots[index_mass[n - 1][j] - 1] != 0:
                        root -= matrix[i][j] * roots[index_mass[n - 1][j] - 1]
                    point = j
            roots[index_mass[n - 1][write_out_cnt] - 1] = root / matrix[i][point]

        write_out_cnt -= 1

    return roots


def get_discrepancy(matrix, x):
    n = len(matrix)


    dis = [0]*size

    for i in range(n):
        r = matrix[i][n]
        for j in range(n):
            r -= matrix[i][j] * x[j]

        dis[i] = r

    return dis


def print_matrix(a):
    size = len(a)


    for row in a:
        for j in range(size + 1):
            print("{:.3f}".format(row[j]), end=' ')
        print()
    print()

#--------------------------------------------------------------------------------------------------

size = 0
perestanovki=0
arrayList = []

size = int(input("Укажите размерность матрицы: "))

if size == 1:
    print("Размерность СЛАУ не может быть равна одному")
elif size == 2:
    print("Формат ввода: 'a11 a12 b1'")
    print("Введите коэффициенты через пробел:")
else:
    print(f"Формат ввода: 'a11 ... a1{size} b1'")
    print("Введите коэффициенты через пробел:")

try:
    for i in range(size):
        inputs = input().split()
        for i in range(len(inputs)):
            arrayList.append(float(inputs[i]))
except ValueError:
    print("Ошибка ввода! Проверьте, что дробные числа записаны через точку")

index_mass = [[0]*(size + 1) for i in range(size)]
mtx = [[0]*(size + 1) for i in range(size)]
A = [[0]*(size) for i in range(size)]
index = 0
for i in range(size):
    for j in range(size + 1):
        mtx[i][j] = arrayList[index]
        index += 1
index=0
for i in range(size):
    for j in range(size):
        A[i][j] = arrayList[index]
        index += 1
    index+=1
print_matrix(mtx)
triangleMtx, perestanovki = get_triangular(mtx)
if triangleMtx != None:
    print("Получена треугольная матрица:")
    print_matrix(triangleMtx)

    print("Количество перестановок:", perestanovki)
    print("Определитель матрицы равен:")
    det = get_determinant(triangleMtx, perestanovki)
    print(det)
    print()

    if det != 0:
        x = get_roots(triangleMtx)
        print("Найдены корни СЛАУ:")
        for v in x:
            if v==-0:
                v=0
            print("{:.3f}".format(v), ' ')
        print()
        print()

        print("Вектор невязки:")
        dis = get_discrepancy(mtx, x)
        for di in dis:
            print(di, ' ')
        print()

    else:
        print("Система не имеет решений! Определитель = 0")
else:
    print("Система не имеет решений! Определитель = 0")
print("Определитель, посчитанный библиотекой = ", np.linalg.det(A))