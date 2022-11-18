import numpy as np # подключаем numpy для работы с векторами
import matplotlib.pyplot as plt # подключаем matplotlib для построения графиков

# delimiter - разделитель в файле
Params = np.genfromtxt("Params.txt", delimiter=",") # читаем файл параметров
mesh = np.genfromtxt("mesh.txt", delimiter=",") # читаем массив точек интерполяции
Fmesh = np.genfromtxt("Fmesh.txt", delimiter=",") # читаем массив значений функции в точках интерполяции


X = np.genfromtxt("X.txt", delimiter=",") # читаем массив точек для графика
L = np.genfromtxt("L.txt", delimiter=",") # читаем массив значений полинома в этих точках
F = np.genfromtxt("F.txt", delimiter=",") # читаем массив значений функции в этих точках
K = np.genfromtxt("K.txt", delimiter=",") # читаем массив значений функции в этих точках
KF = np.genfromtxt("KF.txt", delimiter=",") # читаем массив значений границ

plt.figure(figsize=(16*2, 9*2)) # 16 на 9 - соотношение экрана, умножить на 2 - увеличиваю картинку
plt.xlim([Params[0], Params[1]]) # рисуем график на отрезке [a,b]
plt.plot(X,F, color = 'C0', label = "f(x)") # рисуем график функции
plt.scatter(mesh,Fmesh, color = 'red', s = 8) # выделяем красным точки интерополяции
plt.scatter(K,KF, color = 'blue', s = 16) # выделяем красным точки интерополяции
#label добавляет описание в легенду
plt.plot(X,L, color = 'C1', label = 'L(x)') # рисуем график полинома Лагранжа
plt.ylabel('y') # подписываем ось x
plt.xlabel('x') # подписываем ось y
plt.legend() # выводим легенду
plt.grid(True) # рисуем сетку
plt.savefig('plot.png') # сохраняем график в файл 'plot.png'
plt.show() # рисуем график на экране