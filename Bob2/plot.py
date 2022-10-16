import numpy as np # подключаем numpy для работы с векторами
import matplotlib.pyplot as plt # подключаем matplotlib для построения графиков

Params = np.genfromtxt("Params.txt", delimiter=",")
mesh = np.genfromtxt("mesh.txt", delimiter=",")
Fmesh = np.genfromtxt("Fmesh.txt", delimiter=",")

L = np.genfromtxt("L.txt", delimiter=",") # Вводим в numpy массив последовательность из файла, delimiter - разделитель в файле
X = np.genfromtxt("X.txt", delimiter=",")
F = np.genfromtxt("F.txt", delimiter=",")



#plt.figure() # создаем пустой график
#plt.rcParams['figure.figsize'] = (8, 4)
plt.figure(figsize=(16*2, 9*2))
plt.plot(X,F, color = 'C0', label = "f(x)") # добавляем на него график исходной функции
plt.scatter(mesh,Fmesh, color = 'red', s = 3) # добавляем на него график исходной функции
#label добавляет описание в легенду
plt.plot(X,L, color = 'C1', label = 'L(x)') # теперь график многочлена Лагранжа, o добавляет точки, r красный, "-" между ними говорит, что цвет относится к графику, а не к точкам
plt.ylabel('y') # подписываем ось x
plt.xlabel('x') # подписываем ось y
plt.xlim([Params[0], Params[1]]) # рисуем график на отрезке [a,b]
plt.legend() # выводим легенду
plt.grid(True)
plt.savefig('plot.png')
plt.show() # рисуем график на экране