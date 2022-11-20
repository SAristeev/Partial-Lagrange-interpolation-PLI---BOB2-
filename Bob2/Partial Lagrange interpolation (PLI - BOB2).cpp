#include <iostream>
#include <iomanip>
#include <fstream>

double f(double x) {
	return pow(x - 1, 4) + pow(x - 3, 3) + pow(x - 2, 2) + 200 * sin(10*x);
}
// K интервалов
// На каждом интервале степень полинома Лагранжа N
// Для построения нужно N точек на каждом из интервалов
// С учетом граничных точек общее количество узлов интерполяции 
// M = K * (N - 1) + 1
int main(){
	double a = -2, b = 2;
	const int K = 100; // Количество интервалов
	const int N = 10; // Степень интерполяционных многочленов Лагранжа
	const int M = K * (N - 1) + 1; // Общее количество узлов интерполяции
	double h = (b - a) / (M - 1); // Вычисляем шаг равномерной сетки, узлов у нас на 1 больше, чем интервалов, поэтому M - 1

	double mesh[M]; // Массив равномерной сетки
	for (int i = 0; i < M; i++) {
		mesh[i] = a + i * h; // вычисляем каждый узел равномерной сетки
	}

	const int Mviz = 1920; // количество точек, по которым будем строить график
	double X[Mviz]; // массив этих точек
	double L[Mviz]; // приближение в этих точках
	for (int i = 0; i < Mviz; i++) {
		X[i] = a + i * (b - a) / (Mviz - 1); // заполняем точки, по которым строим график
	}
	// Residual - остаток
	double e1 = 0, e2 = 0, einf = 0; // абсолютная погрешность в L1, L2 и Linf нормах
	double f1 = 0, f2 = 0, finf = 0; // относительная погрешность = абсоютная погрешность \ норма f

	int iviz = 0; // Счетчик для массива-графика
	for (int i = 0; i < K; i++) {
		double D[N];
		for (int j = 0; j < N; j++) { // Вычисление знаменателя для многочлена Лагранжа на интервале
			D[j] = 1;
			for (int k = 0; k < N; k++) {
				if(k!=j){
					// mesh[i * (N - 1) + j] - j-ая точка на i-ом интервале
					// mesh[i * (N - 1) + k] - k-ая точка на i-ом интервале
					D[j] *= mesh[i * (N - 1) + j] - mesh[i * (N - 1) + k]; // Смещенная индексация из-за того, что вычисляем на разных интервалах
				}
			}
		}
		// Массив X[Mviz] и mesh[M] вообще говоря не совпадают
		// Поэтому для каждого из K интервалов
		// начало i-ого инервала mesh[i * (N - 1)]
		// конец i-ого инервала  mesh[(i + 1) * (N - 1)]
		// Тогда будем подставлять в интерполяционный многочлен только те точки,
		// у которых X[iviz] <= mesh[(i + 1) * (N - 1)], те попадают в этот i-ый интервал
		// + ограничение на iviz < Mviz - на случай, если мы можем выйти из массива
		for (; X[iviz] <= mesh[(i + 1) * (N - 1)] && iviz < Mviz; iviz++) {
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(X[iviz]);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						// X[viz] - точка в котором считаем значение полинома
						// mesh[i * (N - 1) + k] - k-ая точка на i-ом интервале
						multiplication *= (X[iviz] - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			L[iviz] = tmp; // подсчет значения полинома Лагранжа в точке X[iviz] - аналогично задаче как в первой задаче
		}

		// подсчет нормы
		for (int l = 0; l < 100; l++) { // шаг h/100
			double x = mesh[i * (N - 1)] + l * h / 100; // здесь x - какая-то фиксированная точка на K-ом отрезке
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(x);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						// аналогично
						// mesh[i * (N - 1) + k] - k-ая точка на i-ом интервале
						multiplication *= (x - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			// в итоге tmp - значение полинома в точке x
			// Вычисляем норму аддитивно
			// То есть я вычисляю норму (в случае L2 только сумму квадратов, корень потом возьму)
			// На каждом из интервалов и складываем результат
			e1 += abs(tmp - f(x)); // L1 норма разницы функции и полинома Лагранжа в точке x
			f1 += abs(f(x)); // L1 норма функции в точке x

			e2 += (tmp - f(x)) * (tmp - f(x)); // Сумма квадратов разницы функции и полинома Лагранжа в точке x
			f2 += f(x) * f(x); // Сумма квадратов функции в точке x

			if (tmp - f(x) > einf) { // Linf - берем максимум, ничего не суммируем
				einf = tmp - f(x); // разница функции и полинома Лагранжа в точке x
			}
			if (f(x) > finf) {
				finf = f(x); // функция 
			}
		}
	}

	e2 = sqrt(e2); // Теперь берем корень из суммы квадратов - L2 норма для разницы функции и полинома Лагранжа в точке x
	f2 = sqrt(f2); // L2 норма для функции

	// относительная погрешность = абсоютная погрешность \ норма f


	// Запись в файлы

	std::ofstream ResudialFile; // файл параметров - в нем a,b и все нормы
	ResudialFile.open("Resudial.txt");
	ResudialFile << std::setiosflags(std::ios_base::scientific);
	ResudialFile << "|----------------|--------------|--------------|--------------|" << std::endl;
	ResudialFile << "|      mesh      |    ||*||1    |    ||*||2    |   ||*||inf   |" << std::endl;
	ResudialFile << "|----------------|--------------|--------------|--------------|" << std::endl;

	ResudialFile << "|  h / 100 | abs | "
		<< e1 << " | "
		<< e2 << " | "
		<< einf << " |" << std::endl;
	ResudialFile << "|----------|-----|--------------|--------------|--------------|" << std::endl;

	ResudialFile << "|          | rel | "
		<< e1 / f1 << " | "
		<< e2 / f2 << " | "
		<< einf / finf << " |" << std::endl;
	ResudialFile << "|----------|-----|--------------|--------------|--------------|" << std::endl;
	ResudialFile.close();


	std::ofstream ParamsFile; // файл параметров - в нем a,b и все нормы
	ParamsFile.open("Params.txt");
	ParamsFile << a << ", " << b << ", "
		<< e1 / f1 << ", " << f2 / e2 << ", " << einf / finf << ", "
		<< e1 << ", " << e2 << ", " << einf << std::endl;
	ParamsFile.close();

	std::ofstream meshFile, FmeshFile; // файлы точек интерполяции и значения функции в этих точках - для выделения красным
	meshFile.open("mesh.txt");
	FmeshFile.open("Fmesh.txt");
	for (int i = 0; i < M - 1; i++) {
		meshFile << mesh[i] << ", ";
		FmeshFile << f(mesh[i]) << ", ";
	}
	meshFile << mesh[M - 1] << std::endl;
	FmeshFile << f(mesh[M - 1]) << std::endl;
	meshFile.close();
	FmeshFile.close();

	std::ofstream XFile, LFile, FFile; // файлы точек построения графиков
	XFile.open("X.txt"); // Точки X
	LFile.open("L.txt"); // Значения полинома Лагранжа в этих точках
	FFile.open("F.txt"); // Значения функции в этих точках
	for (int i = 0; i < Mviz - 1; i++) {
		XFile << X[i] << ", ";
		FFile << f(X[i]) << ", ";
		LFile << L[i] << ", ";
	}
	XFile << X[Mviz - 1] << std::endl;
	FFile << f(X[Mviz - 1]) << std::endl;
	LFile << L[Mviz - 1] << std::endl;
	XFile.close();
	FFile.close();
	LFile.close();

	std::ofstream KFile, KFFile;
	KFile.open("K.txt"); // Значения функции в этих точках
	KFFile.open("KF.txt"); // Значения функции в этих точках
	for (int i = 0; i < K - 1; i++) {
		KFile << mesh[i * (N - 1)] << ", ";
	}
	for (int i = 0; i < K - 1; i++) {
		KFFile << f(mesh[i * (N - 1)]) << ", ";
	}
	KFile << mesh[M - 1] << std::endl;
	KFFile << f(mesh[M - 1]) << std::endl;
	KFile.close();
	KFFile.close();
	
	std::system("python plot.py"); // эта команда вызывает командную строку и включает питоновскую часть задачи
	std::system("del /s /q Params.txt X.txt F.txt L.txt mesh.txt Fmesh.txt K.txt KF.txt");
	return 0;
}