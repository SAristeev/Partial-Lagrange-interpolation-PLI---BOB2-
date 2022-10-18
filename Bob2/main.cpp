#include <iostream>
#include <fstream>

double f(double x) {
	return pow(x - 1, 4) + pow(x - 3, 3) + pow(x - 2, 2) + 100 * sin(x);
}
// K ����������
// �� ������ ��������� ������� �������� �������� N
// ��� ���������� ����� N ����� �� ������ �� ����������
// � ������ ��������� ����� ����� ���������� ����� ������������ 
// M = K * (N - 1) + 1
int main(){
	double a = -2, b = 2;
	const int K = 100; // ���������� ����������
	const int N = 10; // ������� ���������������� ����������� ��������
	const int M = K * (N - 1) + 1; // ����� ���������� ����� ������������
	double h = (b - a) / (M - 1); // ��������� ��� ����������� �����, ����� � ��� �� 1 ������, ��� ����������, ������� M - 1

	double mesh[M]; // ������ ����������� �����
	for (int i = 0; i < M; i++) {
		mesh[i] = a + i * h; // ��������� ������ ���� ����������� �����
	}

	const int Mviz = 1920; // ���������� �����, �� ������� ����� ������� ������
	double X[Mviz]; // ������ ���� �����
	double L[Mviz]; // ����������� � ���� ������
	for (int i = 0; i < Mviz; i++) {
		X[i] = a + i * (b - a) / (Mviz - 1); // ��������� �����, �� ������� ������ ������
	}
	double resABS_1norm = 0, resABS_2norm = 0, resABS_infnorm = 0;
	double resREL_1norm = 0, resREL_2norm = 0, resREL_infnorm = 0;
	double f_1norm = 0, f_2norm = 0, f_infnorm = 0;
	int iviz = 0; // ������� ��� �������-�������
	for (int i = 0; i < K; i++) {
		double D[N];
		for (int j = 0; j < N; j++) { // ���������� ����������� ��� ���������� �������� �� ���������
			D[j] = 1;
			for (int k = 0; k < N; k++) {
				if(k!=j){
					D[j] *= mesh[i * (N - 1) + j] - mesh[i * (N - 1) + k]; // ��������� ���������� ��-�� ����, ��� ��������� �� ������ ����������
				}
			}
		}
	
		for (; X[iviz] <= mesh[(i + 1) * (N - 1)] && iviz < Mviz; iviz++) {
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(X[iviz]);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						multiplication *= (X[iviz] - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			L[iviz] = tmp;
		}

		for (int l = 0; l < 100; l++) {
			double x = mesh[i * (N - 1)] + l * h / 100;
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(x);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						multiplication *= (x - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			resABS_1norm += abs(tmp - f(x));
			f_1norm += abs(f(x));
			resABS_2norm += (tmp - f(x)) * (tmp - f(x));
			f_2norm += f(x) * f(x);
			if (tmp - f(x) > resABS_infnorm) {
				resABS_infnorm = tmp - f(x);
			}
			if (f(x) > f_infnorm) {
				f_infnorm = f(x);
			}
		}
	}
	resABS_2norm = sqrt(resABS_2norm);
	f_2norm = sqrt(f_2norm);
	resREL_1norm = resABS_1norm / f_1norm;
	resREL_2norm = resABS_2norm / f_2norm;
	resREL_infnorm = resABS_infnorm / f_infnorm;


	//
	// ++++++++++++++++++++++++++++++++++++=
	//
	std::ofstream ParamsFile;
	ParamsFile.open("Params.txt");
	ParamsFile << a << ", " << b << ", "
		<< resREL_1norm << ", " << resREL_2norm << ", " << resREL_infnorm << ", "
		<< resABS_1norm << ", " << resABS_2norm << ", " << resABS_infnorm << std::endl;
	ParamsFile.close();

	std::ofstream meshFile, FmeshFile;
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

	std::ofstream XFile, LFile, FFile;
	XFile.open("X.txt");
	LFile.open("L.txt");
	FFile.open("F.txt");
	for (int i = 0; i < Mviz - 1; i++) {
		XFile << X[i] << ", ";
		FFile << f(X[i]) << ", ";
		LFile << L[i] << ", ";
	}
	XFile << X[Mviz - 1] << std::endl;
	FFile << f(X[Mviz - 1]) << std::endl;
	LFile << L[Mviz - 1] << std::endl;
	XFile.close();
	LFile.close();
	
	std::system("python plot.py"); // ��� ������� �������� ��������� ������ � �������� ����������� ����� ������

	return 0;
}