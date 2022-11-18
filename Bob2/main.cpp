#include <iostream>
#include <fstream>

double f(double x) {
	return pow(x - 1, 4) + pow(x - 3, 3) + pow(x - 2, 2) + 200 * sin(10*x);
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
	// Residual - �������
	double resABS_1norm = 0, resABS_2norm = 0, resABS_infnorm = 0; // ���������� ����������� � L1, L2 � Linf ������
	double resREL_1norm = 0, resREL_2norm = 0, resREL_infnorm = 0; // ������������� ����������� � L1, L2 � Linf ������
	double f_1norm = 0, f_2norm = 0, f_infnorm = 0; // ������������� ����������� = ��������� ����������� \ ����� f

	int iviz = 0; // ������� ��� �������-�������
	for (int i = 0; i < K; i++) {
		double D[N];
		for (int j = 0; j < N; j++) { // ���������� ����������� ��� ���������� �������� �� ���������
			D[j] = 1;
			for (int k = 0; k < N; k++) {
				if(k!=j){
					// mesh[i * (N - 1) + j] - j-�� ����� �� i-�� ���������
					// mesh[i * (N - 1) + k] - k-�� ����� �� i-�� ���������
					D[j] *= mesh[i * (N - 1) + j] - mesh[i * (N - 1) + k]; // ��������� ���������� ��-�� ����, ��� ��������� �� ������ ����������
				}
			}
		}
		// ������ X[Mviz] � mesh[M] ������ ������ �� ���������
		// ������� ��� ������� �� K ����������
		// ������ i-��� �������� mesh[i * (N - 1)]
		// ����� i-��� ��������  mesh[(i + 1) * (N - 1)]
		// ����� ����� ����������� � ���������������� ��������� ������ �� �����,
		// � ������� X[iviz] <= mesh[(i + 1) * (N - 1)], �� �������� � ���� i-�� ��������
		// + ����������� �� iviz < Mviz - �� ������, ���� �� ����� ����� �� �������
		for (; X[iviz] <= mesh[(i + 1) * (N - 1)] && iviz < Mviz; iviz++) {
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(X[iviz]);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						// X[viz] - ����� � ������� ������� �������� ��������
						// mesh[i * (N - 1) + k] - k-�� ����� �� i-�� ���������
						multiplication *= (X[iviz] - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			L[iviz] = tmp; // ������� �������� �������� �������� � ����� X[iviz] - ���������� ������ ��� � ������ ������
		}

		// ������� �����
		for (int l = 0; l < 100; l++) { // ��� h/100
			double x = mesh[i * (N - 1)] + l * h / 100; // ����� x - �����-�� ������������� ����� �� K-�� �������
			double tmp = 0;
			for (int j = 0; j < N; j++) {
				double multiplication = f(x);
				for (int k = 0; k < N; k++) {
					if (k != j) {
						// ����������
						// mesh[i * (N - 1) + k] - k-�� ����� �� i-�� ���������
						multiplication *= (x - mesh[i * (N - 1) + k]);
					}
				}
				tmp += multiplication / D[j];
			}
			// � ����� tmp - �������� �������� � ����� x
			// ��������� ����� ���������
			// �� ���� � �������� ����� (� ������ L2 ������ ����� ���������, ������ ����� ������)
			// �� ������ �� ���������� � ���������� ���������
			resABS_1norm += abs(tmp - f(x)); // L1 ����� ������� ������� � �������� �������� � ����� x
			f_1norm += abs(f(x)); // L1 ����� ������� � ����� x

			resABS_2norm += (tmp - f(x)) * (tmp - f(x)); // ����� ��������� ������� ������� � �������� �������� � ����� x
			f_2norm += f(x) * f(x); // ����� ��������� ������� � ����� x

			if (tmp - f(x) > resABS_infnorm) { // Linf - ����� ��������, ������ �� ���������
				resABS_infnorm = tmp - f(x); // ������� ������� � �������� �������� � ����� x
			}
			if (f(x) > f_infnorm) {
				f_infnorm = f(x); // ������� 
			}
		}
	}

	resABS_2norm = sqrt(resABS_2norm); // ������ ����� ������ �� ����� ��������� - L2 ����� ��� ������� ������� � �������� �������� � ����� x
	f_2norm = sqrt(f_2norm); // L2 ����� ��� �������

	// ������������� ����������� = ��������� ����������� \ ����� f

	resREL_1norm = resABS_1norm / f_1norm;
	resREL_2norm = resABS_2norm / f_2norm;
	resREL_infnorm = resABS_infnorm / f_infnorm;


	// ������ � �����

	std::ofstream ParamsFile; // ���� ���������� - � ��� a,b � ��� �����
	ParamsFile.open("Params.txt");
	ParamsFile << a << ", " << b << ", "
		<< resREL_1norm << ", " << resREL_2norm << ", " << resREL_infnorm << ", "
		<< resABS_1norm << ", " << resABS_2norm << ", " << resABS_infnorm << std::endl;
	ParamsFile.close();

	std::ofstream meshFile, FmeshFile; // ����� ����� ������������ � �������� ������� � ���� ������ - ��� ��������� �������
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

	std::ofstream XFile, LFile, FFile, KFile, KFFile; // ����� ����� ���������� ��������
	XFile.open("X.txt"); // ����� X
	LFile.open("L.txt"); // �������� �������� �������� � ���� ������
	FFile.open("F.txt"); // �������� ������� � ���� ������
	KFile.open("K.txt"); // �������� ������� � ���� ������
	KFFile.open("KF.txt"); // �������� ������� � ���� ������
	for (int i = 0; i < Mviz - 1; i++) {
		XFile << X[i] << ", ";
		FFile << f(X[i]) << ", ";
		LFile << L[i] << ", ";
	}
	for (int i = 0; i < K - 1; i++) {
		KFile << mesh[i * (N - 1)] << ", ";
	}
	for (int i = 0; i < K - 1; i++) {
		KFFile << f(mesh[i * (N - 1)]) << ", ";
	}
	XFile << X[Mviz - 1] << std::endl;
	FFile << f(X[Mviz - 1]) << std::endl;
	LFile << L[Mviz - 1] << std::endl;
	KFile << mesh[M - 1] << std::endl;
	KFFile << f(mesh[M - 1]) << std::endl;

	XFile.close();
	FFile.close();
	LFile.close();
	KFile.close();
	KFFile.close();
	
	std::system("python plot.py"); // ��� ������� �������� ��������� ������ � �������� ����������� ����� ������

	return 0;
}