#define _CRT_SECURE_NO_WARNINGS // ��� fopen � visual studio
#include <cmath>
#include <cstdio>

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
	double e1 = 0, e2 = 0, einf = 0; // ���������� ����������� � L1, L2 � Linf ������
	double f1 = 0, f2 = 0, finf = 0; // ������������� ����������� = ��������� ����������� \ ����� f

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
			e1 += abs(tmp - f(x)); // L1 ����� ������� ������� � �������� �������� � ����� x
			f1 += abs(f(x)); // L1 ����� ������� � ����� x

			e2 += (tmp - f(x)) * (tmp - f(x)); // ����� ��������� ������� ������� � �������� �������� � ����� x
			f2 += f(x) * f(x); // ����� ��������� ������� � ����� x

			if (tmp - f(x) > einf) { // Linf - ����� ��������, ������ �� ���������
				einf = tmp - f(x); // ������� ������� � �������� �������� � ����� x
			}
			if (f(x) > finf) {
				finf = f(x); // ������� 
			}
		}
	}

	e2 = sqrt(e2); // ������ ����� ������ �� ����� ��������� - L2 ����� ��� ������� ������� � �������� �������� � ����� x
	f2 = sqrt(f2); // L2 ����� ��� �������

	// ������������� ����������� = ��������� ����������� \ ����� f


	// ������ � �����

	FILE* ResudialFile; // ���� ������������ - ������� txt ���� � ��������
	ResudialFile = fopen("Resudial.txt", "w");

	fprintf(ResudialFile, "|----------------|----------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|      mesh      |     ||*||1     |     ||*||2     |    ||*||inf    |\n");
	fprintf(ResudialFile, "|----------------|----------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|  h / 100 | abs | %14.8e | %14.8e | %14.8e |\n", e1, e2, einf);
	fprintf(ResudialFile, "|          | rel | %14.8e | %14.8e | %14.8e |\n", e1 / f1, e2 / f2, einf / finf);
	fprintf(ResudialFile, "|----------|-----|----------------|----------------|----------------|\n");
	fclose(ResudialFile);

	FILE *ParamsFile; // ���� ���������� - � ��� a,b � ��� �����
	ParamsFile = fopen("Params.txt","w");
	fprintf(ParamsFile,"%f, %f", a, b);
	fclose(ParamsFile);

	FILE *meshFile, *FmeshFile; // ����� ����� ������������ � �������� ������� � ���� ������ - ��� ��������� �������
	meshFile = fopen("mesh.txt", "w");
	FmeshFile = fopen("Fmesh.txt", "w");
	for (int i = 0; i < M - 1; i++) {
		fprintf(meshFile, "%f, ", mesh[i]);
		fprintf(FmeshFile, "%f, ", f(mesh[i]));
	}
	fprintf(meshFile, "%f\n", mesh[M - 1]);
	fprintf(FmeshFile, "%f\n", f(mesh[M - 1]));
	fclose(meshFile);
	fclose(FmeshFile);

	FILE *XFile, *LFile, *FFile; // ����� ����� ���������� ��������
	XFile = fopen("X.txt","w"); // ����� X
	LFile = fopen("L.txt","w"); // �������� �������� �������� � ���� ������
	FFile = fopen("F.txt","w"); // �������� ������� � ���� ������
	for (int i = 0; i < Mviz - 1; i++) {
		fprintf(XFile,"%f, ", X[i]);
		fprintf(FFile,"%f, ", f(X[i]));
		fprintf(LFile,"%f, ", L[i]);
	}
	fprintf(XFile, "%f\n", X[Mviz-1]);
	fprintf(FFile, "%f\n", f(X[Mviz-1]));
	fprintf(LFile, "%f\n", L[Mviz-1]);
	fclose(XFile);
	fclose(FFile);
	fclose(LFile);

	FILE *KFile, *KFFile;
	KFile = fopen("K.txt", "w"); // �������� ������� � ���� ������
	KFFile = fopen("KF.txt", "w"); // �������� ������� � ���� ������
	for (int i = 0; i < K - 1; i++) {
		fprintf(KFile, "%f, ", mesh[i * (N - 1)]);
		fprintf(KFFile, "%f, ", f(mesh[i * (N - 1)]));
	}
	fprintf(KFile, "%f\n", mesh[M - 1]);
	fprintf(KFFile, "%f\n", f(mesh[M - 1]));
	fclose(KFile);
	fclose(KFFile);
	
	std::system("python plot.py"); // ��� ������� �������� ��������� ������ � �������� ����������� ����� ������
	std::system("del /s /q Params.txt X.txt F.txt L.txt mesh.txt Fmesh.txt K.txt KF.txt"); // ������� ������ �����
	return 0;
}