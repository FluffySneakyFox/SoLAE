#include <stdio.h>
#include <windows.h>
#include <time.h>

#include <iostream>
#include <fstream>
using namespace std;

void GenerateMatrixFile(char* filename, int N, int maxElem)
{
	fstream stream(filename, ios::out);
	double* A = new double[N];
	for (int i = 0; i < N; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < N; j++)
		{
			A[j] = (double)(rand() % 1000) / 1000.0 * maxElem + 0.1;
			sum += A[j];
		}
		A[i] = sum + 1.0;

		for (int j = 0; j < N; j++) stream << A[j] << " ";
		stream << endl;

	}
	delete[]A;
	stream.close();
}
void GenerateVectorFile(char* filename, int N, int maxElem)
{
	fstream stream(filename, ios::out);
	for (int i = 0; i < N; i++) stream << (double)(rand() % 1000) / 1000.0 * maxElem + 0.1 << endl;
	stream.close();
}
void ReadMatrixFromFile(char* filename, int N, double** Matrix)
{
	fstream stream(filename, ios::in);
	for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) stream >> Matrix[i][j];
	stream.close();
}
void ReadVectorFromFile(char* filename, int N, double* Vector)
{
	fstream stream(filename, ios::in);
	for (int i = 0; i < N; i++) stream >> Vector[i];
	stream.close();
}
void WriteMatrixInFile(char* filename, int N, double** Matrix)
{
	fstream stream(filename, ios::out);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) stream << Matrix[i][j] << " ";
		stream << endl;
	}
	stream.close();
}
void WriteVectorInFile(char* filename, int N, double* Vector)
{
	fstream stream(filename, ios::out);
	for (int i = 0; i < N; i++) stream << Vector[i] << endl;
	stream.close();
}
void Gauss(int N, double** A, double* B, double* X)
{
	//Transform A into upper triangular matrix
	for (int k = 0; k < N - 1; k++)
	{
		for (int i = k + 1; i < N; i++)
		{
			double mult = A[i][k] / A[k][k];
			for (int j = 0; j < N; j++) A[i][j] -= mult * A[k][j];
			B[i] -= mult * B[k];
		}
	}
	//Calculate X'es
	X[N - 1] = B[N - 1] / A[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = 1; j < N - i; j++)sum += A[i][N - j] * X[N - j];
		X[i] = (B[i] - sum) / A[i][i];
	}
	//Output
	char filenameX[40]; sprintf_s(filenameX, 40, "VectorX%dGauss.txt", N);
	WriteVectorInFile(filenameX, N, X);
}
void Jacobi(int N, double** A, double* B, double* X, double epsilon)
{
	//First approximation of x
	for (int i = 0; i < N; i++) X[i] = B[i] / A[i][i];
	//Data copy of X
	double* XNew = new double[N];
	//Calc
	double maxDelta;
	do
	{
		//New X's
		for (int i = 0; i < N; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < N; j++) if (i != j) sum += A[i][j] * X[j];
			XNew[i] = (B[i] - sum) / A[i][i];
		}
		//Check stop condition by difference |Xnew-Xold|
		maxDelta = -1.0;
		for (int i = 0; i < N; i++) if (maxDelta < fabs(XNew[i] - X[i])) maxDelta = fabs(XNew[i] - X[i]);
		//Copy data between iterations
		for (int i = 0; i < N; i++) X[i] = XNew[i];
	} while (maxDelta >= epsilon);
	//Output
	char filenameX[40]; sprintf_s(filenameX, 40, "VectorX%dJacobi.txt", N);
	WriteVectorInFile(filenameX, N, X);
}
void GaussSeidel(int N, double** A, double* B, double* X, double epsilon)
{
	//First approximation of x
	for (int i = 0; i < N; i++) X[i] = B[i] / A[i][i];
	//Calc
	double maxDelta;
	do
	{
		//if (maxDelta < fabs(XNew[i] - X[i])) maxDelta = fabs(XNew[i] - X[i]);
		maxDelta = -1.0;

		//New X's
		for (int i = 0; i < N; i++)
		{
			double tmp = X[i];

			double sum = 0.0;
			for (int j = 0; j < N; j++) if (i != j) sum += A[i][j] * X[j];
			X[i] = (B[i] - sum) / A[i][i];

			tmp = fabs(X[i] - tmp);
			if (maxDelta < tmp) maxDelta = tmp;
		}
	} while (maxDelta >= epsilon);
	//Output
	char filenameX[40]; sprintf_s(filenameX, 40, "VectorX%dGaussSeidel.txt", N);
	WriteVectorInFile(filenameX, N, X);
}
int main()
{
	system("chcp 1251>nul");
	srand(time(NULL));

	//Settings
	int N = 5;
	int maxElem = 10;
	double epsilon = 0.0001;

	//Filenames
	char filenameA[40]; sprintf_s(filenameA, 40, "MatrixA%d.txt", N);
	char filenameB[40]; sprintf_s(filenameB, 40, "VectorB%d.txt", N);

	//Generate input data
	GenerateMatrixFile(filenameA, N, maxElem);
	GenerateVectorFile(filenameB, N, maxElem);

	//Data arrays
	double** MatrixA = new double* [N]; for (int i = 0; i < N; i++) MatrixA[i] = new double[N];
	double* VectorB = new double[N];
	double* VectorX = new double[N];

	//Read input
	ReadMatrixFromFile(filenameA, N, MatrixA);
	ReadVectorFromFile(filenameB, N, VectorB);

	//Calculations & Output
	cout << "Started calculations..." << endl;
	Gauss(N, MatrixA, VectorB, VectorX);
	Jacobi(N, MatrixA, VectorB, VectorX, epsilon);
	GaussSeidel(N, MatrixA, VectorB, VectorX, epsilon);
	cout << "Complete!" << endl;

	//Free memory
	for (int i = 0; i < N; i++) delete MatrixA[i];
	delete[]MatrixA;
	delete[]VectorB;
	delete[]VectorX;

	system("pause>nul");
	return 0;
}