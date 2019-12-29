#include"linear-sys.h"
// MIT License
// Author: Iman Khademi, December 2019
// http://imankhademi.com
/* main function */

#define TEST_NUMBER 3


// Example 1: LQR Control
void runTest1()
{
	// reference: https://ece.gmu.edu/~gbeale/ece_620/discrete_lqr_01/discrete_lqr_01.html
	int n = 4;
	Matrix a({ {2.f} });
	Matrix b({ {1.0f} });
	Matrix c({ {1.0} });
	ColumnVector x0 = { 100.f };
	Matrix Q({ {10.f} });
	Matrix Qf({ {50.f} });
	Matrix R({ {1.0f} });
	LQR_Control myLQR(a, b, c, x0, Q, Qf, R, n);
	myLQR.runLQRControl();
	myLQR.exportResults();
	return;
}

// Example 2: LQR Control
Matrix identityMatrix(int size);
void runTest2()
{
	// reference: https://stanford.edu/class/ee363/lectures/dlqr.pdf
	int n = 20; // control horizon
	float rho = 1.0f;
	Matrix a({ {1.f,1.f},{0.f,1.f} });
	Matrix b({ {0.f},{1.f} });
	Matrix c({ {1.0f,0.f} });
	ColumnVector x0 = { {1.f,0.f} };
	Matrix Q({ {1.f,0.f},{0.f,0.f} });// C'*C
	Matrix Qf(Q);
	Matrix R;
	R = identityMatrix(1);
	R = rho * R;
	LQR_Control myLQR(a, b, c, x0, Q, Qf, R, n);
	myLQR.runLQRControl();
	myLQR.exportResults();
	return;
}

// Example3: Kalman Filter
ostream& operator<<(ostream& out, Vector2D v);
void runTest3()
{
	Matrix a({ {1.9223f,-0.9604f},{1.f,0.f} });
	Matrix b({ {1.0f},{1.0f} });
	Matrix c({ {1.0f,.0f} });
	Matrix g({ {1.0f},{1.0f} });
	Matrix P0({ {0.01f,0.},{0.f,0.02f} });// must be Symmetric Psitive Semi-Definite
	ColumnVector x0 = { {0.,0.} };
	Matrix wp({ {0.01} });
	Matrix vm({ {0.01} });
	Kalman_Filter filter(a, b, c, g, wp, vm, P0, x0);
	filter.runKalmanFilter(100, { {1.25} });
	superVector K = filter.getKalmanGain();
	filter.exportResults();// kalman_results.csv
}

int main()
{
	unsigned char test = TEST_NUMBER;
	switch (test)
	{
	case 1:
		runTest1();
		break;
	case 2:
		runTest2();
		break;
	case 3:
		runTest3();
		break;
	}
	return 0;
}