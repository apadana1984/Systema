#include"linear-sys.h"
// MIT License
// Author: Iman Khademi, December 2019
// http://imankhademi.com
/* main function */

#define TEST_NUMBER 2

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
	DT_SS_Model mySys(a, b, c, x0);
	LQR_Control myLQR(mySys, Q, Qf, R, n);
	myLQR.runLQRControl();
	myLQR.exportResults();
	return;
}

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
	DT_SS_Model mySys(a, b, c, x0);
	LQR_Control myLQR(mySys, Q, Qf, R, n);
	myLQR.runLQRControl();
	myLQR.exportResults();
	return;
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
	}
	return 0;
}