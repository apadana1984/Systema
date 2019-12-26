#include<fstream>
#include<random>
#include"linear-sys.h"

// MIT License
// Author: Iman Khademi, December 2019
// http://imankhademi.com

enum dtssError { INVALID_A,INVALID_B, INVALID_C, INVALID_D };

vector<float> operator+(vector<float> v1, vector<float> v2);
ostream& operator<<(ostream& out, RowVector v);
ostream& operator<<(ostream& out, Vector2D v);
Matrix Cholesky(Matrix M);
RowVector zeroMeanAWGNSample(Matrix P);
Matrix Inverse(Matrix m);
Vector2D Transpose(Vector2D A);
Matrix identityMatrix(int size);
ColumnVector operator*(Vector2D vl, ColumnVector vr);

DT_SS_Model::DT_SS_Model(Matrix a, Matrix b, Matrix c, ColumnVector x0)
{
	try 
	{
		if (!a.isSquare()) throw INVALID_A;
		if (a.getSize(0) != b.getSize(0)) throw INVALID_B;
		if (c.getSize(1) != a.getSize(0)) throw INVALID_C;
		A = a;
		B = b;
		C = c;
		X.push_back(x0);
		Y.push_back(C * x0);
	}
	catch (dtssError err)
	{
		switch (err)
		{
		case INVALID_A:
			break;
		case INVALID_B:
			break;
		case INVALID_C:
			break;
		}
	}
}

void DT_SS_Model::sysUpdate(ColumnVector ui)
{
	vector<float> x = X.back();
	X.push_back(A*x + B*ui);
	Y.push_back(C * x);
	return;
}

Vector2D DT_SS_Model::getStates()
{
	return X;
}

Vector2D DT_SS_Model::getOutput()
{
	return Y;
}

void DT_SS_Model::stepResponse(int noSamples)
{
	try
	{
		if (noSamples <= 0) throw 1;
	}
	catch (...)
	{
		ofstream logFile("user-log.txt", ios::app);
		logFile << "stepResponse: Number of Samples must be positive!" << endl;
		logFile.close();
		return;
	}
	for (int i = 0; i < noSamples; i++)
		sysUpdate({ 1.0f });
}

void DT_SS_Model::exportSysVariables(string filename)
{
	ofstream results(filename);
	results << "sample,";
	for (int i = 0; i < A.getSize(0); i++)
	{
		results << "x" << i << ",";
	}
	for (int i = 0; i < C.getSize(0); i++)
	{
		results << "y" << i << ",";
	}
	results << endl;
	for (size_t i = 0; i < X.size(); i++)
	{
		results << i << ",";
		results << X[i];
		results << Y[i];
		results << endl;
	}
	return;
}

SS_Noisy_Model::SS_Noisy_Model(Matrix a, Matrix b, Matrix c, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0)
{
	try
	{
		if (!a.isSquare()) throw INVALID_A;
		if (a.getSize(0) != b.getSize(0)) throw INVALID_B;
		if (c.getSize(1) != a.getSize(0)) throw INVALID_C;
		A = a;
		B = b;
		C = c;
		Wc = wp;
		Vc = vm;
		P = P0;
		RowVector xv = zeroMeanAWGNSample(P);
		x0 = x0 + xv;// add uncertainity to the initial conditions
		X.push_back(x0);
		Y.push_back(C * x0);
	}
	catch (dtssError err)
	{
		switch (err)
		{
		case INVALID_A:
			break;
		case INVALID_B:
			break;
		case INVALID_C:
			break;
		}
	}
}

void SS_Noisy_Model::sysUpdate(ColumnVector ui)
{
	vector<float> x = X.back();
	RowVector w = zeroMeanAWGNSample(Wc); // system noise
	RowVector v = zeroMeanAWGNSample(Vc); // measurement noise
	X.push_back(A * x + B * ui + w);
	Y.push_back(C * x + v);
	return;
}

//Qc: Covariance Matrix (Positive Semi-definite Matrix)
RowVector zeroMeanAWGNSample(Matrix Qc)
{
	int n = Qc.getSize(0);
	RowVector v(n);
	Matrix Cc;
	Cc = Cholesky(Qc);
	random_device generator;
	normal_distribution<float> distribution(0.f, 1.f);
	for (int i = 0; i < n; i++)
	{
		v[i] = distribution(generator);
	}
	v = Cc * v;
	return v;
}
/* ****************** LQR Control ******************* */

LQR_Control::LQR_Control(DT_SS_Model sys, Matrix qx, Matrix qf, Matrix ru, int horizon)
{
	Q = qx;
	Qf = qf;
	R = ru;
	system = sys;
	Nh = horizon;
}

void LQR_Control::RicattiSolver()
{
	Pr.push_back(Qf.getData());
	Matrix P_t1; // P(t-1)
	Matrix A,ATr,B,BTr,L,APA,P_t,APB, BPA, R_BPBi;
	A = system.A;
	B = system.B;
	ATr = A.Transpose();
	BTr = B.Transpose();

	// Start the Loop Here
	for (int t = Nh; t > 0; t--) // Dynamic Programming goes backward in time
	{
		P_t = Pr.back();// P(t)
		APA = ATr * P_t;
		APA = APA * A; // A'*P*A
		P_t1 = Q + APA;

		APB = ATr * P_t;
		APB = APB * B;// A'*P*B

		BPA = APB.Transpose();// B'*P*A

		R_BPBi = BTr * P_t;
		R_BPBi = R_BPBi * B;
		R_BPBi = R_BPBi + R;
		R_BPBi = Inverse(R_BPBi);// (R+B'*P*B)^-1

		R_BPBi = R_BPBi * BPA; 
		R_BPBi = -1.f * R_BPBi;// Kt
		FeedbackGain.push_back(R_BPBi.getData());
		R_BPBi = APB * R_BPBi ;

		P_t1 = P_t1 + R_BPBi;// Solution
		Pr.push_back(P_t1.getData()); // attach it to the super vector
	}
	// make them forward in time
	reverse(Pr.begin(),Pr.end());
	reverse(FeedbackGain.begin(), FeedbackGain.end());
	return;
}

void LQR_Control::runLQRControl()
{
	RicattiSolver();// Solve the DP
	for (int i = 0; i < Nh; i++)
	{
		Vector2D Ki = FeedbackGain[i];
		ColumnVector Xt = system.X.back();
		ColumnVector ut = Ki * Xt;
		system.sysUpdate(ut);
	}
	return;
}

void LQR_Control::exportResults(string filename)
{
	ofstream results(filename);
	results << "sample,";
	for (size_t i = 0; i < system.A.getSize(0); i++)
	{
		results << "x" << i << ",";
	}
	for (size_t i = 0; i < system.C.getSize(0); i++)
	{
		results << "y" << i << ",";
	}
	for (size_t i = 0; i < system.B.getSize(1); i++)
	{
		for (size_t j = 0; j < system.A.getSize(1); j++)
		results << "K" << i << j << ",";
	}
	results << endl;
	for (size_t i = 0; i < system.X.size(); i++)
	{
		results << i << ",";
		results << system.X[i];
		results << system.Y[i];
		if (i < system.X.size()-1)
			for (size_t j = 0; j < system.B.getSize(1); j++)
			{
				results << FeedbackGain[i][j];
			}
		results << endl;
	}
	return;
}