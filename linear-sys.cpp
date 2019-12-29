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
vector<float> operator-(vector<float> v1, vector<float> v2);

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

void DT_SS_Model::updateState(ColumnVector u)
{
	vector<float> x = X.back();
	X.push_back(A * x + B * u);
	return;
}

void SS_Noisy_Model::updateState(ColumnVector u, ColumnVector w)
{
	vector<float> x = getStates().back();
	Matrix A, B;
	A = getA();
	B = getB();
	x = A * x + B * u + G * w;
	getStates().push_back(x);
	return;
}

void DT_SS_Model::outputMeasurement()
{
	vector<float> x = X.back();
	Y.push_back(C * x);
}

void SS_Noisy_Model::outputMeasurement(ColumnVector v)
{
	vector<float> x = getStates().back();
	Matrix C;
	C = getC();
	getOutput().push_back( C * x + v );
}

void SS_Noisy_Model::sysUpdate(ColumnVector ui)
{
	ColumnVector w = zeroMeanAWGNSample(Wc);
	ColumnVector v = zeroMeanAWGNSample(Vc);
	updateState(ui,w);
	outputMeasurement(v);
	return;
}

void DT_SS_Model::sysUpdate(ColumnVector ui)
{
	updateState(ui);
	outputMeasurement();
	return;
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
	for (int i = 0; i < getA().getSize(0); i++)
	{
		results << "x" << i << ",";
	}
	for (int i = 0; i < getC().getSize(0); i++)
	{
		results << "y" << i << ",";
	}
	results << endl;
	for (size_t i = 0; i < getStates().size(); i++)
	{
		results << i << ",";
		results << X[i];
		results << Y[i];
		results << endl;
	}
	return;
}

SS_Noisy_Model::SS_Noisy_Model(Matrix a, Matrix b, Matrix c,Matrix g, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0):DT_SS_Model{ a,b,c,x0 }
{
	Wc = wp;
	Vc = vm;
	P = P0;
	G = g;
	RowVector xv = zeroMeanAWGNSample(P);
	x0 = x0 + xv;// add uncertainity to the initial conditions
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

LQR_Control::LQR_Control(Matrix a, Matrix b, Matrix c, ColumnVector x0, Matrix qx, Matrix qf, Matrix ru, int horizon):DT_SS_Model{ a,b,c,x0 }
{
	Q = qx;
	Qf = qf;
	R = ru;
	Nh = horizon;
}

void LQR_Control::RicattiSolver()
{
	Pr.push_back(Qf.getData());
	Matrix P_t1; // P(t-1)
	Matrix A,ATr,B,BTr,L,APA,P_t,APB, BPA, R_BPBi;
	A = getA();
	B = getB();
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
		ColumnVector Xt = getStates().back();
		ColumnVector ut = Ki * Xt;
		sysUpdate(ut);
	}
	return;
}

void LQR_Control::exportResults(string filename)
{
	ofstream results(filename);
	results << "sample,";
	for (int i = 0; i < getA().getSize(0); i++)
	{
		results << "x" << i << ",";
	}
	for (int i = 0; i < getC().getSize(0); i++)
	{
		results << "y" << i << ",";
	}
	for (int i = 0; i < getB().getSize(1); i++)
	{
		for (int j = 0; j < getA().getSize(1); j++)
		results << "K" << i << j << ",";
	}
	results << endl;
	for (size_t i = 0; i < getStates().size(); i++)
	{
		results << i << ",";
		results << getStates()[i];
		results << getOutput()[i];
		if (i < getStates().size()-1)
			for (int j = 0; j < getB().getSize(1); j++)
			{
				results << FeedbackGain[i][j];
			}
		results << endl;
	}
	return;
}

//Kalman Filter for the model given by SS_Noisy_Model
// Reference: "Stochastic Optimal Linear Estimation and Control" by Meditch

Kalman_Filter::Kalman_Filter(Matrix a, Matrix b, Matrix c, Matrix g, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0) :SS_Noisy_Model(a, b, c, g, wp, vm, P0, x0)
{
	ErrorCovariance.push_back(P0.getData());
	ColumnVector x0_hat;
	x0_hat.assign(x0.size(),0); // initial estimate is 0
	Xhat.push_back(x0_hat);
}

void Kalman_Filter::singleStepEstimation(ColumnVector referenceInput, ColumnVector observation)
{
	Matrix Pkk;
	Matrix Pk1_k,Pk1_k1,A,B,C,G,APAt,GWGt,At,Bt,Ct,Gt,W,V,K,CPCt,PCt,KC,I;
	ColumnVector Xhat_kk = Xhat.back();// the latest estimation vector
	ColumnVector Innovation,CAx,Xhat_k1k1,Xhat_k1;
	A = getA();
	B = getB();
	C = getC();
	G = getG();
	W = getWcov();
	V = getVcov();
	At = A.Transpose();
	Ct = C.Transpose();
	Gt = G.Transpose();

	Pkk = ErrorCovariance.back(); // P(k|k)
	APAt = A * Pkk;
	APAt = APAt * At; // A*P(k|k)*A'

	GWGt = G * W;
	GWGt = GWGt * Gt; // G*W*G'

	Pk1_k = APAt + GWGt; // P(k+1|k):Estimation Error Covariance

	PCt = Pk1_k * Ct;
	CPCt = C * PCt;

	CPCt = CPCt + V;

	CPCt = Inverse(CPCt);

	K = PCt * CPCt;// Kalman Gain

	KalmanGain.push_back(K.getData());// attach it

	KC = K * C;

	I = identityMatrix(KC.getSize(0));

	KC = I - KC;

	Pk1_k1 = KC * Pk1_k; // P(k+1|k+1)
	ErrorCovariance.push_back(Pk1_k1.getData());// attach it

	Xhat_k1 = A * Xhat_kk + B * referenceInput;// page: 200, problem 5.8

	Innovation = observation - C * Xhat_k1;

	Xhat_k1k1 = Xhat_k1 + K * Innovation;// Xhat(k+1|k+1): The optimal filtered estimate
	Xhat.push_back(Xhat_k1k1);// attach it

	return;
}

void Kalman_Filter::runKalmanFilter(int samples, ColumnVector referenceInput)
{
	for (int i = 1; i <= samples; i++)
	{
		sysUpdate(referenceInput); // run themain process
		singleStepEstimation(referenceInput, getOutput().back()); // send the latest observation to the Kalman Filter
	}
}

void Kalman_Filter::exportResults(string filename)
{
	ofstream results(filename);
	results << "sample,";
	for (int i = 0; i < getA().getSize(0); i++)
	{
		results << "x" << i << ",";
	}
	for (int i = 0; i < getA().getSize(0); i++)
	{
		results << "x_hat" << i << ",";
	}
	for (int i = 0; i < getA().getSize(0); i++)
	{
		for (int j = 0; j < getC().getSize(0); j++)
			results << "K" << i << j << ",";
	}
	results << endl;
	for (size_t i = 0; i < getStates().size(); i++)
	{
		results << i << ",";
		results << getStates()[i];
		results << getEstimates()[i];
		if (i > 0)
			for (size_t j = 0; j < KalmanGain[0].size(); j++)
			{
				for(int k =0;k<KalmanGain[0][0].size();k++)
					results << KalmanGain [i-1][j][k] << ",";
			}
		results << endl;
	}
	return;
}
