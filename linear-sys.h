#include <random>
#include"matrix.h"

// MIT LicenseF
// Author: Iman Khademi, December 2019
// http://imankhademi.com

typedef vector<Vector2D> superVector;

/*
x(k+1) = A*x(k) + B*u(k)
y(k) = C*x(k)
*/

class DT_SS_Model // Discrete-Time Linear State-Space Model
{
private:
	Matrix A;
	Matrix B;
	Matrix C;
	Vector2D X;
	Vector2D Y;
public:
	DT_SS_Model() {};
	DT_SS_Model(Matrix a, Matrix b, Matrix c, ColumnVector x0);
	Matrix getA() { return A; };
	Matrix getB() { return B; };
	Matrix getC() { return C; };
	Vector2D& getStates() { return X; };
	Vector2D& getOutput() { return Y; };
	void updateState(ColumnVector u);
	void outputMeasurement();
	void sysUpdate(ColumnVector ui);
	void stepResponse(int noSamples);
	void exportSysVariables(string filename="results.csv");
};

/* 
x(k+1) = A*x(k) + B*u(k) + G*w(k)
y(k) = C*x(k) + v(k)

w : process noise, Wc = Cov(w)
v : measurement noise, Vc = Cov(v)
P0: Covariance of Estimation Error at k=0 : P0 = Cov(x(0) - X_hat(0))
x_hat : Estimation of 'x'
*/

class SS_Noisy_Model :public DT_SS_Model
{
private:
	Matrix Vc;// covariance matrix of measurement noise
	Matrix Wc;// covariance matrix of system noise
	Matrix P;// Covariance of estimation error
	Matrix G;
public:
	SS_Noisy_Model(void) {};
	SS_Noisy_Model(Matrix a, Matrix b, Matrix c, Matrix g, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0);
	Matrix getP() { return P; };
	Matrix getVcov() { return Vc; };
	Matrix getWcov() { return Wc; };
	Matrix getG() { return G; };
	void updateState(ColumnVector u, ColumnVector w);
	void outputMeasurement(ColumnVector v);
	void sysUpdate(ColumnVector ui);
};

class LQR_Control:public DT_SS_Model
{
private:
	// Minimize: J = X(N)'*Qf*X(N) + SUM(X'*Q*X + U'*R*U)
	Matrix Q;// State Cost (Positive Semi-Definite)
	Matrix Qf;// Final State Cost (Positive Semi-Definite)
	Matrix R; // Input Cost (Positive Definite)
	int Nh; // Finite Time Horizon
	superVector FeedbackGain;
	superVector Pr; // Ricatti Equation Solution Matrix over time
public:
	void RicattiSolver();
	superVector getFeedbackGain() { return FeedbackGain; };
	superVector getRicattiSolution() { return Pr; };
	LQR_Control(Matrix a, Matrix b, Matrix c, ColumnVector x0, Matrix qx, Matrix qf, Matrix ru, int horizon);
	void runLQRControl();
	void exportResults(string filename = "results.csv");
};

class Kalman_Filter:SS_Noisy_Model
{
private:
	superVector KalmanGain;
	superVector ErrorCovariance;
	Vector2D Xhat;
public:
	Kalman_Filter(Matrix a, Matrix b, Matrix c, Matrix g, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0);
	superVector getKalmanGain() { return KalmanGain; };
	Vector2D getEstimates() { return Xhat; };
	superVector getErrorCovariance() { return ErrorCovariance; };
	void singleStepEstimation(ColumnVector referenceInput,ColumnVector observation);
	void exportResults(string filename = "kalman_results.csv");
	void runKalmanFilter(int samples,ColumnVector uref);
};