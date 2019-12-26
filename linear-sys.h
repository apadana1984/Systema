#include <random>
#include"matrix.h"

// MIT LicenseF
// Author: Iman Khademi, December 2019
// http://imankhademi.com

typedef vector<Vector2D> superVector;

class DT_SS_Model // Discrete-Time Linear State-Space Model
{
public:
	Matrix A;
	Matrix B;
	Matrix C;
	Vector2D X;
	Vector2D Y;

	DT_SS_Model() {};
	DT_SS_Model(Matrix a, Matrix b, Matrix c, ColumnVector x0);
	Vector2D getStates();
	Vector2D getOutput();
	void sysUpdate(ColumnVector ui);
	void stepResponse(int noSamples);
	void exportSysVariables(string filename="results.csv");
};

class SS_Noisy_Model :public DT_SS_Model
{
public:
	SS_Noisy_Model(void) {};
	SS_Noisy_Model(Matrix a, Matrix b, Matrix c, Matrix wp, Matrix vm, Matrix P0, ColumnVector x0);
	Matrix Vc;// covariance matrix of measurement noise
	Matrix Wc;// covariance matrix of system noise
	Matrix P;// Covariance of initial state estimated value
	void sysUpdate(ColumnVector ui);
};

class LQR_Control
{
private:
	// Minimize: J = X(N)'*Qf*X(N) + SUM(X'*Q*X + U'*R*U)
	Matrix Q;// State Cost (Positive Semi-Definite)
	Matrix Qf;// Final State Cost (Positive Semi-Definite)
	Matrix R; // Iput Cost (Positive Definite)
	int Nh; // Finite Time Horizon
	DT_SS_Model system;
	superVector FeedbackGain;
	superVector Pr; // Ricatti Equation Solution Matrix over time
public:
	void RicattiSolver();
	superVector getFeedbackGain() { return FeedbackGain; };
	superVector getRicattiSolution() { return Pr; };
	LQR_Control(DT_SS_Model sys, Matrix qx, Matrix qf, Matrix ru, int horizon);
	void runLQRControl();
	void exportResults(string filename = "results.csv");
};