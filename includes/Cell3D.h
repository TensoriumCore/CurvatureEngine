#pragma once
#include <algorithm>
class Cell3D
{
	public:
		double gamma[3][3]         {};
		double gamma_inv[3][3]     {};
		double K[3][3]             {};
		double Christoffel[3][3][3]{};
		double alpha = 1.0;
		double beta[3] = {0.0, 0.0, 0.0};
		double rho   = 0.0;
		double p     = 0.0;
		double vx    = 0.0;
		double vy    = 0.0;
		double vz    = 0.0;
		double hamiltonian = 0.0;
		double momentum[3] = {0.0, 0.0, 0.0};
		double gamma0[3][3] = {};
		double K0[3][3]     = {};
		double alpha0       = 1.0;
		double beta0[3]     = {0.0, 0.0, 0.0};
		double T[3][3]      = {};
		double dgt[3][3]         = {};
		double dKt[3][3]         = {};
		double alphaStage[4]     = {};
		double betaStage[4][3]   = {};
		double gammaStage[4][3][3] = {};
		double KStage[4][3][3]     = {};
		double beta_con[3] = {0.0, 0.0, 0.0};
		double Gamma3[3][3][3] = {};
		double Ricci[3][3] = {};
		Cell3D() 
			: alpha(1.0), rho(0.0), p(0.0), vx(0.0), vy(0.0), vz(0.0), 
			hamiltonian(0.0), alpha0(1.0)
	{
		std::fill_n(beta, 3, 0.0);
		std::fill_n(momentum, 3, 0.0);
		std::fill_n(beta0, 3, 0.0);
		std::fill_n(beta_con, 3, 0.0);

		std::fill_n(&gamma[0][0], 9, 0.0);
		std::fill_n(&gamma_inv[0][0], 9, 0.0);
		std::fill_n(&K[0][0], 9, 0.0);
		std::fill_n(&Christoffel[0][0][0], 27, 0.0);
		std::fill_n(&gamma0[0][0], 9, 0.0);
		std::fill_n(&K0[0][0], 9, 0.0);
		std::fill_n(&T[0][0], 9, 0.0);
		std::fill_n(&dgt[0][0], 9, 0.0);
		std::fill_n(&dKt[0][0], 9, 0.0);
		std::fill_n(&Gamma3[0][0][0], 27, 0.0);
		std::fill_n(&Ricci[0][0], 9, 0.0);

		std::fill_n(alphaStage, 4, 0.0);
		std::fill_n(&betaStage[0][0], 12, 0.0);
		std::fill_n(&gammaStage[0][0][0], 36, 0.0);
		std::fill_n(&KStage[0][0][0], 36, 0.0);
	}
		inline void storeInitialState()
		{
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					gamma0[a][b] = gamma[a][b];
					K0[a][b]     = K[a][b];
				}
			}
			alpha0 = alpha;
			for (int m = 0; m < 3; m++) {
				beta0[m] = beta[m];
			}
		}

		inline void storeStage(int stage,
				double dtGamma[3][3],
				double dtK[3][3],
				double d_alpha_dt,
				const double d_beta_dt[3])
		{
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					gammaStage[stage][a][b] = dtGamma[a][b];
					KStage[stage][a][b]     = dtK[a][b];
				}
			}
			alphaStage[stage] = d_alpha_dt;
			for (int m = 0; m < 3; m++) {
				betaStage[stage][m] = d_beta_dt[m];
			}
		}

		inline void updateIntermediate(int stage, double dtCoeff)
		{
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					gamma[a][b] = gamma0[a][b] + dtCoeff * gammaStage[stage][a][b];
					K[a][b]     = K0[a][b]     + dtCoeff * KStage[stage][a][b];
				}
			}
			alpha = alpha0 + dtCoeff * alphaStage[stage];
			for (int m = 0; m < 3; m++) {
				beta[m] = beta0[m] + dtCoeff * betaStage[stage][m];
			}
		}

		inline void combineRK4(double dt)
		{
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					gamma[a][b] = gamma0[a][b]
						+ (dt / 6.0) * ( gammaStage[0][a][b]
								+ 2.0 * gammaStage[1][a][b]
								+ 2.0 * gammaStage[2][a][b]
								+       gammaStage[3][a][b] );
					K[a][b] = K0[a][b]
						+ (dt / 6.0) * ( KStage[0][a][b]
								+ 2.0 * KStage[1][a][b]
								+ 2.0 * KStage[2][a][b]
								+       KStage[3][a][b] );
				}
			}
			alpha = alpha0
				+ (dt / 6.0) * ( alphaStage[0]
						+ 2.0 * alphaStage[1]
						+ 2.0 * alphaStage[2]
						+       alphaStage[3] );
			for (int m = 0; m < 3; m++) {
				beta[m] = beta0[m]
					+ (dt / 6.0) * ( betaStage[0][m]
							+ 2.0 * betaStage[1][m]
							+ 2.0 * betaStage[2][m]
							+       betaStage[3][m] );
			}
		}
};

