#ifndef AT_H_
#define AT_H_
#include "../../eigen/Eigen/Dense"
#include <iostream>
//#include </home/roney/Dropbox/bibliotecasCplus/eigen/Eigen/Geometry>
#include "Arduino.h"

// using namespace std;
using namespace Eigen;

class AT {
public:
	AT();
	void kalman(Vector3f giro, Vector3f acel, Vector3f mag, float deltaTempo);
	void observacao();
	void MagQuat();
	void covObs();
	void propagacaoKalman();
	void propagacao(Vector3f giro, float deltaTempo);
	void atualizaMatrixTransicao();
	void atualizaMatrixCovariancia();
	void MagQuat_azNegativo();

	Vector4f multQuat(const Vector4f p, const Vector4f q);
	Vector3f quat2euler(const Vector4f q);
	void rotVetQuat(const Vector4f p_A_B, Vector4f &v_A, const Vector4f v_B);
	// Matrix3f quat2Matrix(const Vector4f &q);
	// Matrix3f crossMatrix(const Vector4f &q);
	float dt = 0.0;
	Matrix4f omega = Matrix4f::Zero();
	const Matrix4f II = Matrix4f::Identity();
	const Matrix3f II3 = Matrix3f::Identity();
	Vector3f g, a, m, l, euler;
	MatrixXf gamaGiro = MatrixXf::Zero(4, 3);
	// quaternion sistema intermediario para local
	Vector4f qAcel = Vector4f::Zero(4);
	// quaternion sitema intermediario -> global
	Vector4f qMag = Vector4f::Zero(4);
	//	quaternion de observacao sistema global -> sistema local
	Vector4f qObs = Vector4f::Zero(4);
	// quaternion propagado atraves dos giroscopios. sistema global -> local
	Vector4f qProp { 1.0f, 0.0f, 0.0f, 0.0f };
	// quaternion global global -> local estimado
	Vector4f qEst { 1.0f, 0.0f, 0.0f, 0.0f };
	// matrizes do filtro de Kalman
	Matrix4f A, R, Q, K = Matrix4f::Zero(), P = Matrix4f::Ones();
	Matrix3f RotAcel;
	//	Matrix4f R = MatrixXf::Zero(6, 6);
	//	Matrix4f Q = MatrixXf::Zero(6, 6);
	//	MatrixXf covAcelMag = MatrixXf::Zero(6, 6);
	//	Matrix3f covGiro = Matrix3f::Zero(3, 3);
	const Matrix3f covGiro =  (Matrix3f() << 2.5e-02, 0.0, 0.0, 0.0, 2.5e-02, 0.0, 0.0, 0.0, 2.5e-02).finished();
	const MatrixXf covAcelMag =
			(MatrixXf(6, 6) << 3.e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0e-02, 0.0,
					        0.0, 0.0, 0.0, 0.0, 0.0, 3.0e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3e-02,
					        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3e-02)
					           .finished();

	float kJac = 0.0f, gamaJac = 0.0f, betaJac1 = 0.0f, betaJac2 = 0.0f;
	MatrixXf jacQuat = MatrixXf::Zero(4, 6);
private:
	// floats temporarios para facilitar computos
	Vector4f qTemp1 = Vector4f::Zero(4), qAcelConj = Vector4f::Zero(4), qTemp3 =
			Vector4f::Zero(4);
	const float craiz_2 = sqrt(2.0);
	// variaveis para calculos das covariancias
	MatrixXf df1_f2 = MatrixXf::Zero(8, 6);
	MatrixXf df2_u = MatrixXf::Identity(6, 6);
	MatrixXf dq_f1 = MatrixXf::Zero(4, 8);
	void dqf1();
	void dqf1_az_negativo();
	void df2u();
	void df2u_az_negativo();
	// para az nao negativo
	void df1f2();
	void df1f2_lx_negativo();
	// para az negativo
	void df1f2_az_negativo();
	void df1f2_az_negativo_lx_negativo();
	Matrix3f crossMatrix(const Vector4f q);
	Matrix3f quat2Matrix(const Vector4f q);
	void atualizaRotAcel();
};

#endif /* AT_H_ */
