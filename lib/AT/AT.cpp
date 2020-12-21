/*
 * Esta classe normaliza o vetor apenas no fim de cada estimacao
 */

#include "AT.h"

AT::AT() {
}
void AT::kalman(Vector3f giro, Vector3f acel, Vector3f mag, float deltaTempo) {
	a << acel.normalized();
	g << giro;
	m << mag.normalized();
	dt = deltaTempo;
	// propaga o modelo
	propagacaoKalman();
	// matriz de covariancia predita
	P = A * P * A.transpose() + Q;
	// computa o quaternio de observacao
	observacao();
	// calcula o ganho de Kalman
	K = P * (P.transpose() + R).inverse();
	// calculo do quaternion estimado
	qEst = qProp + K * (qObs - qProp);
	qEst.normalize();
	// calculo da matriz de covariancia
	P = (II - K) * P;

}
void AT::propagacao(Vector3f giro, float deltaTempo) {
	/* Integracao pura do quaternion de que depende apenas
	 * dos giroscopios */
	g << giro;
	dt = deltaTempo;
	atualizaMatrixTransicao();
	// propaga o modelo
	qProp = A * qProp;
	qProp.normalize();
}

void AT::propagacaoKalman() {
	/* Propaga o modelo baseado nos estados estimados no tempo anterior
	 e nas medidas atuais Atualiza a matrix de covariancia Q do modelo
	 */
	atualizaMatrixCovariancia();
	atualizaMatrixTransicao();
	// propaga o modelo
	qProp = A * qEst;
	qProp.normalize();
}

Vector3f AT::quat2euler(const Vector4f q) {
	/* funcao que extrai os angulos de Euler do
	 * quaternion fornecido */
	// calcula fi
	euler(0) = atan2(2 * (q(0) * q(1) + q(2) * q(3)),
			1.0f - 2.0f * (q(1) * q(1) + q(2) * q(2)));
	// calcula theta
	euler(1) = asin(2 * (q(0) * q(2) - q(3) * q(1)));
	// calcula psi
	euler(2) = atan2(2 * (q(0) * q(3) + q(1) * q(2)),
			1.0f - 2.0f * (q(2) * q(2) + q(3) * q(3)));
	euler *= -180.0f / M_PI;
	return euler;
}

void AT::atualizaMatrixCovariancia() {
	gamaGiro << qEst(1), qEst(2), qEst(3), -qEst(0), -qEst(3), -qEst(2), qEst(
			2), -qEst(0), -qEst(1), -qEst(2), qEst(1), -qEst(0);

	Q = 0.25 * dt * dt * gamaGiro * covGiro * gamaGiro.transpose();

}
void AT::atualizaMatrixTransicao() {
	omega(0, 0) = 0.0f;
	omega(0, 1) = g(0);
	omega(0, 2) = g(1);
	omega(0, 3) = g(2);

	omega(1, 0) = -g(0);
	omega(1, 1) = 0.0f;
	omega(1, 2) = -g(2);
	omega(1, 3) = g(1);

	omega(2, 0) = -g(1);
	omega(2, 1) = g(2);
	omega(2, 2) = 0.0f;
	omega(2, 3) = -g(0);

	omega(3, 0) = -g(2);
	omega(3, 1) = -g(1);
	omega(3, 2) = g(0);
	omega(3, 3) = 0.0f;
	A << II + dt * 0.5 * omega;
}

void AT::observacao() {
	/* atualiza o quaternion de observacao fornecendo-o
	 * normalizado
	 */
	/*
	 * calcula quaternion de rotação que leva um vetor
	 * do sistema intermediario para o sistema do corpo
	 * */
	if (a(2) > 0.0) {
		kJac = sqrt(a(2) + 1.0);
		qAcel(0) = kJac;
		qAcel(1) = -a(1) / kJac;
		qAcel(2) = a(0) / kJac;
		qAcel(3) = 0.0;
		qAcel.normalize();
//		atualizaRotAcel();
//		l = RotAcel*m;
		l = quat2Matrix(qAcel).transpose() * m;
		MagQuat();
	} else {
		kJac = sqrt(1.0 - a(2));
		qAcel(0) = -a(1) / kJac;
		qAcel(1) = kJac;
		qAcel(2) = 0.0;
		qAcel(3) = a(0) / kJac;
		// produto de todas as constantes
		qAcel.normalize();
//		atualizaRotAcel();
//		l = RotAcel*m;
		l = quat2Matrix(qAcel).transpose() * m;
		MagQuat_azNegativo();
	}
	qObs = multQuat(qMag, qAcel);
	qObs.normalize();
	if (qObs.dot(qEst) < 0.0) {
		qObs *= -1.0;
	}
	covObs();
}

void AT::covObs() {
// atualiza a matriz jacobiana para computo das covariancias de observacao
	jacQuat = dq_f1 * df1_f2 * df2_u;
	R = jacQuat * covAcelMag * jacQuat.transpose();
}

void AT::MagQuat() { // atualiza o quaternion
/*leva as medidas do magnetometro escritas no sistema local
 * para o sistema intermediario
 */
	gamaJac = l(0) * l(0) + l(1) * l(1);
	if (l(0) <= 0.0) {
		betaJac2 = sqrt(gamaJac - l(0) * sqrt(gamaJac));
		qMag(0) = l(1) / betaJac2;
		qMag(1) = 0.0;
		qMag(2) = 0.0;
		qMag(3) = betaJac2 / sqrt(gamaJac);
		qMag.normalize();
		dqf1();
		df1f2_lx_negativo();
		df2u();
	} else {
		betaJac1 = sqrt(gamaJac + l(0) * sqrt(gamaJac));
		qMag(0) = betaJac1 / sqrt(gamaJac);
		qMag(1) = 0.0;
		qMag(2) = 0.0;
		qMag(3) = l(1) / betaJac1;
		qMag.normalize();
		dqf1();
		df1f2();
		df2u();
	}
}
void AT::MagQuat_azNegativo() { // atualiza o quaternion
/*leva as medidas do magnetometro escritas no sistema local
 * para o sistema intermediario
 */
	gamaJac = l(0) * l(0) + l(1) * l(1);
	if (l(0) <= 0.0) {
		betaJac2 = sqrt(gamaJac - l(0) * sqrt(gamaJac));
		qMag(0) = l(1) / betaJac2;
		qMag(1) = 0.0;
		qMag(2) = 0.0;
		qMag(3) = betaJac2 / sqrt(gamaJac);
		qMag.normalize();
		dqf1_az_negativo();
		df1f2_az_negativo_lx_negativo();
		df2u_az_negativo();
	} else {
		betaJac1 = sqrt(gamaJac + l(0) * sqrt(gamaJac));
		qMag(0) = betaJac1 / sqrt(gamaJac);
		qMag(1) = 0.0;
		qMag(2) = 0.0;
		qMag(3) = l(1) / betaJac1;
		qMag.normalize();
		dqf1_az_negativo();
		df1f2_az_negativo();
		df2u_az_negativo();
	}
}
Vector4f AT::multQuat(const Vector4f p, const Vector4f q) {
	return Vector4f(p(0) * q(0) - p(1) * q(1) - p(2) * q(2) - p(3) * q(3),
			p(0) * q(1) + p(1) * q(0) + p(2) * q(3) - p(3) * q(2),
			p(0) * q(2) - p(1) * q(3) + p(2) * q(0) + p(3) * q(1),
			p(0) * q(3) + p(1) * q(2) - p(2) * q(1) + p(3) * q(0));
}
void AT::rotVetQuat(const Vector4f p_A_B, Vector4f &v_B, const Vector4f v_A) {
	/* funcao que leva um vetor escrito no sistema A para um sistema B
	 * dado um quaternion que leva de A para B
	 * */
	v_B =
			multQuat(multQuat(p_A_B, (Vector4f(4) << 0.0f, v_A).finished()),
					(Vector4f(4) << p_A_B(0), -p_A_B.block(1, 0, 3, 1)).finished()).block(
					1, 0, 3, 1);
}

void AT::dqf1() {
	/*funcao responsavel por atualizar a matrix dq/df1
	 * quando a componente az nao eh negativa
	 */
	dq_f1 = MatrixXf::Zero(4, 8);
	for (int var = 0; var < 4; ++var) {
		dq_f1(var, var) = qMag(0);
		dq_f1(3 - var, var) = pow(-1, var) * qMag(3);
		dq_f1(var, var + 4) = qAcel(0);
	}
	dq_f1(0, 5) = -qAcel(1);
	dq_f1(1, 4) = qAcel(1);
	dq_f1(2, 7) = -qAcel(1);
	dq_f1(3, 6) = qAcel(1);

	dq_f1(0, 6) = qAcel(2);
	dq_f1(1, 7) = qAcel(2);
	dq_f1(2, 4) = qAcel(2);
	dq_f1(3, 5) = -qAcel(2);
}

void AT::df1f2() {
	df1_f2 = MatrixXf::Zero(8, 6);
	df1_f2(0, 2) = 1.0f / kJac;
	df1_f2(1, 1) = -2.0 * df1_f2(0, 2);
	df1_f2(1, 2) = a(1) / pow(kJac, 3);
	df1_f2(2, 0) = -df1_f2(1, 1);
	df1_f2(2, 2) = -a(0) / pow(kJac, 3);

	df1_f2(4, 3) = pow(l(1), 2) / (betaJac1 * gamaJac);
	df1_f2(4, 4) = l(0) * l(1) / (betaJac1 * gamaJac);
	df1_f2(7, 3) = -l(1) * betaJac1 / pow(gamaJac, 1.5f);
	df1_f2(7, 4) = -l(0) * betaJac1 / pow(gamaJac, 1.5f);
	df1_f2 /= 2 * craiz_2;
}

void AT::df1f2_lx_negativo() {
	df1_f2 = MatrixXf::Zero(8, 6);
	df1_f2(0, 2) = 1.0f / kJac;
	df1_f2(1, 1) = -2.0 * df1_f2(0, 2);
	df1_f2(1, 2) = a(1) / pow(kJac, 3);
	df1_f2(2, 0) = -df1_f2(1, 1);
	df1_f2(2, 2) = -a(0) / pow(kJac, 3);

	df1_f2(4, 3) = l(1) * betaJac2 / pow(gamaJac, 1.5f);
	df1_f2(4, 4) = l(0) * betaJac2 / pow(gamaJac, 1.5f);
	df1_f2(7, 3) = -pow(l(1), 2) / (betaJac2 * gamaJac);
	df1_f2(4, 4) = l(0) * l(1) / (betaJac2 * gamaJac);
	df1_f2 /= 2 * craiz_2;
}

void AT::df2u() {
	/* computo de df2/du para compoente az nao negativa*/
	df2_u(3, 0) = m(2) - (2 * a(0) * m(0) + a(1) * m(1)) / pow(kJac, 2);
	df2_u(3, 1) = -a(0) * m(1) / pow(kJac, 2);
	df2_u(3, 2) = a(0) * (a(0) * m(0) + a(1) * m(1)) / pow(kJac, 4);
	df2_u(3, 3) = 1.0f - pow(a(0) / kJac, 2);
	df2_u(3, 4) = -a(0) * a(1) / pow(kJac, 2);
	df2_u(3, 5) = a(0);
	df2_u(4, 0) = -a(1) * m(0) / pow(kJac, 2);
	df2_u(4, 1) = m(2) - (a(0) * m(0) + 2.0f * a(1) * m(1)) / pow(kJac, 2);
	df2_u(4, 2) = a(1) * (a(0) * m(0) + a(1) * m(1)) / pow(kJac, 4);
	df2_u(4, 3) = -a(0) * a(1) / pow(kJac, 2);
	df2_u(4, 4) = 1.0f - pow(a(1) / kJac, 2);
	df2_u(4, 5) = a(1);
	df2_u(5, 0) = -m(0);
	df2_u(5, 1) = -m(1);
	df2_u(5, 2) = m(2);
	df2_u(5, 3) = -a(0);
	df2_u(5, 4) = -a(1);
	df2_u(5, 5) = a(2);
}

void AT::dqf1_az_negativo() {
	/*funcao responsavel por atualizar a matrix dq/df1
	 * quando a componente az eh negativa
	 */
	dq_f1 = MatrixXf::Zero(4, 8);
	for (int var = 0; var < 4; ++var) {
		dq_f1(var, var) = qMag(0);
		dq_f1(3 - var, var) = pow(-1, var) * qMag(3);
		dq_f1(var, var + 4) = qAcel(0);
	}
	dq_f1(0, 5) = -qAcel(1);
	dq_f1(1, 4) = qAcel(1);
	dq_f1(2, 7) = -qAcel(1);
	dq_f1(3, 6) = qAcel(1);

	dq_f1(0, 7) = -qAcel(3);
	dq_f1(1, 6) = -qAcel(3);
	dq_f1(2, 5) = qAcel(3);
	dq_f1(3, 4) = qAcel(3);

}

void AT::df1f2_az_negativo() {
	df1_f2 = MatrixXf::Zero(8, 6);
	df1_f2(0, 1) = -2.0 / kJac;
	df1_f2(0, 2) = -a(1) / pow(kJac, 3);
	df1_f2(1, 2) = df1_f2(0, 1) / 2.0;
	df1_f2(3, 0) = -df1_f2(0, 1);
	df1_f2(3, 2) = a(0) / pow(kJac, 3);

	df1_f2(4, 3) = pow(l(1), 2) / (betaJac1 * gamaJac);
	df1_f2(4, 4) = -l(0) * l(1) / (betaJac1 * gamaJac);
	df1_f2(7, 3) = -l(1) * betaJac1 / pow(gamaJac, 1.5);
	df1_f2(7, 4) = l(0) * betaJac1 / pow(gamaJac, 1.5);
	df1_f2 /= 2 * craiz_2;
}

void AT::df1f2_az_negativo_lx_negativo() {
	df1_f2 = MatrixXf::Zero(8, 6);
	df1_f2(0, 1) = -2.0 / kJac;
	df1_f2(0, 2) = -a(1) / pow(kJac, 3);
	df1_f2(1, 2) = df1_f2(0, 1) / 2.0;
	df1_f2(3, 0) = -df1_f2(0, 1);
	df1_f2(3, 2) = a(0) / pow(kJac, 3);

	df1_f2(4, 3) = l(1) * betaJac2 / pow(gamaJac, 1.5);
	df1_f2(4, 4) = -l(0) * betaJac2 / pow(gamaJac, 1.5);
	df1_f2(7, 3) = -pow(l(1), 2) / (betaJac2 * gamaJac);
	df1_f2(4, 4) = l(0) * l(1) / (betaJac2 * gamaJac);
	df1_f2 /= 2 * craiz_2;
}

void AT::df2u_az_negativo() {
	/* computo de df2/du para compoente az negativa*/
	df2_u(3, 0) = m(2) - (2 * a(0) * m(0) - a(1) * m(1)) / pow(kJac, 2);
	df2_u(3, 1) = a(0) * m(1) / pow(kJac, 2);
	df2_u(3, 2) = a(0) * (-a(0) * m(0) + a(1) * m(1)) / pow(kJac, 4);
	df2_u(3, 3) = 1.0 - pow(a(0) / kJac, 2);
	df2_u(3, 4) = -a(0) * a(1) / pow(kJac, 2);
	df2_u(3, 5) = a(0);
	df2_u(4, 0) = -a(1) * m(0) / pow(kJac, 2);
	df2_u(4, 1) = m(2) - (a(0) * m(0) - 2.0f * a(1) * m(1)) / pow(kJac, 2);
	df2_u(4, 2) = a(1) * (-a(0) * m(0) + a(1) * m(1)) / pow(kJac, 4);
	df2_u(4, 3) = -a(0) * a(1) / pow(kJac, 2);
	df2_u(4, 4) = -1.0f + pow(a(1) / kJac, 2);
	df2_u(4, 5) = a(1);
	df2_u(5, 0) = m(0);
	df2_u(5, 1) = m(1);
	df2_u(5, 2) = -m(2);
	df2_u(5, 3) = a(0);
	df2_u(5, 4) = -a(1);
	df2_u(5, 5) = -a(2);
}
Matrix3f AT::quat2Matrix(const Vector4f q) {
	return ((2.0 * q(0) * q(0) - 1.0) * II3 + 2.0 * q(0) * crossMatrix(q)
			+ 2.0 * q.block(1, 0, 3, 1) * (q.block(1, 0, 3, 1).transpose()));

}
void AT::atualizaRotAcel() {
	// computa a matrix transposta !!!!! apenas para melhor eficiencia
	RotAcel(0, 0) = qAcel(0) * qAcel(0) + qAcel(1) * qAcel(1)
			- qAcel(2) * qAcel(2) - qAcel(3) * qAcel(3);
	RotAcel(1, 0) = 2 * (qAcel(1) * qAcel(2) - qAcel(0) * qAcel(3));
	RotAcel(2, 0) = 2 * (qAcel(1) * qAcel(3) + qAcel(0) * qAcel(2));

	RotAcel(0, 1) = 2 * (qAcel(1) * qAcel(2) + qAcel(0) * qAcel(3));
	RotAcel(1, 1) = qAcel(0) * qAcel(0) - qAcel(1) * qAcel(1)
			+ qAcel(2) * qAcel(2) - qAcel(3) * qAcel(3);
	RotAcel(2, 1) = 2 * (qAcel(2) * qAcel(3) - qAcel(0) * qAcel(1));

	RotAcel(0, 2) = 2 * (qAcel(1) * qAcel(3) - qAcel(0) * qAcel(2));
	RotAcel(1, 2) = 2 * (qAcel(2) * qAcel(3) + qAcel(0) * qAcel(1));
	RotAcel(2, 2) = qAcel(0) * qAcel(0) - qAcel(1) * qAcel(1)
			- qAcel(2) * qAcel(2) + qAcel(3) * qAcel(3);
}
Matrix3f AT::crossMatrix(const Vector4f q) {
	return ((Matrix3f(3, 3) << 0.0, -q(3), q(2), q(3), 0.0, -q(1), -q(2), q(1), 0.0).finished());
}
