/**
 * @file METODOSCALIBRACAO.h
 * @author Roney D. da Silva
 * @date 8 Apr 2020
 * @copyright 2020 Roney Duarte
 * biblioteca contendo algoritmos de calibracao para acelerometro e
 * magnetometro. Para correto funcionamento a bibioteca EIGEN deve estar
 * indicada no include abaixo
 */
#ifndef METODOSCALIBRACAO_H_
#define METODOSCALIBRACAO_H_
// local da pasta com a biblioteca EIGEN
#include "../../eigen/Eigen/Dense"
using namespace Eigen;
// comente ou descomente a linha abaixo dependendo da plataforma utilizada
#include <Arduino.h>
class METODOSCALIBRACAO {
protected:
public:
  int calibracaoGeometrica(MatrixXf &data, Matrix3f &sF, Vector3f &bias,
                           float moduloCampo) {
    /**
* Metodo geometrico para calibracao de magnetometros
* @param [Matrix com dados] <data> <Matrix com N linha e 3 colunas contendo os
dados de amostra para a calibracao> [Matrix de scala] <sF> Matrix 3x3 com
diagonal contendo os fatores de escalas [bias vector] <bias> Vector 3x1 with the
bias - ellipsoid center [Magnetic field Earth] <moduloCampo> values in uT of
earth field
* @return [Matriz com fatores de escala na diagonal] <sF>
                                   [vetor com o centro da elipsoide] <bias>

* The method fitting an ellipsoid with least squares method. The compensated
measures  M_c are obtained with the mensures M proceeding as follows: M_c = sF *
(M - bias) The method needs a full rotations for a good performance
*/
    // data eh uma matrix com N x 3, sendo N o numero de medidas coletadas
    if (data.rows() == 3) {
      data.transposeInPlace();
    }
    sF = Matrix3f::Identity();
    bias = Vector3f::Zero();
    int N = data.rows();
    MatrixXf H = MatrixXf::Ones(N, 6);
    VectorXf w(N), X(6), tmp(3);
    H.col(0) = data.col(0);
    H.col(1) = data.col(1);
    H.col(2) = data.col(2);
    H.col(3) = -data.col(1).array().pow(2).matrix();
    H.col(4) = -data.col(2).array().pow(2).matrix();
    w = data.col(0).array().pow(2).matrix();
    data.resize(0, 0); // libera memoria
    // X = (H.transpose() * H).inverse() * H.transpose() * w;
    X = H.fullPivHouseholderQr().solve(w);
    H.resize(0, 0); // libera memoria
    bias(0) = X(0) / 2.0;
    bias(1) = X(1) / (2.0 * X(3));
    bias(2) = X(2) / (2.0 * X(4));
    tmp(0) = (X(5) + bias(0) * bias(0) + X(3) * bias(1) * bias(1) +
              X(4) * bias(2) * bias(2)) /
             (moduloCampo * moduloCampo);
    tmp(1) = tmp(0) / X(3);
    tmp(2) = tmp(0) / X(4);
    sF(0, 0) = 1.0 / sqrt(tmp(0));
    sF(1, 1) = 1.0 / sqrt(tmp(1));
    sF(2, 2) = 1.0 / sqrt(tmp(2));
    // checa se ha numero INFINITOS ou INVALIDOS retornoando -2
    // quando encontra numeros NaN e -1 para numeros insfinitos
    for (int i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        if (isnan(sF(i, j)) || isinf(sF(i, j))) {
          return -2;
        }
      }
      if (isnan(bias(i)) || isinf(bias(i))) {
        return -1;
      }
    }
    return 1;
  }
  int calibracaoAcelerometro(MatrixXf &X, MatrixXf data) {
    // data eh uma matrix com N x 3, sendo N o numero de medidas coletadas
    if (data.rows() == 3) {
      data.transposeInPlace();
    }
    int indice;
    int N = data.rows();
    MatrixXf Y = MatrixXf::Zero(N, 3);
    MatrixXf w = MatrixXf::Ones(N, 4);
    for (int var = 0; var < N; var++) {
      data.row(var).maxCoeff(&indice);
      if (data(var, indice) > 1000.) {
        Y(var, indice) = 9.81;
        w.row(var) << data.row(var), 1.0;
      } else {
        data.row(var).minCoeff(&indice);
        Y(var, indice) = -9.81;
        w.row(var) << data.row(var), 1.0;
      }
    }
    X = w.colPivHouseholderQr().solve(Y);
    X.transposeInPlace();
    for (int i = 0; i < X.rows(); ++i) {
      for (int j = 0; j < X.cols(); ++j) {
        if (isnan(X(i, j)) || isinf(X(i, j))) {
          return 0;
        }
      }
    }
    return 1;
  }
};

#endif /* METODOSCALIBRACAO_H_ */
