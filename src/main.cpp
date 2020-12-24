/**
* @file main.cpp
* @author Roney D da Silva
* @date 22 Dec 2020
* @copyright 2020 Roney D da Silva
         Email: roneyddasilva@gmail.com
* @brief O programa fornece dados compensados a partir dos metodos de
compensacao. A atitude eh obtida com o algoritmo denominado AQUA e um filtro de
Kalman linear.
*/

#include "../eigen/Eigen/Dense"
#include "AT.h"
#include "IMU.h"
#include "TELEMETRIA.h"
#include <Arduino.h>

// identificadores de tarefas
TaskHandle_t estimador;
TaskHandle_t comunicacao;

const float dtKalman = 0.025;    // Periodo do filtro de Kalman em segundos
const float dtTelemetria = 25.0; // Periodo para telemetria milissegundos
// instancia da biblioteca de algebra linear Eigen
using namespace Eigen;
// cria vetores para armazenar as medidas dos sensores no sistema do corpo
Vector3f acel = Vector3f::Zero(), gyro = Vector3f::Zero(),
         mag = Vector3f::Zero(), euler = Vector3f::Zero();
// vetor auxiliar para medidas serem transmitidas vias ESPNOW
VectorXf temp = VectorXf::Zero(4 + 3 + 3 + 3 + 3 + 1);
// cria objeto para chamar a imu
IMU IMU(Wire, 0x68);
// cria objeto para inicializar a comunicacao
TELEMETRIA TELEMETRIA;
// instancia a atitude
AT AT;
// funcao que roda dentro do loop do nucleo 0
void estimadorCodigo(void *) {
  // loop infinito no nucleo 0
  TickType_t xUltimoTempoKalman;
  xUltimoTempoKalman = xTaskGetTickCount();
  // xPeriodoKalman contem o perido em tick para computo de um tempo mais
  // preciso pdMS_TO_TICKS(dtKalman) converte as medias de tempo em millisegundo
  // para ticks
  const TickType_t xPeriodoKalman = pdMS_TO_TICKS(1000 * dtKalman);
  for (;;) {
    // le as medidas da IMU, compensa e atualiza os vetores no sistema do corpo
    IMU.readSensor();
    // com as medias atualizado o filtro de Kalman estima o estado
    AT.kalman(gyro, acel, mag, dtKalman);
    /* este metodo faz o loop parar ate que se passe xPeriodo considerando o
     * xUltimoTempoKalman sem consumir ciclos de processamento. Note que a
     * variavel xUltimoTempoKalman eh incrementada devido a referencia em
     * vTaskDelayUntil() quando*/
    vTaskDelayUntil(&xUltimoTempoKalman, xPeriodoKalman);
  }
}
// funcao que roda dentro do loop do nucleo 1
void comunicacaoCodigo(void *) {
  // mesmas informacoes que estao disponiveis no metodo estimadorCodigo
  TickType_t xUltimoTempoTelemetria;
  const TickType_t xPeriodoTelemetria = pdMS_TO_TICKS(dtTelemetria);
  xUltimoTempoTelemetria = xTaskGetTickCount();
  for (;;) {
    // normal
    temp << AT.qEst, AT.quat2euler(AT.qEst), acel, gyro, AT.m, dtKalman;
    TELEMETRIA.envia(temp);
    vTaskDelayUntil(&xUltimoTempoTelemetria, xPeriodoTelemetria);
  }
}
void setup() {
  pinMode(12, OUTPUT); // led branco
  pinMode(16, OUTPUT); // led interno comandos inversos
  digitalWrite(12, LOW);
  digitalWrite(16, HIGH);
  TELEMETRIA.begin();
  Serial.begin(115200);
  IMU.begin(acel, gyro, mag);
  TELEMETRIA.pisca(12);
  vTaskDelay(pdMS_TO_TICKS(3000)); // mesma coisa que delay(3000);
  TELEMETRIA.enviaMensagem((char *)"Coletando dado do Giro");
  /* calibracao do giro e do acelerometro eliminacao de bias */
  digitalWrite(12, LOW);
  // /* calibracao do giro por eliminacao de bias */
  digitalWrite(16, LOW);
  IMU.calibraGyro();
  TELEMETRIA.enviaMensagem((char *)"Fim da calibracao do Giro");
  digitalWrite(16, HIGH);
  TELEMETRIA.pisca(12);
  digitalWrite(16, LOW);
  digitalWrite(16, HIGH);
  TELEMETRIA.enviaMensagem((char *)"Cal. Mag");
  digitalWrite(12, LOW);
  TELEMETRIA.pisca(12);
  delay(2000);
  digitalWrite(12, HIGH);
  TELEMETRIA.enviaMensagem((char *)"Col. dados cal. Mag");
  // calibracao Magnetometro com o metodo geometrico
  while (IMU.calibracaoMagnetometro(22.8968) < 0) {
    TELEMETRIA.enviaMensagem((char *)"Falha na Cal. do Mag");
    delay(100);
    TELEMETRIA.enviaMensagem((char *)"Col. dados Novamente");
    delay(2000);
  }
  digitalWrite(12, LOW);
  TELEMETRIA.enviaMensagem((char *)"Cal. Mag Ok!");
  // Cria a tarefa para executar a funcao estimadoCodigo(), com prioridade 1 e
  // executado no nucleo 0
  xTaskCreatePinnedToCore(estimadorCodigo, /* Nome da funcao */
                          "estimador",     /* Nome da tarefa */
                          10000,           /* tamanho do Stack da tarefa */
                          NULL,            /* paramentros que a funcao recebe */
                          1,               /* prioridade da tarefa */
                          &estimador,      /* Identificador da tarefa */
                          0);              /* numero do nucleo */

  delay(500);
  // criacao da tarefa para telemetria
  xTaskCreatePinnedToCore(comunicacaoCodigo, "comunicacao", 10000, NULL, 1,
                          &comunicacao, 1);
  delay(500);
}
void loop() {}
