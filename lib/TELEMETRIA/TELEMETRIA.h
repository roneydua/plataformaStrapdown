#ifndef TELEMETRIA_h
#define TELEMETRIA_h
#include "../../eigen/Eigen/Dense"

#include "Arduino.h"
#include "WiFi.h"
extern "C" {
#include "esp_now.h"
}

using namespace Eigen;

class TELEMETRIA {

public:
  TELEMETRIA();
  void begin();
  void envia(VectorXf a);
  void enviaMensagem(char *aviso);
  uint8_t _broadcast[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  void pisca(uint8_t porta);
  typedef struct __attribute__((packed)) mensagem {
    float dados[4+3+3+3+3+1];
    char info[35]; // dados para informacao;
  } mensagem;

  typedef struct __attribute__((packed)) controle {
    uint8_t dados;
  } controle;
  static uint8_t comAceito;
  uint8_t _quantDados;

private:
  static void OnDataRecv(const uint8_t *_broadcast, const uint8_t *data,
                         int data_len);
};
#endif
