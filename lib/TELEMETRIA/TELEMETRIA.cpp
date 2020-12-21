#include "TELEMETRIA.h"
#include "Arduino.h"

// cria um objeto que recebe como entrada o EspNow
TELEMETRIA::TELEMETRIA() { }
void TELEMETRIA::begin() {
  // disconeta WIfi para evitar instabilidades
  WiFi.disconnect();
  // inicializa wifi em modo station
  WiFi.mode(WIFI_STA);
  // tenta estabelecer coneccao
  while (esp_now_init() != ESP_OK) {
    ESP.restart();
  };
  Serial.println("ESP Inicializado");
  // Criamos uma variável que irá guardar as informações do slave
  esp_now_peer_info_t slave;
  slave.channel = 1; // Informamos o canal
  slave.encrypt = 0; // 0 para não usar criptografia ou 1 para usar
  // Copia o endereço do array para a estrutura
  memcpy(slave.peer_addr, _broadcast, sizeof(_broadcast));

  esp_now_add_peer(&slave); // Adiciona o slave
  esp_now_register_recv_cb(OnDataRecv);
}

void TELEMETRIA::envia(VectorXf _a) {
  // Array que irá armazenar os valores lidos
  mensagem msg;
  // float temp[_a.size()];
  strcpy(msg.info, "");
  for (size_t i = 0; i < _a.size(); i++) {
    msg.dados[i] = (float)_a(i);
    //Serial.print(msg.dados[i]);
  }
  //Serial.println();
  uint8_t bs[sizeof(mensagem)];

  memcpy(&bs, &msg, sizeof(msg));
  esp_now_send(_broadcast, bs, sizeof(bs));
}
void TELEMETRIA::enviaMensagem(char *aviso) {

  mensagem msg;
  for (size_t i = 0; i < 10; i++) {
    msg.dados[i] = 0.0f;
  }
  strcpy(msg.info, aviso);
  //Serial.println(msg.info);
  uint8_t bs[sizeof(mensagem)];
  memcpy(&bs, &msg, sizeof(msg));
  esp_now_send(_broadcast, bs, sizeof(bs));
}
// funcao responsavel por receber dados
void TELEMETRIA::OnDataRecv(const uint8_t *_broadcast, const uint8_t *data,
                            int data_len) {
  TELEMETRIA TELEMETRIA;
  // Para cada pino
  controle ctrl;
  memcpy(&ctrl, data, sizeof(controle));
  // atualizacontrole(ctrl.dados);
  TELEMETRIA::comAceito = ctrl.dados;
}
uint8_t TELEMETRIA::comAceito;
void TELEMETRIA::pisca(uint8_t porta) {
  digitalWrite(porta, HIGH);
  delay(1000);
  digitalWrite(porta, LOW);
  delay(1000);
}
