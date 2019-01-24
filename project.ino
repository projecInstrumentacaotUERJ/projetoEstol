#include "HX711.h" 
#include "Servo.h"
 
// Criar um Objeto Servo
Servo servo1; 

#define DOUT 2                      
#define CLK  3                       
 
HX711 bau(DOUT, CLK);             // instancia Balança HX711
 
float calibration_factor = 34730;     // fator de calibração aferido na Calibração 
 
void setup()
{
  Serial.begin(9600); 
  bau.set_scale(calibration_factor);             // ajusta fator de calibração
  bau.tare();        
    servo1.attach(5); 
    servo1.write(0);
  // zera a Balança
}
 
void loop()
{
  // Lê o valor do Potenciometro
  int angle = analogRead(0); 
  // Mapeia o valor de 0 a 180 graus
  angle=map(angle, 0, 1023, 0, 180);
  // Repassa o angulo ao ServoWrite
  servo1.write(angle); 
   Serial.print(" - ");
  Serial.print(angle);



  
  Serial.print("Peso: ");                            
  Serial.print(bau.get_units(), 2);         // imprime peso na balança com duas casas decimais 
  Serial.println(" kg");                             // imprime no monitor serial 
  delay(500) ;                                       
  if (Serial.available())                            // se a serial estiver disponivel
  {
    char temp = Serial.read();                       // le carcter da serial 
    if (temp == 't' || temp == 'T')                  // se pressionar t ou T
    {
      bau.tare();                                // zera a balança
      Serial.println(" Balança zerada");             // imprime no monitor serial
    }
  }
}
