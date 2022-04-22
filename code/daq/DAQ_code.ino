#include <OneWire.h>
#include <DallasTemperature.h>
#include <SPI.h>
#include <SD.h>
File myFile;
DateTime now;
int newHour = 0;
int oldHour = 0;
void save_temper

void save_temperature() { 
  myFile = SD.open("temp.txt", FILE_WRITE);
  now = rtc.now();
  myFile.print(now.hour());
  myFile.print(":");
  myFile.print(now.minute());
  myFile.print(",");
  
  for (int i = 2; i <5; i + 2) {
    // Data wire is plugged into pin 2 on the Arduino
    #define ONE_WIRE_BUS i
 
    // Setup a oneWire instance to communicate with any OneWire devices 
    // (not just Maxim/Dallas temperature ICs)
    OneWire oneWire(ONE_WIRE_BUS);
 
    // Pass our oneWire reference to Dallas Temperature.
    DallasTemperature sensors(&oneWire);
 
    void setup(void) {
      // start serial port
      Serial.begin(9600);
      Serial.println("Dallas Temperature IC Control Library Demo");

      // Start up the library
      sensors.begin();
    }

    sensors.requestTemperatures(); // Send the command to get temperatures
    myFile.print("Temperature is: ");
    myFile.print(sensors.getTempCByIndex(0)); 
    }

  myFile.close();
  }

void setup ()
{
Wire.begin();
rtc.begin();
Serial.begin(9600);
Serial.print("Initializing SD card...");
if (!SD.begin(10)) {
Serial.println("initialization failed!");
while (1);
}
Serial.println("initialization done.");
now = rtc.now();
oldHour = now.hour();
}

void loop () {
  now = rtc.now();
  newHour = now.hour();
  if (oldHour != newHour) {
    save_temperature();
    oldHour = newHour;
  }
}
