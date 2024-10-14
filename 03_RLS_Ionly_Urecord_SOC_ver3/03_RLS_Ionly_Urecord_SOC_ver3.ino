#include <Arduino.h>
#include <ADC.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm> 
#include <SD.h>

ADC *adc = new ADC(); // ADC object

// Set the pin numbers
const int readPin0 = A2; // Current sensing pin
const int readPin1 = A3; // Voltage sensing pin
const int TSpin0 =  27; // State signal pin to check the performance in oscilloscope. Voltage signal.
const int TSpin1 =  39; // Negative voltage control => IGBT control -
const int TSpin2 =  40;  // Positive voltage control => IGBT control +
const int redPin = A0; // to know if teensy works well or not in visual
const int greenPin = A1;
const int bluePin = A4;

// Basic settings
const int highDuration = 5000;  // Duration of high state in microseconds => IGBT turned on (current flows) for this time. Generally 2~3ms. need to discuss with Saitama uni.
const uint32_t freq = 50000; // ADC working frequency in Hz. If it's too big, there are over-calculation, which makes delay. 

// Variables for ADC values and status
volatile uint16_t adc_val0 = 0;
volatile uint16_t adc_val1 = 0;
volatile bool adc0_ready = false;
volatile bool adc1_ready = false;

const float sigma = 1.5f; // change based on the matlab (Online_exp.mlx file)
int sz = round(3 * sigma); // change based on the matlab (Online_exp.mlx file)
std::vector<float> gaussKernel; // Gaussian kernel for smoothing
std::vector<float> window(sz, 0.0); // Filtering window
//std::vector<float> voltWindow(sz, 0.0);

// State machine variables
int state = 0;

// Sliding window parameters
const int window_size_1 = 50; // Window length (ex. if freq 1e5 and window_size_1=50, it means one data set in the window_size_1 is treated in 0.5ms)
const int window_size_2 = 1000; // adjustable value
const int window_size = 20; // adjustable value
const double threshold1 = 50; // fault current threshold. adjustable value
double delta_T = 0.00002; // Time interval in seconds // adjustable, 1/freq
double ms_T = 0.02; // Time interval in milliseconds // to record in SD card => delta_T*1000 = 1/freq*1000, adjustable

// sliding window 찾아보기. 
std::vector<float> data_window(window_size_1, 0.0); // data processing function
std::vector<float> voltage_window(window_size_1, 0.0); // same as data_window but for voltage data
std::vector<float> peak_window(window_size, 0.0); // to detect peak current
std::vector<float> post_surge_data; // Stores current data after surge
std::vector<float> voltage_data; // stores voltage after surge
std::vector<float> time_since_surge; // Time tracking since surge
std::vector<float> A_estimates_rls; // RLS filtered data

// RLS parameters
double A_estimate_rls = 0;
double prev_A = 0.0f; // Previous estimate of parameter A
double P_rls = 1e6; // Initial error covariance
double lambda = 0.97; // Forgetting factor(0.95~1, if it is going to be closer 1, more accurate but delayed) => only adjustable value in RLS parameters
double z; // Observation
double K_rls; // Gain vector

// Other global variables
double H; // Current time
double T; // Millisecond time
double I; // Current
double U; // Voltage
double It0; // Peak Current
double current; // Raw Data from the sensor
double voltage; // Raw Data from the sensor
int k = 0;
double start_time; // Start time of the event
double peak_time; // Time at peak
double commutation; // Time of commutation event
volatile bool readyToWrite = false;
volatile bool dataWritten = false;
int mark = 1;

// SD card writing parameters
const size_t WRITE_BATCH_SIZE = 1000; // Number of data points to write at once
size_t dataIndex = 0; // Index for data writing

// Function declarations (implementation not shown in this snippet)
void updateDataWindow(float newData, int k);
void updateVoltageWindow(float newData, int k);
bool allGreaterOrEqual(const std::vector<float>& window, float threshold);
void writeDataToSDCard2();
void processData();

// Function to create Gaussian kernel
std::vector<float> createGaussKernel(float sigma, int sz) {
    std::vector<float> gaussKernel(sz);
    float sum = 0.0;
    for(int i = 0; i < sz; ++i) {
        float x = i - sz / 2;
        gaussKernel[i] = exp(-x * x / (2 * sigma * sigma));
        sum += gaussKernel[i];
    }
    for(auto& value : gaussKernel) {
        value /= sum;
    }
    return gaussKernel;
}

// Function to update the data window with a new data point
void updateDataWindow(float newData, int k) {
    // Shift the current data window to make space for the new data
    if (k <= window_size_1) {
        data_window[k - 1] = newData;
    } else {
        std::rotate(data_window.begin(), data_window.begin() + 1, data_window.end());
        data_window.back() = newData;
    }
}

void updateVoltageWindow(float newData, int k) {
    // Shift the voltage data window to make space for the new data
    if (k <= window_size_1) {
        voltage_window[k - 1] = newData;
    } else {
        std::rotate(voltage_window.begin(), voltage_window.begin() + 1, voltage_window.end());
        voltage_window.back() = newData;
    }
}

void updatePeakWindow(float newData, int k) {
    // Similar to updateDataWindow, but for the peak detection window
    if (k <= window_size) {
        peak_window[k - 1] = newData;
    } else {
        std::rotate(peak_window.begin(), peak_window.begin() + 1, peak_window.end());
        peak_window.back() = newData;
    }
}

// Function to check if all elements in a window are greater than or equal to a threshold
bool allGreaterOrEqual(const std::vector<float>& window, float threshold) {
    // Iterate through the window and check each value against the threshold
    for (float value : window) {
        if (value < threshold) {
            return false; // If any value is below the threshold, return false
        }
    }
    return true; // All values are above the threshold
}

void writeDataToSDCard2() {
  // Creates a .csv file and records data
  String fileName = "data_" + String(millis()) + ".csv";
  File dataFile = SD.open(fileName.c_str(), FILE_WRITE);

  if (dataFile) {
    dataFile.print("Start Time:,");
    dataFile.println(start_time);
    dataFile.print("Peak Time:,");
    dataFile.println(peak_time);
    dataFile.print("Commutation Time:,");
    dataFile.println(commutation);

    dataFile.println("Post Surge Data,Voltage Data,Time Since Surge,A Estimates RLS");

    while (dataIndex < post_surge_data.size()) {
        for (size_t i = dataIndex; i < std::min(dataIndex + WRITE_BATCH_SIZE, post_surge_data.size()); ++i) {
            dataFile.print(post_surge_data[i]); dataFile.print(",");
            dataFile.print(voltage_data[i]); dataFile.print(",");
            dataFile.print(time_since_surge[i]); dataFile.print(",");
            if (A_estimates_rls.size() > i) { 
                dataFile.print(A_estimates_rls[i]);
            }
            dataFile.println();
        }
        dataIndex += WRITE_BATCH_SIZE;
    }
    dataFile.close();
    Serial.println("CSV data is written on SD");
  } else {
    Serial.println("Failed to open file, retrying...");
    delay(100); // Check again after delay
  }
}

void setup() {
  // Pin setup and initial states
  pinMode(TSpin0, OUTPUT); // oscilloscope
  pinMode(TSpin1, OUTPUT); // IGBT -
  pinMode(TSpin2, OUTPUT); // IGBT +

  pinMode(redPin, OUTPUT); // LED red
  pinMode(greenPin, OUTPUT); // LED green
  pinMode(bluePin, OUTPUT); // LED blue (not used)

  digitalWrite(TSpin0, LOW);
  digitalWrite(TSpin1, HIGH);
  digitalWrite(TSpin2, LOW);

  digitalWrite(redPin, HIGH);
  digitalWrite(greenPin, HIGH);
  digitalWrite(bluePin, HIGH);

  // Serial communication setup
  Serial.begin(38400);
  while (!Serial && millis() < 5000); // Wait for serial monitor

  // Create the Gaussian kernel for filtering
  gaussKernel = createGaussKernel(sigma, sz);

  // ADC setup for both channels with high speed and resolution settings
  adc->adc0->setAveraging(8); // the number of data to find the average (do not change the value) 
  adc->adc0->setResolution(10); // 
  adc->adc0->setConversionSpeed(ADC_CONVERSION_SPEED::VERY_HIGH_SPEED);
  adc->adc0->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_HIGH_SPEED);

  adc->adc1->setAveraging(8); 
  adc->adc1->setResolution(10); 
  adc->adc1->setConversionSpeed(ADC_CONVERSION_SPEED::VERY_HIGH_SPEED); 
  adc->adc1->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_HIGH_SPEED);

// Start single read and enable interrupts for ADC channels
  adc->adc0->startSingleRead(readPin0);
  adc->adc0->enableInterrupts(adc0_isr);
  adc->adc0->startTimer(freq); 

  adc->adc1->startSingleRead(readPin1);
  adc->adc1->enableInterrupts(adc1_isr);
  adc->adc1->startTimer(freq); 

  adc->startSynchronizedSingleRead(readPin0, readPin1);

  // SD card
  if (!SD.begin(BUILTIN_SDCARD)) {
        Serial.println("SD Error");
        return;
    }
    Serial.println("SD Ready");
    Serial.println("Start");
}

void loop() {
  // Continuously checks for ADC data and processes it when ready
  // Manages state transitions and executes RLS algorithm
  // Handles data writing to SD card when required
  if (adc0_ready && adc1_ready) {

    std::rotate(window.begin(), window.begin() + 1, window.end());
    window.back() = current;

    //std::rotate(voltWindow.begin(),voltWindow.begin() + 1, voltWindow.end());
    //voltWindow.back() = voltage;


    float filteredValue = std::inner_product(window.begin(), window.end(), gaussKernel.begin(), 0.0f);
    I = filteredValue;

    //float filteredVoltValue = std::inner_product(voltWindow.begin(), voltWindow.end(), gaussKernel.begin(), 0.0f);
    //U = filteredVoltValue;

    processData();

    adc0_ready = false;
    adc1_ready = false;

    adc->startSynchronizedSingleRead(readPin0, readPin1);
  }

  if(state==1 && mark ==1 ){
    Serial.print("Fault come and time is");
    Serial.println(start_time);
    digitalWrite(TSpin0, HIGH);
    mark = mark+1;
  }
  if(state==2 && mark ==2 ){
    Serial.print("Peak come and time is");
    Serial.println(peak_time);
    digitalWrite(TSpin0, LOW);
    mark = mark+1;
  }
  if (readyToWrite) {
    Serial.print("Commutation time is");
    Serial.println(commutation);
        writeDataToSDCard2(); 
        readyToWrite = false; // remark
    }

}

void adc0_isr() {
  //digitalWrite(TSpin0, HIGH);
  adc_val0 = adc->adc0->readSingle();
  current = (adc_val0 * (3.3 / adc->adc0->getMaxValue()) - 0.001)* 2000 ; // 2000 is the value for Rogowski coil in current sensor
  
  if(current < 0){
    current = 0;
  }

  adc0_ready = true;
  
  #if defined(__IMXRT1062__)  // Teensy 4.0
    asm("DSB");
  #endif
}

void adc1_isr() {
    adc_val1 = adc->adc1->readSingle();
    voltage = (adc_val1 * (3.3 / adc->adc1->getMaxValue()) - 0.001)* 100 ;

    if(voltage < 0){
    voltage = 0;
  }

    adc1_ready = true;
}

// The processData function contains the core logic for analyzing the incoming data.
void processData() {
    // Increment the data index (k) for each call
    k = k + 1;

    // Update windows with the latest current and voltage value (I and U)
    updateDataWindow(I, k);
    updatePeakWindow(I, k);
    updateVoltageWindow(voltage, k);
    
    // State machine for processing data
    // finding fault current
    if (state == 0) {
        // In state 0, the system is waiting for a surge in current
        // Extract a subset of data_window for checking against threshold1
        std::vector<float> window_2(data_window.end() - window_size_2, data_window.end()); // make another window because data window is too long, so not accurate.

        // Check if the last few readings are above the threshold1
        if (allGreaterOrEqual(window_2, threshold1)) {
            // Loop through the data window to find the point of sudden increase
            digitalWrite(TSpin0, HIGH);
            post_surge_data = std::vector<float>(data_window.end() - window_size_2, data_window.end());
            voltage_data = std::vector<float>(voltage_window.end() - window_size_2, voltage_window.end());

            // Initialize the time_since_surge array with time data
            time_since_surge.resize(post_surge_data.size());
            for (size_t i = 0; i < post_surge_data.size(); ++i) {
            time_since_surge[i] = i * ms_T;
            }
            // Calculate the total time duration H and convert to milliseconds (T)
            H = (post_surge_data.size() - 1) * delta_T;
            T = H * 1e3;
            start_time = T;

            // Update state and control pins
            state = 1;
            Serial.print("Voltage1 is");
            Serial.println(voltage);
            Serial.print("Current1 is");
            Serial.println(current);
        }
    }

    // In state 1, the system has detected a surge and is now looking for a peak current
    if (state == 1) {
        // Check if the start of peak_window has a lower sum than the end, indicating a peak
        if (std::accumulate(peak_window.begin(), peak_window.begin() + 5, 0.0f) >
            std::accumulate(peak_window.end() - 5, peak_window.end(), 0.0f)) {
            state = 2; // Change to state 2 for post-peak processing
            Serial.print("Voltage2 is");
            Serial.println(voltage);
            Serial.print("Current2 is");
            Serial.println(current);
            digitalWrite(TSpin0, LOW);
            digitalWrite(redPin, LOW);
            // Determine peak current
            It0 = *std::max_element(peak_window.begin(), peak_window.end());
            peak_time = T; // Record the peak time
            // Resize the RLS estimates array
            A_estimates_rls.resize(post_surge_data.size(), 0.0f);
        } else {
            // Continue accumulating data post-surge
            post_surge_data.push_back(I);
            voltage_data.push_back(voltage);

            // Update time variables
            H += delta_T;
            T += ms_T;
            time_since_surge.push_back(T);
        }
    }

    // State for handling RLS algorithm and commutation detection
    // when the trigger signal occurs
    if (state == 2){
        // Append current reading to post-surge data and update time variables
        post_surge_data.push_back(I);
        voltage_data.push_back(voltage);
        H += delta_T;
        T += ms_T;
        time_since_surge.push_back(T);

        // RLS algorithm implementation
        z = std::log(It0) - std::log(I + 1e-6); // Current observation
        K_rls = P_rls * H / (lambda + H * P_rls * H); // Calculate gain vector
        A_estimate_rls += K_rls * (z - H * A_estimate_rls); // Update estimate of A
        P_rls = (1 - K_rls * H) * P_rls / lambda; // Update error covariance matrix

        // Check conditions for commutation and update state
        if (H > 1e-3 && A_estimate_rls > 50 && (A_estimate_rls - prev_A) < 0) {
            commutation = T; // Commutation time in milliseconds
            digitalWrite(TSpin0, HIGH);
            digitalWrite(TSpin1, LOW);
            digitalWrite(TSpin2, HIGH);
            digitalWrite(redPin, HIGH);
            Serial.print("IGBT ON by RLS");
            delayMicroseconds(highDuration);
            digitalWrite(TSpin0, LOW);
            digitalWrite(TSpin1, HIGH);
            digitalWrite(TSpin2, LOW);
            digitalWrite(greenPin, LOW);
            Serial.println("IGBT OFF");
            Serial.print("Voltage3 is");
            Serial.println(voltage);
            Serial.print("Current3 is");
            Serial.println(current);
            state = 3;
        } else {
            prev_A = A_estimate_rls;
        }

        // Store the current estimate of A
        A_estimates_rls.push_back(A_estimate_rls);
    }

    // State for finalizing data after commutation
    if (state == 3 && !dataWritten){
        // Check if the required data has been collected
        if (H >= 0.05) {
            state = 4;
        } else {
            // Continue collecting data and applying RLS algorithm
            post_surge_data.push_back(I);
            voltage_data.push_back(voltage);
            H += delta_T;
            T += ms_T;
            time_since_surge.push_back(T);

            // Continue RLS algorithm with the latest data
            z = std::log(It0) - std::log(I + 1e-6);
            K_rls = P_rls * H / (lambda + H * P_rls * H);
            A_estimate_rls += K_rls * (z - H * A_estimate_rls);
            P_rls = (1 - K_rls * H) * P_rls / lambda;

            // Store the current estimate of A
            A_estimates_rls.push_back(A_estimate_rls);
        }
    }

    // State to indicate readiness for data writing
    if (state == 4) {
        readyToWrite = true;
        dataWritten = true;
        state = 5;
    }
    
  if(state == 5){
  }

  }