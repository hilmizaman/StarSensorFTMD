
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

int main() {
  
  auto start = high_resolution_clock::now();
  //cout << "Method 2: Weighted Center of Gravity\n";
  // READING IMAGE NUMBER ON PYTHON

    string imgNumberStr;
    ifstream myFile;
    myFile.open("./program_starcentroiding/imgNumber.csv");
    
    while(myFile.good()) {
      string line;
      getline(myFile, line, ',');
      imgNumberStr = line;
    }
    myFile.close();
    int imgNumber = stoi(imgNumberStr);
    //cout << "The number of the images is " << imgNumber << ".\n";

  // CREATE NECESSARY ARRAYS FOR ITERATIONS

    // ARRAYS OF IMAGE PATHS

      string imgFolderPath = "./result_starcentroiding/starImages2/";
      string imgFileNameTitle = "stars";
      string imgExtension = ".png";
      int imgIndexInt = 0;
      string imgIndexStr = "";
      string imgPath = "";
      string imgPathArray[imgNumber] = {};

      for (int i = 0; i < imgNumber; i++) {
        imgIndexInt = i + 1;
        imgIndexStr = to_string(imgIndexInt);
        imgPath = imgFolderPath + imgFileNameTitle + imgIndexStr + imgExtension;
        imgPathArray[i] = imgPath;
      }

    // ARRAYS OF RESULT PATHS
      string resultFolderPath = "./result_starcentroiding/starPositionCalculated/secondMethod/WCG";
      string resultExtension = ".csv";
      int resultIndexInt = 0;
      string resultIndexStr = "";
      string resultPath = "";
      string resultPathArray[imgNumber] = {};

      for (int i = 0; i < imgNumber; i++) {
        resultIndexInt = i + 1;
        resultIndexStr = to_string(resultIndexInt);
        resultPath = resultFolderPath + resultIndexStr + resultExtension;
        resultPathArray[i] = resultPath;
      }

  for (int bigIndex = 0; bigIndex < imgNumber; bigIndex++) {

    //cout << "\nSTARS IMAGE NUMBER " << bigIndex + 1 << "\n";

    // 1. BASIC DECLARATION AND DEFINITION
    
        // Defining the Image and Matrix to Store the Pixel Value
        cv::Mat img = cv::imread(imgPathArray[bigIndex]);      // Loading a Star Image
        cv::Mat stars;                                  // Declaring an Empty Matrix to Store Converted Image
        cv::cvtColor(img, stars, cv::COLOR_BGR2GRAY);   // Converting loaded image to grayscale image
        
        int rows = img.rows;
        int columns = img.cols;
        int toCheck[rows][columns] = {};
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                toCheck[i][j] = (int)stars.at<uchar>(i,j);
                //if (bigIndex == 5) {cout << toCheck[i][j] << " ";}
            }   //if (bigIndex == 5) {cout << endl;}
        } 

        //if (bigIndex == 7) {cout << toCheck[353][324] << endl;}     

    // 2. FINDING THE STAR CENTROID

        // Declaring Region of Interest
        int ROI_size_default = 8;                     // The size of ROI rows and columns
        int ROI_size = ROI_size_default;

        // Declaring Signal and Noise Threshold
        int signalThreshold = 50;  // Minimum brightness of a pixel to be identified
        int noiseThreshold = 15;    // Arbitrary number significantly lower than the signal threshold

        // Declaring the Required Values for Centroiding Calculation
        float DN = 0.0;        
        float xPix = 0.0;        
        float xCent = 0.0;        
        float yPix = 0.0;        
        float yCent = 0.0;
              
        vector<float> centroidsX;
        vector<float> centroidsY;
        int centroidsNum = 0;

        // Declaring Required Array and Vectors
        int isChecked[rows][columns];       // Array to store binary values to check whether a pixel has been checked or not
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                isChecked[i][j] = 0;
            }
        }

        // Iterating Through Each Pixels
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) { 
                
                /*if (bigIndex == 7 && i > 320 && j > 320 && i < 360 && j < 360) {
                            cout << i+1 << ", " << j+1 << ". " << toCheck[i][j] << endl;
                          }*/

                if (isChecked[i][j] == 0) {                 // Checking whether the pixel has been checked before
                    if (toCheck[i][j] >= signalThreshold) {  // Checking whether the pixel value is above the signal threshold

                          //if (bigIndex == 5) {cout << i+1 << ", " << j+1 << ". " << toCheck[i][j] << endl;}

                          if (i <= 6 || j <= 6 || i >= 497 || i >= 663) {
                            ROI_size = 6;
                            if (i <= 4 || j <= 4 || i >= 499 || i >= 665) {
                              ROI_size = 4;
                              if (i <= 2 || j <= 2 || i >= 501 || i >= 667) {
                                ROI_size = 2;
                              }
                            } 
                          }   else {ROI_size = ROI_size_default;}

                           for (int i2 = 0; i2 < (ROI_size_default); i2++) {
                              for (int j2 = 0; j2 < (ROI_size_default); j2++) {
                                isChecked[i - (ROI_size_default/2) + i2][j - (ROI_size_default/2) + j2] = 1;
                                }
                            }

                      int ROI[ROI_size][ROI_size] = {};     // Region of Interest
                      for (int i = 0; i < ROI_size; i++) {
                        for (int j = 0; j < ROI_size; j++) {
                          ROI[i][j] = 0;
                        }
                      }


                        // AT THIS POINT, THE ROI HAS BEEN CREATED, BUT IT IS STILL EMPTY


                        for (int i2 = 0; i2 < ROI_size; i2++) {      // Iterating through the premade ROI
                            for (int j2 = 0; j2 < ROI_size; j2++) {
                                ROI[i2][j2] = toCheck[i - ROI_size/2 + i2][j - ROI_size/2 + j2];    // Assigning the values of the pixels to the ROI
                                isChecked[i - ROI_size/2 + i2][j - ROI_size/2 + j2] = 1;
                            }
                        }  

                        // AT THIS POINT, THE ROI HAS BEEN FILLED WITH NUMBERS FROM THE stars        
                        // LET'S CHECK IT OUT FIRST!!     
                        // CHECKED! IT IS GOOD!

                        // Calculating the Sum of ROI Brightness
                        for (int i2 = 0; i2 < ROI_size; i2++) {
                            for (int j2 = 0; j2 < ROI_size; j2++) {
                                if (ROI[i2][j2] >= noiseThreshold) {
                                    DN = DN + ((ROI[i2][j2] - noiseThreshold)*(ROI[i2][j2] - noiseThreshold));  
                                }
                            }
                        }

                        // Calculating Total Horizontal (xPix) and Vertical (yPix) Pixels Value
                        for (int i2 = 0; i2 < ROI_size; i2++) {
                            for (int j2 = 0; j2 < ROI_size; j2++) {
                                if (ROI[i2][j2] >= noiseThreshold) {
                                    xPix = xPix + (j - ROI_size/2 + j2) * ((ROI[i2][j2] - noiseThreshold)*(ROI[i2][j2] - noiseThreshold));
                                    yPix = yPix + (i - ROI_size/2 + i2) * ((ROI[i2][j2] - noiseThreshold)*(ROI[i2][j2] - noiseThreshold));
                                }                             
                            }  
                        }

                        // Calculating Centroids
                        xCent = (xPix/DN) + 0.5;
                        yCent = (yPix/DN) + 0.5;
                        
                        DN = 0.0;
                        xPix = 0.0;
                        yPix = 0.0;

                        centroidsX.push_back(xCent);
                        centroidsY.push_back(yCent);
                        centroidsNum++;              
                    }
                    isChecked[i][j] = 1;
                }
            }
        }           

        // Rescaling the Centroids
        float xfactor = 1;
        float yfactor = 1;

        for (int i = 0; i < centroidsX.size(); i++) {
            centroidsX[i] = (centroidsX[i] * xfactor);
            centroidsX[i] = centroidsX[i] - (3280.0/2.0);
            centroidsX[i] = centroidsX[i] - 3.6;
            centroidsY[i] = (centroidsY[i] * yfactor);
            centroidsY[i] = (centroidsY[i] * -1.0) + (2464.0/2.0);
            centroidsY[i] = centroidsY[i] - 1.36;
        }

        //cout << centroidsNum << "\n";

        // Displaying the Centroids
        float centroids[centroidsNum][2] = {};
        for (int i = 0; i < centroidsX.size(); i++) {
          centroids[i][0] = centroidsX[i];
          centroids[i][1] = centroidsY[i];
          //if (bigIndex == 15) {cout << i+1 << ". (" << centroids[i][0] << ", " << centroids[i][1] << ")\n";}                                                                          
        }
        
        // Saving the Calculation Results into CSV File
        ofstream out(resultPathArray[bigIndex]);
        for (auto& row: centroids) {
          for (auto col : row)
            out << col << ',';
          out << '\n';
        }
  }

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  int timeConsumption = duration.count();

  cout << "2. Weighted Center of Gravity Time Consumption: " << timeConsumption/imgNumber << " microseconds\n";

    return 0;
}
