#include <iostream>
#include <fstream>
#include <stdexcept>

int main() {
    std::string imgNumberStr;
    std::ifstream myFile;
    myFile.open("./data/imgNumber.csv");

    if (myFile.is_open()) {
        getline(myFile, imgNumberStr);
        myFile.close();

        if (!imgNumberStr.empty() && std::isdigit(imgNumberStr[0])) {
            int imgNumber = std::stoi(imgNumberStr);
            std::cout << "The number of images is " << imgNumber << ".\n";
        } else {
            std::cerr << "Invalid or empty content in the file.\n";
        }
    } else {
        std::cerr << "Unable to open the file.\n";
    }

    return 0;
}
