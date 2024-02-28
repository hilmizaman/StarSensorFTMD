#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::string imgNumberStr;
    std::ifstream myFile("./program_starcentroiding/imgNumber.csv");

    if (myFile.is_open()) {
        if (std::getline(myFile, imgNumberStr)) {
            try {
                int imgNumber = std::stoi(imgNumberStr);
                std::cout << "The number of images is: " << imgNumber << std::endl;
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid integer format in file" << std::endl;
            }
        } else {
            std::cerr << "Error: Failed to read from file" << std::endl;
        }
        myFile.close();
    } else {
        std::cerr << "Error: Failed to open file" << std::endl;
    }

    return 0;
}
