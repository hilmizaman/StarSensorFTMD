centroid: main/starcent.cpp
	g++ -fdiagnostics-color=always -g main/starcent.cpp -o main/centroid -I/usr/local/include/opencv4 -I/usr/include/python3.10 -L/usr/local/lib -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc -lpython3.10 
clear:
	rm -f main/centroid