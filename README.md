# Kontributor
Program Star Sensor FTMD. Dibuat oleh mahasiswa FTMD yang dibimbing oleh dosen sebagai tugas akhir.

## Penyusun program:
1. 13614063 - Nugraha Setya Ardi
2. 13616019 - Alif Ramadhana
3. 13616040 - I Wayan Fajar Surya Negara
4. 13617115 - Brian Mohammed Catraguna
5. 13618012 - David Carlos Xaverius Pardede
6. 13619021 - Ahmad Izza Falahuddin
7. 13620033 - Hilmi Nuruzzaman

## Pembimbing:
1. Dr. Eng. Ridanto Eko Poetro, S.T, M.Sc.
2. M. Arif Saifudin, S.T, M.Kom.
3. Ir. Toto Indriyanto, M.Sc., Ph.D.
4. Luqman Fathurrohim, S.T, M.T.

# Cara Menjalankan Program
## Cara set up program di WSL
1. Install Windows Subsystem for Linux.
2. Buat user WSL kemudian install anaconda.
3. Clone repository menuju folder di WSL.
4. Install library yang digunakan dalam program.
    - OpenCV
    - Python untuk C++
5. Pastikan setiap path yang ada dalam program sudah benar, terutama di bagian **main** dan **.vscode**.

## Cara menjalankan program star sensor di WSL
1. Perhatikan setiap path yang dipanggil di setiap baris kode pada folder **main** dan **.vscode**. Ganti menjadi path yang digunakan atau ganti menggunakan *relative path*.
    > Contoh kode yang perlu diperhatikan, ganti path menjadi path yang benar.
    ```
    string imgNumberStr;
    ifstream myFile;
    myFile.open("/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/main/imgNumber.csv");
    ```
2. Ambil satu gambar bintang yang sudah di-generate oleh program **star_simulator-main** kemudian pindahkan ke folder **main** dan rename file gambar menjadi *stars1.jpg*
3. Compile program **main/starcent.cpp** dengan menggunakan perintah `"build centroid"` di program **.vscode/tasks.json**.
    - Bisa juga melakukan compile dengan run command `make clear` kemudian run `make` di terminal.
4. Panggil program yang sudah selesai di build di terminal dengan menggunakan perintah `./main/centroid`.
5. Program akan mengeluarkan data ID, RA, DE, Roll, Error, dan total waktu di terminal.

## Cara membuat simulasi gambar bintang dari attitude yang diketahui
1. Buka program **star_simulator-main/starsim.py**
2. Edit attitude simulasi bintang menjadi attitude yang ingin di-generate pada baris kode berikut
    ```
    #Right ascension, declination and roll input prompt from user
    ra0 = 240
    de0 = 60
    roll0 = 30
    ```
3. Pada baris paling bawah, perhatikan penamaan dan output gambar dari program. Ubah variabel *file_name* menjadi nama file yang diinginkan. Atau uncomment bagian `#cv2.imwrite("main/stars1.jpg",background)` untuk langsung generate gambar di folder **main**.
```
file_name = f"ra{ra0}_de{de0}_roll{roll0}.jpg"
cv2.imwrite(file_name,background)
#cv2.imwrite("main/stars1.jpg",background)
```

# Cara menjalankan Program di Raspberry Pi
## Setup Raspberry Pi
1. Install OS Raspberry Pi pada SD card yang digunakan sebagai storage utama.
2. Clone repository ini ke Raspberry Pi.
3. Install semua library yang dibutuhkan, terutama OpenCV.
4. Pastikan setiap path yang dipanggil di setiap baris kode benar.

## Cara install OpenCV di Raspi
Install OpenCV di Raspi kadang membuat Raspinya error atau freeze, salah satu solusinya adalah dengan menambah swap space di Raspi.
1. Naikkan swap space di Raspi dengan menjalankan kode
    ```
    sudo nano /etc/dphys-swapfile
    ```
    Kemudian naikkan swap size menjadi 2048 MB pada bagian `CONF_SWAPSIZE=2048`. Lalu restart swap service di Raspi
    ```
    sudo /etc/init.d/dphys-swapfile stop
    sudo /etc/init.d/dphys-swapfile start
    ```
2. Install dependencies pada Raspi.
    ```
    sudo apt update
    sudo apt install -y build-essential cmake git pkg-config libjpeg-dev libtiff-dev libpng-dev \
        libavcodec-dev libavformat-dev libswscale-dev libv4l-dev libxvidcore-dev libx264-dev \
        libgtk2.0-dev libgtk-3-dev libatlas-base-dev gfortran python3-dev

    ```
3. Clone Repository OpenCV ke Raspi.
    ```
    cd ~
    git clone https://github.com/opencv/opencv.git
    git clone https://github.com/opencv/opencv_contrib.git
    ```
4. Buat direktori untuk melakukan build OpenCV.
    ```
    cd ~/opencv
    mkdir build
    cd build
    ```
5. Configure instalasi menggunakan CMake.
    ```
    cmake -D CMAKE_BUILD_TYPE=RELEASE \
        -D CMAKE_INSTALL_PREFIX=/usr/local \
        -D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib/modules \
        -D ENABLE_NEON=ON \
        -D ENABLE_VFPV3=ON \
        -D WITH_TBB=ON \
        -D BUILD_TBB=ON \
        -D BUILD_TESTS=OFF \
        -D WITH_EIGEN=OFF \
        -D WITH_GSTREAMER=OFF \
        -D WITH_V4L=ON \
        -D WITH_LIBV4L=ON \
        -D WITH_OPENGL=ON ..
    ```
6. Build OpenCV menggunakan make
    ```
    make -j2
    ```
7. Install OpenCV setelah melakukan build.
    ```
    sudo make install
    sudo ldconfig
    ```

## Cara menjalankan program di Raspi
1. Perhatikan setiap path yang dipanggil di setiap baris kode pada folder **main** dan **.vscode**. Ganti menjadi path yang digunakan atau ganti menggunakan *relative path*.
    > Contoh kode yang perlu diperhatikan, ganti path menjadi path yang benar.
    ```
    string imgNumberStr;
    ifstream myFile;
    myFile.open("/home/hilmi/star-sensor-ftmd/StarSensorFTMD_hilmi/StarSensorFTMD/main/imgNumber.csv");
    ```
2. Ambil satu gambar bintang yang sudah di-generate oleh program **star_simulator-main** kemudian pindahkan ke folder **main** dan rename file gambar menjadi *stars1.jpg*
3. Compile program **main/starcent.cpp** dengan run command `make clear` kemudian run `make` di terminal.
4. Panggil program yang sudah selesai di build di terminal dengan menggunakan perintah `./main/centroid`.
5. Program akan mengeluarkan data ID, RA, DE, Roll, Error, dan total waktu di terminal.

## Cara menjalankan program di Raspi dengan melakukan pengambilan gambar
1. Pastikan Raspi sudah mendukung program pengambilan gambar.
2. Buka program **main/starcent.cpp** kemudian uncomment bagian fungsi pengambilan gambar di raspi.
    ```
    void takePicture(const std::string& fileName) {
        std::string command = "libcamera-still -o " + fileName;
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Error: Failed to take picture" << std::endl;
        }
    }
    ```
3. Uncomment juga bagian `takePicture(main/stars1.jpg);` di program **main/starcent.cpp**.
4. Compile ulang program menggunakan `make clear` dan `make` ddi terminal.
5. Panggil program yang sudah di build dan program akan mengambil gambar di Raspi.