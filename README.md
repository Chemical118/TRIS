TRIS : Telomeric Repeat motif Identification tool with Short-read sequencing
================
TRIS is tool that can identify Telomeric repeat motif (TRM) with short-read sequencing data. This tool looks for repeated sequences in a single read to find candidates for TRMs, iterating through them to finally find a TRM.

### Install
You can download a binary from release or build from source.

[Windows (x86_64)](https://github.com/Chemical118/TRIS/releases/latest/download/tris-windows-x86_64.tar.gz)  
[Linux (x86_64)](https://github.com/Chemical118/TRIS/releases/latest/download/tris-linux-x86_64.tar.gz)  
[MacOS (x86_64, arm64)](https://github.com/Chemical118/TRIS/releases/latest/download/tris-macos-universal.tar.gz)  

### Install from source
Windows (Visual Studio)
```sh
git clone https://github.com/Chemical118/TRIS.git
cd TRIS

git clone https://github.com/Microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg\scripts\buildsystems\vcpkg.cmake
cmake --build build --config Release

.\build\Release\tris_test.exe -- test\test.fastq.gz test\test.fastq
.\build\Release\tris.exe -h
```

Unix
```sh
git clone https://github.com/Chemical118/TRIS.git
cd TRIS

git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build

./build/tris_test.exe -- test/test.fastq.gz test/test.fastq
./build/tris -h
```

### Quick Start
```sh
tris MIN_MER MAX_MER <short_read_data1.fastq.gz> <short_read_data2.fastq>... -t <number of threads> 
```
MIN_MER : minimum length of sequence to find telomere [MIN_MER >= 3]  
MAX_MER : maximum length of sequence to find telomere [MAX_MER <= 64]  

The following is a recommended command line to run TRIS.
```sh
tris 5 32 <short_read_data1.fastq.gz> <short_read_data2.fastq>... -t <number of threads> 
```

### Output
```
><short_read_data1.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>
...

><short_read_data2.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>
...
```

### Author
Hyunwoo Ryu <wowo0118@korea.ac.kr>

*Special thanks to*  
Jiho Choi <sdatoli@korea.ac.kr>  
Kyungmo Ku <kyungmoku7141@gmail.com>