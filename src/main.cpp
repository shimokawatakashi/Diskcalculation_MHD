#include "dust_disk.hpp"
#include "utils.hpp"

#include <iostream>
// #include <omp.h>
#include <fstream>
#include <chrono>


int main(int argc, char* argv[])
{
    bool restart = false;
    String input_file;
    String restart_gasfile;
    String restart_dustfile;
    String restart_Bfile;
    double restart_year;

    // if (argc != 2) {
    //     std::cerr << "# ERROR : main, no input file" << std::endl;
    //     std::exit(1);
    // }
    // input_file = argv[1];
    // std::cout << "input file : " << input_file << std::endl;
    
    if(std :: string(argv[1]) == "-r"){
        restart = true;
        input_file = argv[2];
        restart_dustfile = argv[3];
        restart_gasfile = argv[4];
        restart_Bfile = argv[5];
        restart_year = std :: stod(argv[6]);
    }
    else{
        input_file = argv[1];
    }
    
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    // ######################################################
    Utils::InputConfigure input;
    input.ReadFile(input_file);

    DustDisk dustdisk;
    dustdisk.InputParameter(input);
    
    if(!restart){
        dustdisk.Initialization();
    }
    else{
        dustdisk.Initialization_Restart(restart_dustfile,restart_gasfile,restart_Bfile,restart_year);
    }
    dustdisk.CalculationRun();

    // ######################################################

    end = std::chrono::system_clock::now();
    const double time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "time = " << (time*1.0e-3) << "sec" << std::endl;
    
    return 0;
}