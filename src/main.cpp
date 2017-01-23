#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include "Simulator.hpp"
#include "Pattern.hpp"

#define SCALE 64
#define BWRATIO 100


int main(int argc, char** argv)
{
    Application turbulence1(70*SCALE,128.2*BWRATIO,32768/SCALE); //in GB
    Application astrophysics(240*SCALE,423.4*BWRATIO,8192/SCALE);
    Application turbulence2(1.2*SCALE,235.8*BWRATIO,4096/SCALE);
    Application plasmaphysics(7554*SCALE,34304*BWRATIO,32768/SCALE);

    Machine machine(3*BWRATIO,1); //the total number of cores is not important, it is determined by the applications used

    std::string outname = "results.txt";
    if (argc > 1)
    {
        outname = std::string(argv[1]);
    }
    std::ofstream out_data(outname,std::ios::out);

    for (int num_set=11; num_set<=110; num_set++)
    {
        std::ostringstream stm;
        stm << num_set;
        std::ifstream app_file("sets/set"+stm.str(),std::ios::in);
        double w,io; int procs;
        double Tmin = 0;
        std::vector<Application> apps;
        while (app_file >> w >> io >> procs)
        {
            apps.push_back(Application(w*SCALE,io*BWRATIO,procs/SCALE));
            Tmin = std::max(Tmin,w*SCALE+(io*BWRATIO)/(procs/SCALE));
        }
        app_file.close();

        //Tmin = 495000;
        double Tmax = 10*Tmin, epsilon = 0.01;

        Simulator simu(Pattern(0,apps),machine,"insert","splitv2");

        simu.run(Tmin,Tmax,epsilon,"temp.txt");
        auto res = simu.getResult();
        out_data << num_set << " ";
        for (auto d:res)
            out_data << d << " ";
        out_data << "\n";
    }
    out_data.close();
}
