#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include "Simulator.hpp"
#include "Pattern.hpp"

#define SCALE 64
#define BWRATIO 100

void createOutput(std::string filename,std::string outname);

int main(int argc, char** argv)
{
    Application turbulence1(70*SCALE,128.2*BWRATIO,32768/SCALE); //in GB
    Application astrophysics(240*SCALE,423.4*BWRATIO,8192/SCALE);
    Application turbulence2(1.2*SCALE,235.8*BWRATIO,4096/SCALE);
    Application plasmaphysics(7554*SCALE,34304*BWRATIO,32768/SCALE);

    Machine machine(3*BWRATIO,1); //300/640 cores sending to have congestion

    std::vector<Application> real_apps = {plasmaphysics,turbulence2,turbulence2};

    std::ofstream out_data("results.txt",std::ios::out);

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
            std::cout << w << " " << io << " " << procs << "\n";
            apps.push_back(Application(w*SCALE,io*BWRATIO,procs/SCALE));
            Tmin = std::max(Tmin,w*SCALE+(io*BWRATIO)/(procs/SCALE));
        }
        app_file.close();

        //double Tmin = 494902;
        double Tmax = 10*Tmin, epsilon = 0.01;

        Simulator simu(Pattern(0,apps),machine,"insert","splitv2",0);

        std::string outname = "temp.txt";
        if (argc > 1)
        {
            outname = std::string(argv[1]);
        }
        if (argc > 2)
            if (std::string(argv[2]) == "opt")
                simu.setOpt(true);
        simu.run(Tmin,Tmax,epsilon,outname);
        auto res = simu.getResult();
        out_data << num_set-10 << " ";
        for (auto d:res)
            out_data << d << " ";
        out_data << "\n";
    }
    out_data.close();
    //simu.print_pattern();
    //simu.exportSchedule("expe_input_5_v3.txt");
}

std::vector<std::pair<int,int> > getZeroInt(std::string logname)
{
    std::ifstream in(logname,std::ios::in);
    std::string buff;
    std::vector<std::pair<int,int> > switches;
    while (buff != "Value")
        in >> buff;
    int currentzero = 1, currentnonzero = 0;
    in >> buff;
    while (buff != "All")
    {
        std::istringstream ss(buff.substr(2));
        ss >> currentnonzero;
        if (currentnonzero > currentzero)
        {
            switches.push_back(std::make_pair(currentzero,currentnonzero-currentzero));
        }
        currentzero = currentnonzero+1;
        in >> buff;
        in >> buff;
    }
    in >> buff; //other
    in >> buff; //variables
    in >> buff; //in
    in >> buff; //the
    in >> buff; //range
    in >> buff; //1-x
    std::istringstream ss(buff.substr(2));
    int nbVars; ss >> nbVars;
    if (currentzero != nbVars+1) //meaning last vars are 0
    {
        if (switches.size() > 0)
        {
            if (switches[0].first == 1)
                switches[0] = std::make_pair(currentzero,switches[0].second+nbVars+1-currentzero);
            else
                switches.push_back(std::make_pair(currentzero,nbVars+1-currentzero));
        }
        else
            switches.push_back(std::make_pair(currentzero,nbVars+1-currentzero));
    }
    return switches;
}

double readPeriod(std::string filename)
{
    std::string buff;
    std::ifstream in(filename,std::ios::in);
    while (buff != "Optimal:")
        in >> buff;
    double period;
    in >> buff >> buff >> period;
    return period;
}

void createOutput(std::string filename,std::string outname)
{
    std::ifstream lp(filename+".lp",std::ios::in);
    if (!lp)
    {
        std::cout << "File " << filename << ".lp does not exist.\n";
        return;
    }
    system(std::string("cplex -c \"r "+filename+".lp\" \"opt\" \"display solution variables -\" > cplex_out").c_str());
    std::ifstream sol("cplex_out",std::ios::in);

    std::vector<double> times;
    std::string buff;
    while (buff != "Value")
        sol >> buff;
    int index=0;
    double value=0,totalValue=0;
    int currentIndex = 0;
    while (buff != "All")
    {
        sol >> buff;
        std::istringstream convert(buff.substr(2));
        convert >> index;
        index--;
        for (int i=currentIndex; i<index; i++)
            times.push_back(totalValue);
        sol >> value;
        totalValue += value;
        times.push_back(totalValue);
        currentIndex = index+1;
    }
    //We complete with the end
    sol >> buff; sol >> buff; sol >> buff; sol >> buff; sol >> buff; sol >> buff;
    std::istringstream convert(buff.substr(2));
    convert >> index;
    for (int i=currentIndex; i<index; i++)
        times.push_back(totalValue);

    /* PROCESSING THE LP FILE TO OUTPUT INTERVALS */
    std::ofstream out(outname,std::ios::out);

    while (buff != "TO")
        lp >> buff;
    lp.ignore();
    lp.ignore();
    getline(lp,buff);
    while (buff != "BOUNDS")
    {
        out << "[ ";
        while (buff != "")
        {
            getline(lp,buff); //get the IO line
            std::stringstream ss(buff);
            double bw;
            std::string temp;
            int varIndex;
            out << "[ ";
            bool first = true;
            while (1)
            {
                ss >> bw >> temp; //Get bw and t_i
                std::stringstream convert(temp.substr(2));
                std::cout << bw << " " << temp << "\n";
                convert >> varIndex;
                varIndex--;
                double startT,endT;
                if (varIndex == 0)
                    startT = 0;
                else
                    startT = times[varIndex-1];
                endT = times[varIndex];
                if (endT > startT && bw > 0)
                {
                    if (!first)
                        out << ", ";
                    out << "( " << startT << " , ";
                    out << endT << " , ";
                    out << bw << " )";
                    first = false;
                }
                ss >> temp; //Get + or >=
                if (temp != "+")
                    break;
            }
            out << " ] ";
            getline(lp,buff); //get next Work or blank line for end of app
        }
        getline(lp,buff); //next app, first work
        out << " ]\n";
    }
    sol.close();
    lp.close();
    out.close();
}
