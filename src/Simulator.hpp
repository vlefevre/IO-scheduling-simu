#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <functional>
#include <vector>
#include <chrono>
#include "Pattern.hpp"

#define NB_LP 20

class Machine
{
public:
    Machine(double B, double b) : B(B), b(b) {}
    double getBW() const { return B; }
    double getProcBW() const { return b; }
private:
    double B; //bqndwidth between I/O node and file system
    double b; //bandwidth per processor on I/O network
};

//try to insert a new instance of application k in the pattern pat, IO operations can be split
std::vector<Interval> splitSchedulable(Pattern& pat,const Machine& m, int k,std::chrono::microseconds& totalDuration);
std::vector<Interval> splitSchedulableV2(Pattern& pat,const Machine& m, int k,std::chrono::microseconds& totalDuration);
//one heuristic
void insertInSchedule(Pattern& pat, const Machine& m, std::function<std::vector<Interval>(Pattern&, const Machine&,int, std::chrono::microseconds&)>, double);

class Simulator
{
public:
    //construct a simulator object given a pattern (mostly empty), a machine and given heuristics
    Simulator(const Pattern &patt, const Machine &m, std::string heur = "insert", std::string mod="split") : pattern(patt), machine(m) {
        if (heur == "insert")
            heuristic = insertInSchedule;

        if (mod == "split")
            model = splitSchedulable;
        else if (mod == "splitv2")
            model = splitSchedulableV2;
    }
    void setHeuristic(std::string heur);
    void setModel(std::string mod);
    void run(double Tmin, double Tmax, double epsilon, std::string filename); //executes a heuristic on several periods and computes SysEfficiency and Diation
    void print_pattern() { pattern.print(); }
    void exportSchedule(std::string out) { pattern.exportSchedule(out); }
    double getSlowdown(); //Dilation
    double getPerf(); //SysEfficiency
    std::vector<double> getResult() { return {result_performance, result_dilation, result_period}; }
private:
    std::function<std::vector<Interval>(Pattern&, const Machine&,int,std::chrono::microseconds&)> model;
    std::function<void(Pattern&, const Machine&, std::function<std::vector<Interval>(Pattern&, const Machine&,int,std::chrono::microseconds&)>, double)> heuristic;
    Pattern pattern;
    Machine machine;

    double result_dilation;
    double result_performance;
    double result_period;
};

#endif
