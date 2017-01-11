#include "Simulator.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <sstream>
#define DEBUG(X) std::cout << X << "\n"; std::cout.flush();

std::vector<std::pair<int,int> > getZeroInt(std::string logname);
double readPeriod(std::string filename);
double myFmod(double nb, double mod);

double Simulator::getPerf()
{
    std::vector<double> Perf(pattern.getNbApps());
    for (int k=0; k<Perf.size(); k++)
        Perf[k] = pattern.getProcs(k)*(pattern.getNbInstances(k)*pattern.getWork(k)/pattern.getPeriod())/640;
    return std::accumulate(Perf.begin(),Perf.end(),0.);
}

double Simulator::getSlowdown()
{
    std::vector<double> Sd(pattern.getNbApps());
    for (int k=0; k<Sd.size(); k++)
    {
        Sd[k] = pattern.getNbInstances(k)*(pattern.getWork(k)+pattern.getIOVol(k)/std::min(machine.getProcBW()*pattern.getProcs(k),machine.getBW()))/pattern.getPeriod();
    }
    return *std::min_element(Sd.begin(),Sd.end());
}

void Simulator::run(double Tmin, double Tmax, double epsilon, std::string filename)
{
    std::ofstream file(filename,std::ios::out);
	//Initialized variables to memorize best values
    double T = Tmin;
    double TmaxPerf,TminSd;
    double maxPerf = -1, minSd = -1;
    double SdmaxPerf = -1;
    std::chrono::microseconds total_duration(0);
    std::chrono::microseconds total_duration_IIS(0);
    while (T <= Tmax)
    {
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::high_resolution_clock::time_point t2;
        //DEBUG("T = " << std::setprecision(10) << T) //console display
        file << "T = " << T << "\n";
        pattern.clear(T); //initialize the pattern
        t1 = std::chrono::high_resolution_clock::now();
        heuristic(pattern,machine,model,epsilon); //run the insertion heuristic
        t2 = std::chrono::high_resolution_clock::now();
        total_duration_IIS += std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
        //print_pattern();
        file << "Performance : " << getPerf() << "\n";
        if (getPerf() > maxPerf + 0.000001)
        {
            maxPerf = getPerf();
            SdmaxPerf = getSlowdown();
            TmaxPerf = T;
        }
        file << "Slowdown : " << 1/getSlowdown() << "\n";
        if (getSlowdown() > minSd + 0.000001)
        {
            minSd = getSlowdown();
            TminSd = T;
        }
        for (int k=0; k<pattern.getNbApps(); k++)
        {
            double yield = pattern.getNbInstances(k)*pattern.getWork(k)/pattern.getPeriod();
           file << "\tApp " << k << " has " << pattern.getNbInstances(k) << " instances, for a yield of " << yield << ".\n";
        }
        T = round(T*(1+epsilon)*1000)/1000.;
    }
    T = TmaxPerf;
    double minReducedPeriod = T;
    while (true)
    {
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::high_resolution_clock::time_point t2;
        //DEBUG("T = " << std::setprecision(10) << T) //console display
        file << "T = " << T << "\n";
        pattern.clear(T); //initialize the pattern
        t1 = std::chrono::high_resolution_clock::now();
        heuristic(pattern,machine,model,epsilon); //run the insertion heuristic
        t2 = std::chrono::high_resolution_clock::now();
        total_duration_IIS += std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
        //print_pattern();
        file << "Performance : " << getPerf() << "\n";
        file << "Slowdown : " << 1/getSlowdown() << "\n";
        for (int k=0; k<pattern.getNbApps(); k++)
        {
            double yield = pattern.getNbInstances(k)*pattern.getWork(k)/pattern.getPeriod();
           file << "\tApp " << k << " has " << pattern.getNbInstances(k) << " instances, for a yield of " << yield << ".\n";
        }
        //This, in practice, checks if the number of instances changes in the pattern
        if (getPerf() > maxPerf*TmaxPerf/T - PRECISION && getPerf() < maxPerf*TmaxPerf/T + PRECISION)
        {
            minReducedPeriod = T;
            T = T - (TmaxPerf - TmaxPerf/(1+epsilon))/(1/epsilon);
        } else {
            break;
        }
    }
    file << "------------------------\n";
    file << "Max Performance of " << maxPerf << " for T = " << TmaxPerf << "\n";
    file << "Min Slowdown of " << 1/minSd << " for T = " << TminSd << "\n";
    file << "Max Reduced Performance of " << maxPerf*TmaxPerf/minReducedPeriod << " with dilation of " << 1/(SdmaxPerf*TmaxPerf/minReducedPeriod) << ", for T = " << minReducedPeriod << "\n";
    result_dilation = 1/(SdmaxPerf*TmaxPerf/minReducedPeriod);
    result_performance = maxPerf*TmaxPerf/minReducedPeriod;
    result_period = minReducedPeriod;
    file.close();
    std::cout << "Total time : " << total_duration_IIS.count() << " us.\n";
}

void Simulator::setHeuristic(std::string heur)
{
    if (heur == "insert")
        heuristic = insertInSchedule;
}

void Simulator::setModel(std::string mod)
{
	if (mod == "split")
        model = splitSchedulable;
}

std::vector<Interval> splitSchedulable(Pattern& pat, const Machine& m, int k, std::chrono::microseconds& totalDuration)
{
    double T = pat.getPeriod();
    double Tmin = T+1;
    std::vector<Interval> IOintList;
    std::vector<Interval> FinalIOintList;
    int i_ref = pat.getNbInstances(k)+1;
    auto itInstance = pat.getInstances(k).begin();
    for (int i=0; i<pat.getNbInstances(k) || (pat.getNbInstances(k)==0 && i==0); i++)
    {

        //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //We compute the first moment where the new instance can send IO
        double dateBeginIo = 0;
        //And the last moment = beginning of next instance
        double dateEndIo = T;

        if (pat.getNbInstances(k) > 0)
        {
            dateBeginIo = pat.getEndIOTransferDate(k,itInstance);
            itInstance++;
            if (itInstance == pat.getInstances(k).end())
                itInstance = pat.getInstances(k).begin();
            dateEndIo = *itInstance;
        }
        bool placeToWork = false;
        if (dateEndIo > dateBeginIo+PRECISION)
        {
            if (dateEndIo >= dateBeginIo + pat.getWork(k) - PRECISION)
                placeToWork = true;
        } else if (dateEndIo+PRECISION < dateBeginIo) {
            if (dateEndIo+T >= dateBeginIo + pat.getWork(k) - PRECISION)
                placeToWork = true;
        }
		//Special case when there is no instance
        if (pat.getNbInstances(k) > 0)
            dateBeginIo = myFmod(dateBeginIo+pat.getWork(k),T);

        //std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        //std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );
        //totalDuration += duration;
        if (placeToWork || (pat.getNbInstances(k) == 0 && T >= pat.getWork(k)) )
        {
            //DEBUG(dateBeginIo << " " << dateEndIo)
            pat.addDate(dateBeginIo);//we need to split the interval
            auto itBegin = pat.getDates().find(dateElement(dateBeginIo));
            auto itEnd = pat.getDates().find(dateElement(dateEndIo));
            ///Compute the bandwidths available
            int nbIntervals = (std::distance(itBegin,pat.getDates().end())+std::distance(pat.getDates().begin(),itEnd))%pat.getDates().size();
            IOintList = std::vector<Interval>(nbIntervals);
            double dataLeft = pat.getIOVol(k);
            double Time = 0;
            int cpt = 0; //to stop when used all possible intervals
            double prevMax = m.getBW()+1000;
            int prevI = -1;
            while (dataLeft > 0+PRECISION && cpt<nbIntervals)
            {
                ///We select the interval with max bandwidth available
                double currentMax = 0;
                int currentI = 0;
                int iMax = 0;
                auto itMax = pat.getDates().begin();
                //DEBUG(k << ": " << itBegin->date << " " << itEnd->date)
                for (auto it=itBegin; it!=itEnd; )
                {
                    //DEBUG("it is at " << it->date)
                    if (it->endPeriod)
                    {
                        it = pat.getDates().begin();
                    } else {
                        double availableBw = std::min(m.getBW() - it->totalBW,pat.getProcs(k)*m.getProcBW());
                        //DEBUG(availableBw)
                        if (availableBw > currentMax)
                        {
                            if ((availableBw < prevMax) || (availableBw == prevMax && currentI > prevI))
                            {
                                itMax = it;
                                currentMax = availableBw;
                                iMax = currentI;
                                //DEBUG("max chosen " << prevMax << " " << currentMax << " " << prevI << " " << currentI)
                            }
                        }
                        it++;
                        currentI++;
                    }
                }
                prevMax = currentMax;
                prevI = iMax;
                auto itMax2 = std::next(itMax);
                double maxAddedTime = myFmod(T+itMax2->date-itMax->date,T);
                if (pat.getDates().size()==2)
                    maxAddedTime = T;
                double TimeAdded = std::min(maxAddedTime,dataLeft/currentMax);
                //DEBUG(k << ": " << TimeAdded << " " << currentMax)
                Time += TimeAdded;
                dataLeft -= currentMax*TimeAdded;
                Interval IOint;
                IOint.ts = itMax->date;
                IOint.te = IOint.ts+TimeAdded;
                IOint.bw = currentMax;
                IOintList[prevI] = IOint;
                cpt++;
            }
			///NOT ENOUGH PLACE TO INSERT A NEW INSTANCE
            if (dataLeft > 0+PRECISION || (pat.getNbInstances(k) == 0 && Time > T-pat.getWork(k)+PRECISION))
                continue;
            ///Then mark this as the new best schedule
            if (Time < Tmin)
            {
                Tmin = Time;
                i_ref = i;
                FinalIOintList = IOintList;
            }
        }
    }
	///If we found a solution, return it
    if (Tmin < T+1)
    {
        for (auto it = FinalIOintList.begin(); it != FinalIOintList.end();)
        {
            if (it->bw == 0)
                it = FinalIOintList.erase(it);
            else
                it++;
        }
        return FinalIOintList;
    } else {
        return {};
    }
}

std::vector<Interval> splitSchedulableV2(Pattern& pat, const Machine& m, int k, std::chrono::microseconds& totalDuration)
{
    if (pat.getNbInstances(k) == 0)
    {
        return splitSchedulable(pat,m,k,totalDuration);
    } else {
    double T = pat.getPeriod();
    double Tmin = T+1;
    std::vector<Interval> IOintList;
    std::vector<Interval> FinalIOintList;
    auto itInstance = pat.getLastInstanceIterator(k);
    //for (int i=0; i<pat.getNbInstances(k) || (pat.getNbInstances(k)==0 && i==0); i++)
    //{

        //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //We compute the first moment where the new instance can send IO
        double dateBeginIo = 0;
        //And the last moment = beginning of next instance
        double dateEndIo = T;

        if (pat.getNbInstances(k) > 0)
        {
            dateBeginIo = pat.getEndIOTransferDate(k,itInstance);
            itInstance++;
            if (itInstance == pat.getInstances(k).end())
                itInstance = pat.getInstances(k).begin();
            dateEndIo = *itInstance;
        }
        bool placeToWork = false;
        if (dateEndIo > dateBeginIo+PRECISION)
        {
            if (dateEndIo >= dateBeginIo + pat.getWork(k) - PRECISION)
                placeToWork = true;
        } else if (dateEndIo+PRECISION < dateBeginIo) {
            if (dateEndIo+T >= dateBeginIo + pat.getWork(k) - PRECISION)
                placeToWork = true;
        }
		//Special case when there is no instance
        if (pat.getNbInstances(k) > 0)
            dateBeginIo = myFmod(dateBeginIo+pat.getWork(k),T);

        //std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        //std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );
        //totalDuration += duration;
        if (placeToWork || (pat.getNbInstances(k) == 0 && T >= pat.getWork(k)) )
        {
            //DEBUG(dateBeginIo << " " << dateEndIo)
            pat.addDate(dateBeginIo);//we need to split the interval
            auto itBegin = pat.getDates().find(dateElement(dateBeginIo));
            auto itEnd = pat.getDates().find(dateElement(dateEndIo));
			///Compute the bandwidths available
            ///ALREADY DONE NOW AHAHAH
            IOintList = std::vector<Interval>();
            double dataLeft = pat.getIOVol(k);
            double Time = 0;
            auto itMax = itBegin;
            while (dataLeft > 0+PRECISION && itMax!=itEnd)
            {
                double availableBw = std::min(m.getBW() - itMax->totalBW,pat.getProcs(k)*m.getProcBW());
                if (availableBw > 0)
                {
                    auto itMax2 = std::next(itMax);
                    double maxAddedTime = myFmod(T+itMax2->date-itMax->date,T);
                    if (pat.getDates().size()==2)
                        maxAddedTime = T;
                    double TimeAdded = std::min(maxAddedTime,dataLeft/availableBw);
                    //DEBUG(k << ": " << TimeAdded << " " << currentMax)
                    Time += TimeAdded;
                    dataLeft -= availableBw*TimeAdded;
                    Interval IOint;
                    IOint.ts = itMax->date;
                    IOint.te = IOint.ts+TimeAdded;
                    IOint.bw = availableBw;
                    IOintList.push_back(IOint);
                }
                itMax++;
                if (itMax == pat.getDates().end())
                    itMax = pat.getDates().begin();
            }
			///NOT ENOUGH PLACE TO INSERT A NEW INSTANCE
            if (dataLeft > 0+PRECISION || (pat.getNbInstances(k) == 0 && Time > T-pat.getWork(k)+PRECISION))
                Time = Tmin+1;
            ///Then mark this as the new best schedule
            if (Time < Tmin)
            {
                Tmin = Time;
                FinalIOintList = IOintList;
            }
        }
    //}
	///If we found a solution, return it
    if (Tmin < T+1)
    {
        for (auto it = FinalIOintList.begin(); it != FinalIOintList.end();)
        {
            if (it->bw == 0)
                it = FinalIOintList.erase(it);
            else
                it++;
        }
        return FinalIOintList;
    } else {
        return {};
    }
    }
}

void insertInSchedule(Pattern& pat, const Machine& m, std::function<std::vector<Interval>(Pattern&, const Machine&,int,std::chrono::microseconds&)> model, double epsilon)
{
    double T = pat.getPeriod();
    std::vector<bool> schedulable(pat.getNbApps(),true);
    int remainingApps = pat.getNbApps(); //To know if we are done
    /*for (int i=0; i < schedulable.size(); i++)
    {
        Pattern pat_copy = pat;
        if (model(pat_copy,m,i,).size() > 0)
            schedulable[i] = true;
        else
        {
            schedulable[i] = false;
            remainingApps--;
        }
    }*/
	///Compute W/TimeIO for each app
    std::vector<double> ratioWIO(pat.getNbApps());
    for (int i=0; i < ratioWIO.size(); i++)
        ratioWIO[i] = pat.getWork(i)*std::min(pat.getProcs(i)*m.getProcBW(),m.getBW())/pat.getIOVol(i);
	///Then insert while we can insert new instances of some apps
	int cpt=0;
    std::chrono::microseconds totalDuration(0);
    while (remainingApps > 0)
    {
        int nextApp = -1;
        double maxSd = 2;
        double maxRatio = -1;
		///We find the maximal app w.r.t (Dilation, W/TimeIO)
        for (int k=0; k<schedulable.size(); k++)
        {
            if (!schedulable[k])
                continue;
            double slowdown = pat.getNbInstances(k)*(pat.getWork(k)+pat.getIOVol(k)/std::min(m.getProcBW()*pat.getProcs(k),m.getBW()))/pat.getPeriod();
            if (slowdown < maxSd-PRECISION)
            {
                maxSd = slowdown;
                maxRatio = ratioWIO[k];
                nextApp = k;
            } else if (slowdown == maxSd && ratioWIO[k] > maxRatio+PRECISION) {
                maxRatio = ratioWIO[k];
                nextApp = k;
            }
        }
        assert(nextApp > -1);
        auto intv = model(pat,m,nextApp,totalDuration);
        //DEBUG(nextApp);
		///We test if nextApp is schedulable or not
        if (intv.size() > 0)
        {
            pat.addInstance(nextApp,intv);
        } else {
            schedulable[nextApp] = false;
            remainingApps--;
        }
        cpt++;
    }
    //std::cout << "Total time spent in finding dateBeginIo and dateEndIo : " <<  totalDuration.count() << "\n";
}
