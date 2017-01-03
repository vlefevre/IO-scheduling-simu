#include "Pattern.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <set>

#define DEBUG(X) std::cout << X << "\n"; std::cout.flush();

double myFmod(double nb, double mod)
{
    double res = fmod(nb,mod);
    if (res > mod-PRECISION && res < mod+PRECISION)
        res = 0;
    return res;
}

std::set<double>::iterator circular_next(std::set<double>& l, std::set<double>::iterator it, int n=1 )
{
    auto it2 = it;
    for (int i=0; i<n; i++)
    {
        it2++;
        if (it2 == l.end())
            it2 = l.begin();
    }
    return it2;
}

double Pattern::getBeginWDate(int app, int instance)
{
    if (app >= 0 && app < instances.size())
    {
        if (instance >= 0 && instance < instances[app].size())
        {
            auto it = instances[app].begin();
            std::advance(it,instance);
            return *it;
        }
    }
    return -1;
}

double Pattern::getEndWDate(int app, int instance)
{
    if (app >= 0 && app < instances.size())
    {
        if (instance >= 0 && instance < instances[app].size())
        {
            auto it = instances[app].begin();
            std::advance(it,instance);
            return myFmod(*it+apps[app].getWork(),getPeriod());
        }
    }
    return -1;
}

double Pattern::getBeginIODate(int app, int instance)
{
    return Pattern::getEndWDate(app, instance);
}

double Pattern::getEndIODate(int app, int instance)
{
    return Pattern::getBeginWDate(app, (instance+1)%instances[app].size());
}

double Pattern::getBeginIOTransferDate(int app, int instance)
{
    double beginIO = getBeginIODate(app, instance);
    if (beginIO > -1)
    {
        if (beginIO == T)
            beginIO = 0;
        auto iter = dates.find(dateElement(beginIO));
        while (iter->appBW[app] == 0)
        {
            iter++;
            if (iter->endPeriod)
            {
                iter = dates.begin();
            }
        }
        return iter->date;
    }
    return -1;
}

double Pattern::getEndIOTransferDate(int app, int instance)
{
    double endIO = getEndIODate(app, instance);
    if (endIO > -1)
    {
        if (endIO == 0)
            endIO += T;
        auto iter = dates.find(dateElement(endIO));
        while (iter->appBW[app] == 0)
        {
            if (iter == dates.begin())
            {
                iter = dates.end();
            }
            iter--;
        }
        return std::next(iter)->date;
    }
    return -1;
}

double Pattern::getBeginIOTransferDate(int app, std::set<double>::iterator it)
{
    double beginIO = *it;
    if (beginIO > -1)
    {
        if (beginIO == T)
            beginIO = 0;
        auto iter = dates.find(dateElement(beginIO));
        while (iter->appBW[app] == 0)
        {
            iter++;
            if (iter->endPeriod)
            {
                iter = dates.begin();
            }
        }
        return iter->date;
    }
    return -1;
}

double Pattern::getEndIOTransferDate(int app, std::set<double>::iterator it)
{
    auto nextIt = std::next(it);
    if (nextIt == instances[app].end())
        nextIt = instances[app].begin();
    double endIO = *nextIt;
    if (endIO > -1)
    {
        if (endIO == 0)
            endIO += T;
        auto iter = dates.find(dateElement(endIO));
        while (iter->appBW[app] == 0)
        {
            if (iter == dates.begin())
            {
                iter = dates.end();
            }
            iter--;
        }
        return std::next(iter)->date;
    }
    return -1;
}

int Pattern::getNbApps()
{
    return apps.size();
}

int Pattern::getNbInstances(int app)
{
    return instances[app].size();
}

double Pattern::getWork(int app)
{
    return apps[app].getWork();
}

double Pattern::getIOVol(int app)
{
    return apps[app].getIOVol();
}

double Pattern::getProcs(int app)
{
    return apps[app].getProcs();
}

void Pattern::addDate(double newDate)
{
    newDate = myFmod(newDate+T,T);
    auto returnValue = dates.insert(dateElement(newDate));
    if (!returnValue.second)
        return;
    auto iter = returnValue.first;
    auto iter2 = iter; //Element inserted
    iter--; //Element before
    iter2->totalBW = iter->totalBW;
    for (int i=0; i<iter->appBW.size(); i++)
    {
        iter2->appBW.push_back(iter->appBW[i]);
    }
}

void Pattern::addInstanceDate(int app, double date)
{
    date = myFmod(date+T,T);
    if (app >= 0 && app < apps.size())
    {
        addDate(date);
        instances[app].insert(date);
    }
}

void Pattern::addInstance(int app, const std::vector<Interval> &intList)
{
    addInstanceDate(app, intList[0].ts-apps[app].getWork());
    for (auto i:intList)
    {
        addDate(i.ts);
        addDate(i.te);
    }
    auto iter = dates.begin();
    auto iterInt = intList.begin();
    int cpt = 0;
    do
    {
        double ts = myFmod(iterInt->ts+T,T);
        double te = myFmod(iterInt->te+T,T);
        iter = dates.find(dateElement(ts));
        while (iter->date < te-PRECISION || iter->date > te+PRECISION)
        {
            iter->appBW[app] = iterInt->bw;
            iter->totalBW += iterInt->bw;
            iter++;
            if (iter->endPeriod)
            {
                    iter = dates.begin();
            }
        }
        iterInt++;
        cpt++;
    } while (iterInt != intList.end());
    assert(cpt == intList.size());
}

/*void Pattern::setPeriod(double newT)
{
}*/

double Pattern::getBWUsed(double t)
{
    t = myFmod(t,T);
    auto it = dates.find(dateElement(t));
    assert(it != dates.end());
    return it->totalBW;
}

void Pattern::print()
{
    std::cout << std::fixed << std::setprecision(3)<< "Dates :\t";
    for (auto d:dates)
        std::cout << std::setfill('0') << std::setw(9) << d.date << " \t";
    std::cout << "\nTotal :\t";
    for (auto d:dates)
        std::cout << std::setfill('0') << std::setw(9) << d.totalBW << " \t";
    std::cout << "\n";
    for (int k=0; k<getNbApps(); k++)
    {
        std::cout << "App " << k << " :\t";
        for (auto d:dates)
        {
            std::cout << std::setfill('0') << std::setw(9) << d.appBW[k] << " \t";
        }
        std::cout << "\n";
    }
    for (int k=0; k<apps.size(); k++)
    {
        std::cout << "Application " << k << " (" << getNbInstances(k) << "):\n";
        //auto iter = dates.begin();
        std::cout << "Instances :\n";
        for (auto d:instances[k])
            std::cout << d << " ";
        std::cout << "\n";
        for (auto i=instances[k].begin(); i != instances[k].end(); i++)
        {
            auto index = std::distance(instances[k].begin(),i);
            std::cout << "Works between " << getBeginWDate(k,index) << " and " << getEndWDate(k,index) << ".\n";
            std::cout << "Performs I/O between " << getBeginIOTransferDate(k,index) << " and " << getEndIOTransferDate(k,index) << ".\n";
        }
        std::cout << "\n";
    }
}

///For CPLEX optimization

void Pattern::print_cplex(std::string filename,double B,double b,std::vector<std::vector<int> >& priority,std::string inputLP = "",std::vector<std::pair<int,int> > inverts = std::vector<std::pair<int,int> >())
{
    std::ofstream out(filename,std::ios::out);
    std::vector<dateApp> events;
    int prioSize = 0;
    for(int i = 0; i < priority.size(); i++)
        prioSize += priority[i].size();
    //DEBUG("print cplex");
    //if (priority[0].size() > 3)
    //    DEBUG("priority " << priority[0][3])
    if (inputLP == "")
    {
        events = std::vector<dateApp>();
        for (auto e:eventsSet)
        {
            events.push_back(dateApp(e.date,e.app,e.instance,e.start));
        }
    } else {
        std::ifstream in(inputLP,std::ios::in);
        std::string buff;
        while (buff != "END")
            in >> buff;
        in.ignore();
        in.ignore();
        in.ignore();//Go to the line of events
        double readDate;
        int readApp;
        int readInstance;
        bool readStart;
        while (in >> readDate >> readApp >> readInstance >> readStart)
        {
            events.push_back(dateApp(readDate,readApp,readInstance,readStart));
        }
        for (auto p:inverts)
        {
            if (((p.first-1+p.second)%events.size() < p.first-1 || p.first==1)  && events[p.first-1].start && (p.first-1+p.second)%events.size() > 0)
            {
                priority[events[p.first-1].app][(events[p.first-1].instance-1+getNbInstances(events[p.first-1].app))%getNbInstances(events[p.first-1].app)] -= prioSize;
                //DEBUG("decreased priority of " << events[p.first-1].app << " " << events[p.first-1].instance-1)
            }
            dateApp temp = events[p.first-1];
            for (auto i=0; i<p.second; i++)
            {
                if ((p.first+i)%events.size() == 1 && events[(p.first+i)%events.size()].start){
                    priority[events[(p.first+i)%events.size()].app][(events[(p.first+i)%events.size()].instance-1+getNbInstances(events[(p.first+i)%events.size()].app))%getNbInstances(events[(p.first+i)%events.size()].app)] += prioSize;
                    //DEBUG("increased priority of " << events[(p.first+i)%events.size()].app << " " << events[(p.first+i)%events.size()].instance-1)
                }
                events[(p.first-1+i)%events.size()] = events[(p.first+i)%events.size()];
            }
            events[(p.first-1+p.second)%events.size()] = temp;
        }
    }

    //DEBUG("interversions done");
    /*DEBUG("EVENTS :")
    for (auto e:events)
    { DEBUG(e.date << " " << e.app << " " << e.instance << " " << e.start) }
    DEBUG("EVENTS SET :")
    for (auto e:eventsSet)
    { DEBUG(e.date << " " << e.app << " " << e.start) }*/

    out << "MINIMIZE\n\n";
    for (int i=0; i<events.size(); i++)
    {
        out << "t_" << (i+1);
        if (i < events.size()-1)
            out << " + ";
        else
            out << "\n";
    }
    out << "\nSUBJECT TO\n\n";
    std::vector<std::vector<double> > lpBW(getNbApps());

    std::set<PrioInstance,PriorityComp> IOapps;
    //First find the events of beginning of IO for all apps sending at the beginning
    for (auto k=0; k<getNbApps(); k++)
    {
        auto it = events.begin();
        while (true)
        {
            if (it->app == k)
            {
                if (!it->start)
                    IOapps.insert(PrioInstance(k,it->instance,priority[k][it->instance]));
                break;
            }
            if (it == events.begin())
                it = events.end();
            it--;
        }
    }

    //DEBUG("first applications found");
    std::vector<std::vector<bool> > seen(getNbApps());
    for (int k=0; k<getNbApps(); k++)
        seen[k] = std::vector<bool>(getNbInstances(k),false);
    int num_interval = 1;
    for (auto it=events.begin(); it!=events.end(); it++)
    {
        double availableBW = B;
        for (int k=0; k<getNbApps(); k++)
            lpBW[k].push_back(0);
        if (it != events.begin())
        {if (it->start) //one app start working so we loose one instance (app) in IOapps
        {
            //DEBUG("Removing " << it->app << " " << it->instance-1 << " " << priority[it->app][(it->instance-1+getNbInstances(it->app))%getNbInstances(it->app)])
            IOapps.erase(PrioInstance(0,0,priority[it->app][(it->instance-1+getNbInstances(it->app))%getNbInstances(it->app)]));
            seen[it->app][(it->instance-1+getNbInstances(it->app))%getNbInstances(it->app)] = true;
        } else { //one app start sending I/O
            //if (seen[it->app][it->instance])
            //    DEBUG("Already seen " << it->app << " " << it->instance)
            //DEBUG("Adding " << it->app << " " << it->instance << " " << priority[it->app][it->instance])
            if (seen[it->app][it->instance])
                IOapps.insert(PrioInstance(it->app,it->instance,priority[it->app][it->instance]+prioSize));
            else
                IOapps.insert(PrioInstance(it->app,it->instance,priority[it->app][it->instance]));
        }}

        /*if (num_interval == 2)
        {
            DEBUG("t_" << num_interval)
            for (auto rgi:IOapps)
                DEBUG(rgi.app << " " << rgi.instance << " " << rgi.priority);
        }*/
        for (auto itIO = IOapps.begin(); itIO != IOapps.end(); itIO++)
        {
            lpBW[itIO->app][lpBW[itIO->app].size()-1] = std::min(availableBW,b*getProcs(itIO->app));
            availableBW -= lpBW[itIO->app].back();
        }
        num_interval++;
    }
    //DEBUG("bw remplies")
    for (int k=0; k<getNbApps(); k++)
    {
        auto it = events.begin();
        while (it->app != k || !it->start) it++;
        auto it2 = it;
        bool working = true;
        do
        {
            int index=std::distance(events.begin(),it2);
            if (working)
            {
                out << "t_" << (index+1);
            } else {
                out << lpBW[k][index] << " t_" << (index+1);
            }
            it2++;
            if (it2 == events.end()) it2 = events.begin();
            if (working && !it2->start && it2->app == k) //end of working phase
            {
                out << " >= " << getWork(k) << "\n";
                working = false;
            } else if (!working && it2->start && it2->app == k) { //end of io phase
                out << " >= " << getIOVol(k) << "\n";
                working = true;
            } else {
                out << " + ";
            }
        } while (it2 != it);
        out << "\n";
    }
    //DEBUG("constraints written")
    out << "BOUNDS\n\n";
    for (int i=0; i<events.size(); i++)
    {
        out << "t_" << (i+1) << " >= 0\n";
    }
    out << "\nEND\n\n\n";

    for (dateApp e:events)
    {
        out << std::setprecision(10) << e.date << " " << e.app << " " << e.instance << " " << e.start << std::setprecision(3) << "\n";
    }
    out.close();
}

std::vector<std::vector<int> > Pattern::init_cplex()
{
    std::vector<std::vector<int> > priority;
    for (int k=0; k<getNbApps(); k++)
    {
        priority.push_back(std::vector<int>(getNbInstances(k),0));
        for (int i=0; i<getNbInstances(k); i++)
        {
            //workingTimes[k].push_back(std::make_pair(getBeginWDate(k,i),getEndWDate(k,i)));
            eventsSet.insert(dateApp(getBeginWDate(k,i),k,i,true));
            eventsSet.insert(dateApp(getEndWDate(k,i),k,i,false));
        }
    }
    int prio = 0;
    for (auto e:eventsSet)
    {
        if (e.start)
        {
            priority[e.app][(e.instance-1+getNbInstances(e.app))%getNbInstances(e.app)] = prio;
            prio++;
        }
    }

    return priority;
}

void Pattern::exportSchedule(std::string filename)
{
    std::ofstream out(filename, std::ios::out);
    for (int k=0; k<getNbApps(); k++)
    {
        out << "[ ";
        for (int i=0; i<getNbInstances(k); i++)
        {
            out << "[ ";
            double dateBegin = getBeginIOTransferDate(k,i);
            auto it = dates.find(dateElement(dateBegin));
            double dateEnd = getEndIOTransferDate(k,i);
            auto itEnd = dates.find(dateElement(dateEnd));
            bool first = true;
            while (it != itEnd)
            {
                if (it->appBW[k] > 0+PRECISION)
                {
                    if (!first)
                        out << ", ";
                    out << "( " << it->date << " , " << std::next(it)->date << " , " << it->appBW[k] << " )";
                    first = false; 
                }
                it++;
                if (it->endPeriod)
                    it = dates.begin();
            }
            out << " ] ";
        }
        out << "]\n";
    }
}
