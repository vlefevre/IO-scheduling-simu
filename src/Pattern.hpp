#ifndef PATTERN_HPP
#define PATTERN_HPP

#include <vector>
#include <set>
#include <list>
#include <string>

#define PRECISION 0.00002

struct dateElement
{
    dateElement(double d) : date(d) { totalBW = 0; appBW = std::vector<double>(); endPeriod = false; }
    double date;
    mutable double totalBW;
    mutable std::vector<double> appBW;
    mutable bool endPeriod;
};

struct dateComparison {
    bool operator()(const dateElement& lhs, const dateElement& rhs) const { return (lhs.date+PRECISION<rhs.date); }
};

struct dateApp
{
    double date;
    int app;
    int instance;
    bool start;
    dateApp(double d, int a, int i, bool s) : date(d), app(a), instance(i), start(s) {}

    bool operator==(const dateApp& d) const
    {
        return (instance == d.instance && app==d.app && start==d.start);
    }

    bool operator<(const dateApp& d) const
    {
        if (date<d.date-PRECISION) return true;
        if (date>=d.date-PRECISION && date <= d.date+PRECISION && app<d.app) return true;
        return false;
    }
};

struct PrioInstance
{
    int app;
    int instance;
    int priority;
    PrioInstance(int k,int i, int p) : app(k), instance(i), priority(p) {}
};

struct PriorityComp {
    bool operator()(const PrioInstance& lhs, const PrioInstance& rhs) const { return (lhs.priority<rhs.priority); }
};

struct Interval
{
    double ts; //beginning time of interval
    double bw; //bandwidth during the interval
    double te; //ending time of interval
};

class InstanceList
{
    public:
        std::set<double> instances;
        std::set<double>::iterator getLastInstanceIterator() { return lastInstance; }
        void insert(double instance) { lastInstance = instances.insert(instance).first; }
        int size() { return instances.size(); }
        std::set<double>::iterator begin() { return instances.begin(); }
        std::set<double>::iterator end() { return instances.end(); }
    private:
        std::set<double>::iterator lastInstance;
};

class Application
{
public:
    Application(double Work, double IO, int Procs) : W(Work), IOVol(IO), NbProcs(Procs) {}
    double getWork() { return W; }
    double getIOVol() { return IOVol; }
    int getProcs() { return NbProcs; }
private:
    double W; //work
    double IOVol; //number of bytes exchanged during I/O transfers
    int NbProcs;
};

class Pattern
{
public:
    Pattern(double T, const std::vector<Application> &apps) : T(T), apps(apps) {
        dateElement d(0);
        d.appBW = std::vector<double>(apps.size(),0);
        dates.insert(d);
        d.date = T;
        d.endPeriod = true;
        dates.insert(d);
        instances = std::vector<InstanceList >(apps.size());
    }
    double getPeriod() { return T; }
    //void setPeriod(double newT);
    double getBeginWDate(int app, int instance);
    double getEndWDate(int app, int instance);
    double getBeginIODate(int app, int instance);
    double getEndIODate(int app, int instance);
    double getBeginIOTransferDate(int app, int instance);
    double getEndIOTransferDate(int app, int instance);
    double getBeginIOTransferDate(int app, std::set<double>::iterator);
    double getEndIOTransferDate(int app, std::set<double>::iterator);
    int getNbApps();
    int getNbInstances(int app);
    double getWork(int app);
    double getIOVol(int app);
    double getProcs(int app);
    std::set<dateElement,dateComparison>& getDates() { return dates; }
    std::set<double>& getInstances(int k) { return instances[k].instances; }
    std::set<double>::iterator getLastInstanceIterator(int k) { return instances[k].getLastInstanceIterator(); }
    void addDate(double newDate); //add a new date in the timeline
    void addInstance(int app, const std::vector<Interval>& listInt); //add an instance of application app with I/O sent according to the intervals listInt
    double getBWUsed(double t); //returns the bandwidth used at instant t
    void print(); //prints detailed information about the pattern
    void exportSchedule(std::string filename);
    void clear(double newT) { *this = Pattern(newT,apps); }
    std::vector<std::vector<int> > init_cplex();
    void print_cplex(std::string filename,double B,double b,std::vector<std::vector<int> >& ,std::string inputLP,std::vector<std::pair<int,int> >); //writes a LP which minimize the period of the pattern with given bandwidths
private:
    double T;
    std::vector<Application> apps;
    std::set<dateElement,dateComparison> dates; //instants during the period, with corresponding bandwidths
    std::vector<InstanceList> instances; //for each app, the beginning of the working periods i.e. the beginning of the instances

    std::set<dateApp> eventsSet;
    void addInstanceDate(int app, double date); //add a special date defining the beginning of an instance
};

#endif
