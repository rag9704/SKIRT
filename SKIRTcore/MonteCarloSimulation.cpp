/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustDistribution.hpp"
#include "DustMix.hpp"
#include "DustSystem.hpp"
#include "FatalError.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "Log.hpp"
#include "MonteCarloSimulation.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "ProcessAssigner.hpp"
#include "Random.hpp"
#include "StellarSystem.hpp"
#include "TextOutFile.hpp"
#include "TimeLogger.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"
// The following includes are used for dynamic casting and should not be included in the final build
#include "VoronoiDustGrid.hpp"
#include "TreeDustGrid.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////

MonteCarloSimulation::MonteCarloSimulation()
    : _is(0), _packages(0), _minWeightReduction(1e4),
      _minfs(0), _xi(0.5), _continuousScattering(false),
      _lambdagrid(0), _ss(0), _ds(0)
{
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();

    // protect implementation limit
    if (_packages < 0)
        throw FATALERROR("Number of photon packages is negative");
    if (_packages > 1e15)
        throw FATALERROR("Number of photon packages is larger than implementation limit of 1e15");
    if (_minWeightReduction < 1e3)
        throw FATALERROR("The minimum weight reduction factor should be larger than 1000");
    if (_minfs < 0)
        throw FATALERROR("The minimum number of scattering events is negative");
    if (_minfs > 1000)
        throw FATALERROR("The minimum number of scattering events should be smaller than 1000");
    if (_xi < 0 || _xi > 1)
        throw FATALERROR("The scattering bias should be between 0 and 1");
    if (!_lambdagrid)
        throw FATALERROR("Wavelength grid was not set");
    if (!_ss)
        throw FATALERROR("Stellar system was not set");
    if (!_is)
        throw FATALERROR("Instrument system was not set");
    if (_prePackages < 0)
        throw FATALERROR("Number of prepackages is negative");
    if (_dynamicIterations < 0)
        throw FATALERROR("Number of dynamic iterations is negative");
    if ((_prePackages != 0) != (_dynamicIterations!=0))
        throw FATALERROR("All dynamic grid parameters should be either on or off.");
    // dust system is optional; nr of packages has a valid default
    // Cache the total amount of packages
    _totalPackages = _packages;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setChunkParams(double packages)
{
    // Cache the number of wavelengths
    _Nlambda = _lambdagrid->Nlambda();

    // Determine the number of chunks and the corresponding chunk size
    if (packages <= 0)
    {
        _Nchunks = 0;
        _chunksize = 0;
        _Npp = 0;
    }
    else
    {
        // Get the number of processes and threads per process
        int Nprocs = _comm->size();
        int Nthreads = _parfac->maxThreadCount();

        // Step 1: consider threading and determine the total number of chunks
        int totalChunks = 0;
        if (Nthreads == 1) totalChunks = 1;
        else totalChunks = ceil( std::max({(10.*double(Nthreads*Nprocs)/_Nlambda), packages/1e7}) );

        // Step 2: consider the work division and determine the number of chunks per process (_Nchunks)
        if (_comm->dataParallel())  // Do some wavelengths for all chunks
        {
            _chunksize = ceil(packages/totalChunks);
            _Nchunks = totalChunks;
            _myTotalNpp = _lambdagrid->assigner()->assigned() * _Nchunks * _chunksize;
        }
        else                        // Do all wavelengths for some chunks
        {
            if ((totalChunks % Nprocs)) totalChunks = totalChunks + Nprocs - (totalChunks % Nprocs);
            _chunksize = ceil(packages/totalChunks);
            _Nchunks = totalChunks/Nprocs;
            _myTotalNpp = _Nlambda * _Nchunks * _chunksize;
        }

        // Calculate the the definitive number of photon packages per wavelength
        _Npp = totalChunks * _chunksize;

        _log->info("Using " + QString::number(totalChunks) + " chunks per wavelength");
    }

    // Determine the log frequency; continuous scattering is much slower!
    _logchunksize = _continuousScattering ? 5000 : 50000;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setInstrumentSystem(InstrumentSystem* value)
{
    if (_is) delete _is;
    _is = value;
    if (_is) _is->setParent(this);
}

////////////////////////////////////////////////////////////////////

InstrumentSystem* MonteCarloSimulation::instrumentSystem() const
{
    return _is;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setPackages(double value)
{
    _packages = value;
}

////////////////////////////////////////////////////////////////////

double MonteCarloSimulation::packages() const
{
    return _packages;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setMinWeightReduction(double value)
{
    _minWeightReduction = value;
}

////////////////////////////////////////////////////////////////////

double MonteCarloSimulation::minWeightReduction() const
{
    return _minWeightReduction;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setMinScattEvents(double value)
{
    _minfs = value;
}

////////////////////////////////////////////////////////////////////

double MonteCarloSimulation::minScattEvents() const
{
    return _minfs;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setScattBias(double value)
{
    _xi = value;
}

////////////////////////////////////////////////////////////////////

double MonteCarloSimulation::scattBias() const
{
    return _xi;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setContinuousScattering(bool value)
{
    _continuousScattering = value;
}

////////////////////////////////////////////////////////////////////

bool MonteCarloSimulation::continuousScattering() const
{
    return _continuousScattering;
}

//////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setPrePackages(double value)
{
    _prePackages = value;
}

//////////////////////////////////////////////////////////////////////

double MonteCarloSimulation::prePackages() const
{
    return _prePackages;
}

//////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setDynamicIterations(int value)
{
    _dynamicIterations = value;
}

//////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dynamicIterations() const
{
    return _dynamicIterations;
}

////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dimension() const
{
    return qMax(_ss->dimension(), _ds ? _ds->dimension() : 1);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::initprogress(QString phase)
{
    _phase = phase;
    _Ndone = 0;

    _log->info(QString::number(_Npp) + " photon packages for "
               + (_Nlambda==1 ? QString("a single wavelength") : QString("each of %1 wavelengths").arg(_Nlambda)));

    _timer.start();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logprogress(quint64 extraDone)
{
    // accumulate the work already done
    _Ndone.fetch_add(extraDone);

    // space the messages at least 3 seconds apart; in the interest of speed,
    // we do this without locking, so once in a while two consecutive messages may slip through
    if (_timer.elapsed() > 3000)
    {
        _timer.restart();
        double completed = _Ndone * 100. / (_myTotalNpp);
        _log->info("Launched " + _phase + " photon packages: " + QString::number(completed,'f',1) + "%");
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runstellaremission()
{
    TimeLogger logger(_log, "the stellar emission phase");
    setChunkParams(_packages);
    initprogress("stellar emission");
    Parallel* parallel = find<ParallelFactory>()->parallel();

    if (_lambdagrid->assigner())
        parallel->call(this, &MonteCarloSimulation::dostellaremissionchunk, _lambdagrid->assigner(), _Nchunks);
    else
        parallel->call(this, &MonteCarloSimulation::dostellaremissionchunk, _Nlambda, _Nchunks);

    // Wait for the other processes to reach this point
    _comm->wait("the stellar emission phase");
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::dostellaremissionchunk(size_t index)
{
    int ell = index % _Nlambda;
    double L = _ss->luminosity(ell)/_Npp;
    if (L > 0)
    {
        double Lthreshold = L / minWeightReduction();
        PhotonPackage pp,ppp;

        quint64 remaining = _chunksize;
        while (remaining > 0)
        {
            quint64 count = qMin(remaining, _logchunksize);
            for (quint64 i=0; i<count; i++)
            {
                _ss->launch(&pp,ell,L);
                if (pp.luminosity()>0)
                {
                    peeloffemission(&pp,&ppp);
                    if (_ds) while (true)
                    {
                        _ds->fillOpticalDepth(&pp);
                        if (_continuousScattering) continuouspeeloffscattering(&pp,&ppp);
                        simulateescapeandabsorption(&pp,_ds->storeabsorptionrates());
                        if (pp.luminosity()<=0 || (pp.luminosity()<=Lthreshold && pp.nScatt()>=_minfs)) break;
                        simulatepropagation(&pp);
                        if (!_continuousScattering) peeloffscattering(&pp,&ppp);
                        simulatescattering(&pp);
                    }
                }
            }
            logprogress(count);
            remaining -= count;
        }
    }
    else logprogress(_chunksize);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peeloffemission(const PhotonPackage* pp, PhotonPackage* ppp)
{
    Position bfr = pp->position();

    foreach (Instrument* instr, _is->instruments())
    {
        Direction bfknew = instr->bfkobs(bfr);
        ppp->launchEmissionPeelOff(pp, bfknew);
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::peeloffscattering(const PhotonPackage* pp, PhotonPackage* ppp)
{
    int Ncomp = _ds->Ncomp();
    int ell = pp->ell();
    Position bfr = pp->position();

    // Determine the weighting factors of the phase functions corresponding to
    // the different dust components: each component h is weighted by kappasca(h)*rho(m,h)
    QVarLengthArray<double,4> wv(Ncomp);
    if (Ncomp==1)
        wv[0] = 1.0;
    else
    {
        int m = _ds->whichcell(bfr);
        if (m==-1) return; // abort peel-off
        for (int h=0; h<Ncomp; h++) wv[h] = _ds->mix(h)->kappasca(ell) * _ds->density(m,h);
        double sum = 0;
        for (int h=0; h<Ncomp; h++) sum += wv[h];
        if (sum<=0) return; // abort peel-off
        for (int h=0; h<Ncomp; h++) wv[h] /= sum;
    }

    // Now do the actual peel-off
    foreach (Instrument* instr, _is->instruments())
    {
        Direction bfkobs = instr->bfkobs(bfr);
        Direction bfkx = instr->bfkx();
        Direction bfky = instr->bfky();
        double I = 0, Q = 0, U = 0, V = 0;
        for (int h=0; h<Ncomp; h++)
        {
            DustMix* mix = _ds->mix(h);
            double w = wv[h] * mix->phaseFunctionValue(pp, bfkobs);
            StokesVector sv;
            mix->scatteringPeelOffPolarization(&sv, pp, bfkobs, bfkx, bfky);
            I += w * sv.stokesI();
            Q += w * sv.stokesQ();
            U += w * sv.stokesU();
            V += w * sv.stokesV();
        }
        ppp->launchScatteringPeelOff(pp, bfkobs, I);
        ppp->setPolarized(I, Q, U, V, pp->normal());
        instr->detect(ppp);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::continuouspeeloffscattering(const PhotonPackage *pp, PhotonPackage *ppp)
{
    int ell = pp->ell();
    Position bfr = pp->position();
    Direction bfk = pp->direction();

    int Ncomp = _ds->Ncomp();
    QVarLengthArray<double,4> kappascav(Ncomp);
    QVarLengthArray<double,4> kappaextv(Ncomp);
    for (int h=0; h<Ncomp; h++)
    {
        DustMix* mix = _ds->mix(h);
        kappascav[h] = mix->kappasca(ell);
        kappaextv[h] = mix->kappaext(ell);
    }

    int Ncells = pp->size();
    for (int n=0; n<Ncells; n++)
    {
        int m = pp->m(n);
        if (m!=-1)
        {
            QVarLengthArray<double,4> wv(Ncomp);
            double ksca = 0.0;
            double kext = 0.0;
            for (int h=0; h<Ncomp; h++)
            {
                double rho = _ds->density(m,h);
                wv[h] = rho*kappascav[h];
                ksca += rho*kappascav[h];
                kext += rho*kappaextv[h];
            }
            if (ksca>0.0)
            {
                for (int h=0; h<Ncomp; h++) wv[h] /= ksca;
                double albedo = ksca/kext;
                double tau0 = (n==0) ? 0.0 : pp->tau(n-1);
                double dtau = pp->dtau(n);
                double s0 = (n==0) ? 0.0 : pp->s(n-1);
                double ds = pp->ds(n);
                double factorm = albedo * exp(-tau0) * (-expm1(-dtau));
                double s = s0 + _random->uniform()*ds;
                Position bfrnew(bfr+s*bfk);
                foreach (Instrument* instr, _is->instruments())
                {
                    Direction bfkobs = instr->bfkobs(bfrnew);
                    Direction bfkx = instr->bfkx();
                    Direction bfky = instr->bfky();
                    double I = 0, Q = 0, U = 0, V = 0;
                    for (int h=0; h<Ncomp; h++)
                    {
                        DustMix* mix = _ds->mix(h);
                        double w = wv[h] * mix->phaseFunctionValue(pp, bfkobs);
                        StokesVector sv;
                        mix->scatteringPeelOffPolarization(&sv, pp, bfkobs, bfkx, bfky);
                        I += w * sv.stokesI();
                        Q += w * sv.stokesQ();
                        U += w * sv.stokesU();
                        V += w * sv.stokesV();
                    }
                    ppp->launchScatteringPeelOff(pp, bfrnew, bfkobs, factorm*I);
                    ppp->setPolarized(I, Q, U, V, pp->normal());
                    instr->detect(ppp);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulateescapeandabsorption(PhotonPackage* pp, bool storeabsorptionrates)
{
    double taupath = pp->tau();
    int ell = pp->ell();
    double L = pp->luminosity();
    bool ynstellar = pp->isStellar();
    int Ncomp = _ds->Ncomp();

    // Easy case: there is only one dust component
    if (Ncomp==1)
    {
        double albedo = _ds->mix(0)->albedo(ell);
        double expfactor = -expm1(-taupath);
        if (storeabsorptionrates)
        {
            int Ncells = pp->size();
            for (int n=0; n<Ncells; n++)
            {
                int m = pp->m(n);
                if (m!=-1)
                {
                    double taustart = (n==0) ? 0.0 : pp->tau(n-1);
                    double dtau = pp->dtau(n);
                    double expfactorm = -expm1(-dtau);
                    double Lintm = L * exp(-taustart) * expfactorm;
                    double Labsm = (1.0-albedo) * Lintm;
                    _ds->absorb(m,ell,Labsm,ynstellar);
                }
            }
        }
        double Lsca = L * albedo * expfactor;
        pp->setLuminosity(Lsca);
    }

    // Difficult case: there are different dust components.
    // The absorption/scattering in each cell is weighted by the density contribution of the component.
    else
    {
        Array kappascav(Ncomp);
        Array kappaextv(Ncomp);
        for (int h=0; h<Ncomp; h++)
        {
            DustMix* mix = _ds->mix(h);
            kappascav[h] = mix->kappasca(ell);
            kappaextv[h] = mix->kappaext(ell);
        }
        int Ncells = pp->size();
        double Lsca = 0.0;
        for (int n=0; n<Ncells; n++)
        {
            int m = pp->m(n);
            if (m!=-1)
            {
                double ksca = 0.0;
                double kext = 0.0;
                for (int h=0; h<Ncomp; h++)
                {
                    double rho = _ds->density(m,h);
                    ksca += rho*kappascav[h];
                    kext += rho*kappaextv[h];
                }
                double albedo = (kext>0.0) ? ksca/kext : 0.0;
                double taustart = (n==0) ? 0.0 : pp->tau(n-1);
                double dtau = pp->dtau(n);
                double expfactorm = -expm1(-dtau);
                double Lintm = L * exp(-taustart) * expfactorm;
                double Lscam = albedo * Lintm;
                Lsca += Lscam;
                if (storeabsorptionrates)
                {
                    double Labsm = (1.0-albedo) * Lintm;
                    _ds->absorb(m,ell,Labsm,ynstellar);
                }
            }
        }
        pp->setLuminosity(Lsca);
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulatepropagation(PhotonPackage* pp)
{
    double taupath = pp->tau();
    if (taupath==0.0) return;
    double tau = 0.0;
    if (_xi==0.0)
        tau = _random->exponcutoff(taupath);
    else
    {
        double X = _random->uniform();
        tau = (X<_xi) ? _random->uniform()*taupath : _random->exponcutoff(taupath);
        double p = -exp(-tau)/expm1(-taupath);
        double q = (1.0-_xi)*p + _xi/taupath;
        double weight = p/q;
        pp->setLuminosity(pp->luminosity()*weight);
    }
    double s = pp->pathlength(tau);
    pp->propagate(s);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::simulatescattering(PhotonPackage* pp)
{
    // Randomly select a dust mix; the probability of each dust component h is weighted by kappasca(h)*rho(m,h)
    DustMix* mix = _ds->randomMixForPosition(pp->position(), pp->ell());

    // Now perform the scattering using this dust mix
    Direction bfknew = mix->scatteringDirectionAndPolarization(pp, pp);
    pp->scatter(bfknew);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::dynamicGrid()
{
    Log* log = find<Log>();
    // Check if we have a voronoi dustgrid
    if(dynamic_cast<VoronoiDustGrid*>(_ds->dustGrid()) != NULL)
    {
        // Draw extra voronoi points from temperature distribution
        VoronoiDustGrid* vorogrid = dynamic_cast<VoronoiDustGrid*>(_ds->dustGrid());
        vorogrid->drawFromTemperatureDistribution();
        // Make sure all DustSystem variables are reinitialized correctly (eg the absorption table)
        _ds->reinitialiseGrid();
        log->info("Reinitialized voronoi grid");
    }
    else if(dynamic_cast<TreeDustGrid*>(_ds->dustGrid()) != NULL)
    {
        TreeDustGrid* treeGrid = dynamic_cast<TreeDustGrid*>(_ds->dustGrid());
        treeGrid->dynamicGrid();
        // Make sure all DustSystem variables are reinitialized correctly (eg the absorption table)
        _ds->reinitialiseGrid();
        log->info("Reinitialized tree grid");
    }
    else
    {
        throw FATALERROR("Use a voronoi grid or tree dust grid when using a dynamic grid!");
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::write(QString filenameSuffix)
{
    TimeLogger logger(_log, "writing results");
    if (_is && filenameSuffix=="") _is->write(); // Only write out results at end of simulation (no suffix)
    if (_ds) _ds->write(filenameSuffix);
}

////////////////////////////////////////////////////////////////////
