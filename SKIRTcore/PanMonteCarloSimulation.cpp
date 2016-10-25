/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Log.hpp"
#include "NR.hpp"
#include "PanDustSystem.hpp"
#include "PanMonteCarloSimulation.hpp"
#include "PanWavelengthGrid.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "PhotonPackage.hpp"
#include "Random.hpp"
#include "SED.hpp"
#include "StellarSystem.hpp"
#include "TimeLogger.hpp"
#include "Units.hpp"
#include "FatalError.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////

PanMonteCarloSimulation::PanMonteCarloSimulation()
    : _pds(0)
{
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::setupSelfAfter()
{
    MonteCarloSimulation::setupSelfAfter();

    // properly size the array used to communicate between rundustXXX() and the corresponding parallel loop
    _Ncells = _pds ? _pds->Ncells() : 0;
    if (_pds && _pds->dustemission()) _Labsbolv.resize(_Ncells);
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::setWavelengthGrid(PanWavelengthGrid* value)
{
    if (_lambdagrid) delete _lambdagrid;
    _lambdagrid = value;
    if (_lambdagrid) _lambdagrid->setParent(this);
}

////////////////////////////////////////////////////////////////////

PanWavelengthGrid* PanMonteCarloSimulation::wavelengthGrid() const
{
    return dynamic_cast<PanWavelengthGrid*>(_lambdagrid);
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::setStellarSystem(StellarSystem* value)
{
    if (_ss) delete _ss;
    _ss = value;
    if (_ss) _ss->setParent(this);
}

////////////////////////////////////////////////////////////////////

StellarSystem* PanMonteCarloSimulation::stellarSystem() const
{
    return _ss;
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::setDustSystem(PanDustSystem* value)
{
    if (_ds) delete _ds;
    _ds = value;
    _pds = value;
    if (_ds) _ds->setParent(this);
}

////////////////////////////////////////////////////////////////////

PanDustSystem* PanMonteCarloSimulation::dustSystem() const
{
    return _pds;
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::runSelf()
{
    // If there are prepackages (used to do a coarse simulation to determine the dynamic grid)
    if(_prePackages > 0)
    {
        setPackages(_prePackages);
        runstellaremission();
        dynamicGrid();
        // properly resize the array used to communicate between rundustXXX() and the corresponding parallel loop
        _Ncells = _pds ? _pds->Ncells() : 0;
        if (_pds && _pds->dustemission()) _Labsbolv.resize(_Ncells);
        setPackages(_totalPackages); // Use the cached variable to go back to the normal amount of packages
    }
    runstellaremission();

    if (_pds && _pds->dustemission())
    {
        if (_pds && _pds->selfAbsorption()) rundustselfabsorption();
        rundustemission();
    }

    write();
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::rundustselfabsorption()
{
    TimeLogger logger(_log, "the dust self-absorption phase");

    // Initialize the total absorbed luminosity in the previous cycle
    double prevLabsdusttot = 0.;

    // Perform three "stages" of max 100 cycles each; the first stage uses 10 times less photon packages
    const int Nstages = 3;
    const char* stage_name[] = {"first-stage", "second-stage", "last-stage"};
    const double stage_factor[] = {1./10., 1./3., 1.};
    const double stage_epsmax[] = {0.010, 0.007, 0.005};
    for (int stage=0; stage<Nstages; stage++)
    {
        bool fixedNcycles = _pds->cycles();
        const int Ncyclesmax = fixedNcycles ? _pds->cycles() : 100;
        bool convergence = false;
        int cycle = 1;
        while (cycle<=Ncyclesmax && (!convergence || fixedNcycles))
        {
            TimeLogger logger(_log, "the " + QString(stage_name[stage]) + " dust self-absorption cycle "
                              + QString::number(cycle));

            // Construct the dust emission spectra
            _log->info("Calculating dust emission spectra...");
            _pds->calculatedustemission();
            _log->info("Dust emission spectra calculated.");

            // Determine the bolometric luminosity that is absorbed in every cell (and that will hence be re-emitted).
            _Labsbolv = _pds->Labsbolv();
            // Set the absorbed dust luminosity to zero in all cells
            _pds->rebootLabsdust();

            // Perform dust self-absorption, using the appropriate number of packages for the current stage
            setChunkParams(packages()*stage_factor[stage]);
            initprogress(QString(stage_name[stage]) + " dust self-absorption cycle " + QString::number(cycle));

            Parallel* parallel = find<ParallelFactory>()->parallel();
            if (_lambdagrid->assigner())
                parallel->call(this, &PanMonteCarloSimulation::dodustselfabsorptionchunk,
                               _lambdagrid->assigner(), _Nchunks);
            else
                parallel->call(this, &PanMonteCarloSimulation::dodustselfabsorptionchunk, _Nlambda, _Nchunks);

            // Wait for the other processes to reach this point
            _comm->wait("this self-absorption cycle");
            _pds->sumResults();

            // Determine and log the total absorbed luminosity in the vector Labstotv.
            double Labsdusttot = _pds->Labsdusttot();
            _log->info("The total absorbed stellar luminosity is "
                       + QString::number(_units->obolluminosity(_pds->Labsstellartot())) + " "
                       + _units->ubolluminosity() );
            _log->info("The total absorbed dust luminosity is "
                       + QString::number(_units->obolluminosity(Labsdusttot)) + " "
                       + _units->ubolluminosity() );

            // Check the criteria to terminate the self-absorption cycle:
            // - the total absorbed dust luminosity should change by less than epsmax compared to the previous cycle;
            // - the last stage must perform at least 2 cycles (to make sure that the energy is properly distributed)
            double eps = fabs((Labsdusttot-prevLabsdusttot)/Labsdusttot);
            prevLabsdusttot = Labsdusttot;
            if ( (stage<Nstages-1 || cycle>1) && eps<stage_epsmax[stage])
            {
                _log->info("Convergence reached; the last increase in the absorbed dust luminosity was "
                           + QString::number(eps*100, 'f', 2) + "%");
                convergence = true;
            }
            else
            {
                _log->info("Convergence not yet reached; the increase in the absorbed dust luminosity was "
                           + QString::number(eps*100, 'f', 2) + "%");
            }
            cycle++;
        }
        if (!convergence)
        {
            _log->error("Convergence not yet reached after " + QString::number(Ncyclesmax) + " "
                        + QString(stage_name[stage]) + " cycles!");
        }
    }
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::dodustselfabsorptionchunk(size_t index)
{
    // Determine the wavelength index for this chunk
    int ell = index % _Nlambda;

    // Determine the luminosity to be emitted at this wavelength index
    Array Lv(_Ncells);
    for (int m=0; m<_Ncells; m++)
    {
        double Labsbol = _Labsbolv[m];
        if (Labsbol>0.0) Lv[m] = Labsbol * _pds->dustluminosity(m,ell);
    }
    double Ltot = Lv.sum();

    // Emit photon packages
    if (Ltot > 0)
    {
        Array Xv;
        NR::cdf(Xv, Lv);

        PhotonPackage pp;
        double L = Ltot / _Npp;
        double Lthreshold = L / minWeightReduction();

        quint64 remaining = _chunksize;
        while (remaining > 0)
        {
            quint64 count = qMin(remaining, _logchunksize);
            for (quint64 i=0; i<count; i++)
            {
                double X = _random->uniform();
                int m = NR::locate_clip(Xv,X);
                Position bfr = _pds->randomPositionInCell(m);
                Direction bfk = _random->direction();
                pp.launch(L,ell,bfr,bfk);
                while (true)
                {
                    _pds->fillOpticalDepth(&pp);
                    simulateescapeandabsorption(&pp,true);
                    double L = pp.luminosity();
                    if (L==0.0) break;
                    if (L<=Lthreshold && pp.nScatt()>=minScattEvents()) break;
                    simulatepropagation(&pp);
                    simulatescattering(&pp);
                }
            }
            logprogress(count);
            remaining -= count;
        }
    }
    else logprogress(_chunksize);
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::rundustemission()
{
    TimeLogger logger(_log, "the dust emission phase");

    // Construct the dust emission spectra
    _log->info("Calculating dust emission spectra...");
    _pds->calculatedustemission();
    _log->info("Dust emission spectra calculated.");

    // Determine the bolometric luminosity that is absorbed in every cell (and that will hence be re-emitted).
    _Labsbolv = _pds->Labsbolv();

    // Perform the actual dust emission, possibly using more photon packages to obtain decent resolution
    setChunkParams(packages()*_pds->emissionBoost());
    initprogress("dust emission");
    Parallel* parallel = find<ParallelFactory>()->parallel();
    if (_lambdagrid->assigner())
        parallel->call(this, &PanMonteCarloSimulation::dodustemissionchunk, _lambdagrid->assigner(), _Nchunks);
    else
        parallel->call(this, &PanMonteCarloSimulation::dodustemissionchunk, _Nlambda, _Nchunks);

    // Wait for the other processes to reach this point
    _comm->wait("the dust emission phase");
}

////////////////////////////////////////////////////////////////////

void PanMonteCarloSimulation::dodustemissionchunk(size_t index)
{
    // Determine the wavelength index for this chunk
    int ell = index % _Nlambda;

    // Determine the luminosity to be emitted at this wavelength index
    Array Lv(_Ncells);
    for (int m=0; m<_Ncells; m++)
    {
        double Labsbol = _Labsbolv[m];
        if (Labsbol>0.0) Lv[m] = Labsbol * _pds->dustluminosity(m,ell);
    }
    double Ltot = Lv.sum();  // the total luminosity to be emitted at this wavelength index

    // Emit photon packages
    if (Ltot > 0)
    {
        // We consider biasing in the selection of the cell from which the photon packages are emitted.
        // A fraction of the cells is selected from the "natural" distribution, in which each cell is
        // weighted according to its total luminosity (Lv[em]). The other cells are selected from
        // a uniform distribution in which each cell has a equal probability.
        double xi = _pds->emissionBias();    // the fraction to be selected from a uniform distribution

        // the cumulative distribution of the natural pdf (Lv does not need to be normalized before)
        Array cumLv;
        NR::cdf(cumLv, Lv);

        PhotonPackage pp,ppp;
        double Lmean = Ltot/_Ncells;
        double Lem = Ltot / _Npp;
        double Lthreshold = Lem / minWeightReduction();

        quint64 remaining = _chunksize;
        while (remaining > 0)
        {
            quint64 count = qMin(remaining, _logchunksize);
            for (quint64 i=0; i<count; i++)
            {
                int m;
                double X = _random->uniform();
                if (X<xi)
                {
                    // rescale the deviate from [0,xi[ to [0,Ncells[
                    m = max(0,min(_Ncells-1,static_cast<int>(_Ncells*X/xi)));
                }
                else
                {
                    // rescale the deviate from [xi,1[ to [0,1[
                    m = NR::locate_clip(cumLv,(X-xi)/(1-xi));
                }
                double weight = 1.0/(1-xi+xi*Lmean/Lv[m]);
                Position bfr = _pds->randomPositionInCell(m);
                Direction bfk = _random->direction();
                pp.launch(Lem*weight,ell,bfr,bfk);
                peeloffemission(&pp,&ppp);
                while (true)
                {
                    _pds->fillOpticalDepth(&pp);
                    if (continuousScattering()) continuouspeeloffscattering(&pp,&ppp);
                    simulateescapeandabsorption(&pp,false);
                    double L = pp.luminosity();
                    if (L==0.0) break;
                    if (L<=Lthreshold && pp.nScatt()>=minScattEvents()) break;
                    simulatepropagation(&pp);
                    if (!continuousScattering()) peeloffscattering(&pp,&ppp);
                    simulatescattering(&pp);
                }
            }
            logprogress(count);
            remaining -= count;
        }
    }
    else logprogress(_chunksize);
}

////////////////////////////////////////////////////////////////////
