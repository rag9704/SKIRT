/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Log.hpp"
#include "OligoDustSystem.hpp"
#include "OligoMonteCarloSimulation.hpp"
#include "OligoWavelengthGrid.hpp"
#include "StellarSystem.hpp"
#include "InstrumentSystem.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////

OligoMonteCarloSimulation::OligoMonteCarloSimulation()
{
}

////////////////////////////////////////////////////////////////////

void OligoMonteCarloSimulation::setWavelengthGrid(OligoWavelengthGrid* value)
{
    if (_lambdagrid) delete _lambdagrid;
    _lambdagrid = value;
    if (_lambdagrid) _lambdagrid->setParent(this);
}

////////////////////////////////////////////////////////////////////

OligoWavelengthGrid* OligoMonteCarloSimulation::wavelengthGrid() const
{
    return dynamic_cast<OligoWavelengthGrid*>(_lambdagrid);
}

////////////////////////////////////////////////////////////////////

void OligoMonteCarloSimulation::setStellarSystem(StellarSystem* value)
{
    if (_ss) delete _ss;
    _ss = value;
    if (_ss) _ss->setParent(this);
}

////////////////////////////////////////////////////////////////////

StellarSystem* OligoMonteCarloSimulation::stellarSystem() const
{
    return _ss;
}

////////////////////////////////////////////////////////////////////

void OligoMonteCarloSimulation::setDustSystem(OligoDustSystem* value)
{
    if (_ds) delete _ds;
    _ds = value;
    if (_ds) _ds->setParent(this);
}

////////////////////////////////////////////////////////////////////

OligoDustSystem* OligoMonteCarloSimulation::dustSystem() const
{
    return dynamic_cast<OligoDustSystem*>(_ds);
}

////////////////////////////////////////////////////////////////////

void OligoMonteCarloSimulation::runSelf()
{
    // If there are prepackages (used to do a coarse simulation to determine the dynamic grid)
    while(_dynamicIterations > 0)
    {
        _log->info("Starting prepackage iteration "+QString::number(_dynamicIterations));
        setPackages(_prePackages);
        runstellaremission();
        dynamicGrid();
        find<InstrumentSystem>()->reset(); // Reset the instruments
        _dynamicIterations--;
    }
    setPackages(_totalPackages); // Use the cached variable to go back to the normal amount of packages
    runstellaremission();

    write();
}

////////////////////////////////////////////////////////////////////
