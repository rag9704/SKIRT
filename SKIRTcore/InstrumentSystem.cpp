/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FatalError.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////

InstrumentSystem::InstrumentSystem()
{
}

//////////////////////////////////////////////////////////////////////

void InstrumentSystem::addInstrument(Instrument* value)
{
    if (!value) throw FATALERROR("Instrument pointer shouldn't be null");
    value->setParent(this);
    _instruments << value;
}

//////////////////////////////////////////////////////////////////////

QList<Instrument*> InstrumentSystem::instruments() const
{
    return _instruments;
}

//////////////////////////////////////////////////////////////////////

void InstrumentSystem::write()
{
    foreach (Instrument* instrument, _instruments) instrument->write();
}

//////////////////////////////////////////////////////////////////////
