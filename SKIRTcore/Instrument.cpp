/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
//////////////////////////////////////////////////////////////////*/

#include "Instrument.hpp"
#include "DustSystem.hpp"
#include "FatalError.hpp"
#include "PhotonPackage.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////

Instrument::Instrument()
    : _ds(0)
{
}

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    try
    {
        _ds = find<DustSystem>();
    }
    catch (FatalError)
    {
        _ds = 0;
    }
}

////////////////////////////////////////////////////////////////////

void Instrument::setInstrumentName(QString value)
{
    _instrumentname = value;
}

////////////////////////////////////////////////////////////////////

QString Instrument::instrumentName() const
{
    return _instrumentname;
}

////////////////////////////////////////////////////////////////////

double Instrument::opticalDepth(const PhotonPackage* pp, double distance) const
{
    return _ds ? _ds->opticaldepth(pp->ell(),pp->position(),pp->direction(),distance) : 0;
}

////////////////////////////////////////////////////////////////////
