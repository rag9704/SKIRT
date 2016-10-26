/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTSYSTEM_HPP
#define INSTRUMENTSYSTEM_HPP

#include <vector>
#include <QPair>
#include "SimulationItem.hpp"
class Instrument;
class ParallelFactory;

//////////////////////////////////////////////////////////////////////

/** An InstrumentSystem instance keeps a list of zero or more instruments. The instruments can be
    of various nature (e.g. photometric, spectroscopic,...) and do not need to be located at the
    same observing position. */
class InstrumentSystem : public SimulationItem
{
    Q_OBJECT
    Q_CLASSINFO("Title", "an instrument system")

    Q_CLASSINFO("Property", "instruments")
    Q_CLASSINFO("Title", "the instruments")
    Q_CLASSINFO("Optional", "true")
    Q_CLASSINFO("Default", "SimpleInstrument")

    //============= Construction - Setup - Destruction =============

public:
    /** The default constructor; creates an empty instrument system. */
    Q_INVOKABLE InstrumentSystem();

    //======== Setters & Getters for Discoverable Attributes =======

public:
    /** This function inserts an instrument into the instrument system at the specified index. */
    Q_INVOKABLE void insertInstrument(int index, Instrument* value);

    /** This function removes the instrument with the specified index from the instrument system. */
    Q_INVOKABLE void removeInstrument(int index);

    /** This function returns the list of instruments in the instrument system. */
    Q_INVOKABLE QList<Instrument*> instruments() const;

    //======================== Other Functions =======================

public:
    /** This function writes down the results of the instrument system. It calls the write()
        function for each of the instruments. */
    void write();

    /** This function resets all detectors, putting all fluxes to zero */
    void reset();

    //======================== Data Members ========================

private:
    // discoverable attributes
    QList<Instrument*> _instruments;
};

////////////////////////////////////////////////////////////////////

#endif // INSTRUMENTSYSTEM_HPP
