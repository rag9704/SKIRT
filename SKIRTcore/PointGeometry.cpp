/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include <limits>
#include "PointGeometry.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////

PointGeometry::PointGeometry()
{
}

//////////////////////////////////////////////////////////////////////

int
PointGeometry::dimension()
const
{
    return 0;
}

//////////////////////////////////////////////////////////////////////

double
PointGeometry::density(Position bfr)
const
{
    return bfr.radius() == 0 ? numeric_limits<double>::infinity() : 0;
}

//////////////////////////////////////////////////////////////////////

Position
PointGeometry::generatePosition()
const
{
    return Position();
}

//////////////////////////////////////////////////////////////////////

double
PointGeometry::SigmaX()
const
{
    return numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double
PointGeometry::SigmaY()
const
{
    return numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double
PointGeometry::SigmaZ()
const
{
    return numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////
