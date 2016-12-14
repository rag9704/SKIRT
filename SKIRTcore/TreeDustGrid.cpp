/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include <cfloat>
#include <cmath>
#include "DustDistribution.hpp"
#include "DustGridPath.hpp"
#include "DustGridPlotFile.hpp"
#include "DustSystem.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "Random.hpp"
#include "TreeDustGrid.hpp"
#include "TreeNode.hpp"
#include "TreeNodeBoxDensityCalculator.hpp"
#include "TreeNodeSampleDensityCalculator.hpp"
#include "Units.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////

TreeDustGrid::TreeDustGrid()
    : _minlevel(0), _maxlevel(0),
      _search(TopDown), _Nrandom(100),
      _maxOpticalDepth(0), _maxMassFraction(0), _maxDensDispFraction(0),
      _assigner(0), _random(0), _parallel(0), _dd(0), _dmib(0),
      _totalmass(0), _eps(0),
      _Nnodes(0), _highestWriteLevel(0),
      _useDmibForSubdivide(false)
{
}

//////////////////////////////////////////////////////////////////////

TreeDustGrid::~TreeDustGrid()
{
    for (int l=0; l<_Nnodes; l++)
        delete _tree[l];
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setupSelfBefore()
{
    BoxDustGrid::setupSelfBefore();
    Log* log = find<Log>();

    // Validate attribute values

    if (_minlevel < 0) throw FATALERROR("The minimum tree level should be at least 0");
    if (_maxlevel < 2) throw FATALERROR("The maximum tree level should be at least 2");
    if (_maxlevel <= _minlevel) throw FATALERROR("Maximum tree level should be larger than minimum tree level");
    if (_Nrandom < 1) throw FATALERROR("Number of random samples must be at least 1");
    if (_maxOpticalDepth < 0.0) throw FATALERROR("The maximum mean optical depth should be positive");
    if (_maxMassFraction < 0.0) throw FATALERROR("The maximum mass fraction should be positive");
    if (_maxDensDispFraction < 0.0) throw FATALERROR("The maximum density dispersion fraction should be positive");
    if (_maxTempVolFraction < 0.0) throw FATALERROR("The maximum temperature * volume fraction should be positive");
    if (_maxTempGradVolFraction < 0.0) throw FATALERROR("The maximum temperature graident * volume fraction should "
                                                        "be positive");

    // Cache some often used values
    // A Parallel instance is created with a limited amount of threads (4) for performance reasons
    // or with a single thread in case of multiprocessing, to ensure consistency.
    size_t nthreads = find<PeerToPeerCommunicator>()->isMultiProc() ? 1 : 4;
    _random = find<Random>();
    _parallel = find<ParallelFactory>()->parallel(nthreads);
    _dd = find<DustDistribution>();
    _ds = find<DustSystem>();
    _dmib = _dd->interface<DustMassInBoxInterface>();
    _useDmibForSubdivide = _dmib && !_maxDensDispFraction;
    _totalmass = _dd->mass();
    _eps = 1e-12 * extent().widths().norm();
    _totalvolume = boundingbox().volume();

    // Create the root node

    _tree.push_back(createRoot(extent()));

    // Recursively subdivide the root node until all nodes satisfy the
    // necessary criteria. When finished, set the number _Nnodes.

    int currentlevel = -1;
    unsigned int l = 0;
    while (true)
    {
        TreeNode* node = _tree[l];
        int level = node->level();
        if (level>currentlevel)
        {
            log->info("Starting subdivision of level " + QString::number(level) + "...");
            currentlevel = level;
        }
        if (l%50000 == 0)
            log->info("Subdividing node number " + QString::number(l) + "...");
        if (node->ynchildless())
            subdivide(node);
        l++;
        if (l>=_tree.size()) break;
    }
    _Nnodes = _tree.size();

    // Construction of a vector _idv that contains the node IDs of all
    // leaves. This is the actual dust cell vector (only the leaves will
    // eventually become valid dust cells). We also create a vector
    // _cellnumberv with the cell numbers of all the nodes (i.e. the
    // rank m of the node in the vector _idv if the node is a leaf, and
    // -1 if the node is not a leaf).

    int m = 0;
    _cellnumberv.resize(_Nnodes,-1);
    for (int l=0; l<_Nnodes; l++)
    {
        if (_tree[l]->ynchildless())
        {
            _idv.push_back(l);
            _cellnumberv[l] = m;
            m++;
        }
    }
    int Ncells = _idv.size();

    // Log the number of cells

    log->info("Construction of the tree finished.");
    log->info("  Total number of nodes: " + QString::number(_Nnodes));
    log->info("  Total number of leaves: " + QString::number(Ncells));
    vector<int> countv(_maxlevel+1);
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = _tree[_idv[m]];
        int level = node->level();
        countv[level]++;
    }
    log->info("  Number of leaf cells of each level:");
    for (int level=0; level<=_maxlevel; level++)
        log->info("    Level " + QString::number(level) + ": " + QString::number(countv[level]) + " cells");

    // Determine the number of levels to be included in 3D grid output (if such output is requested)

    if (writeGrid())
    {
        int cumulativeCells = 0;
        for (_highestWriteLevel=0; _highestWriteLevel<=_maxlevel; _highestWriteLevel++)
        {
            cumulativeCells += countv[_highestWriteLevel];
            if (cumulativeCells > 1500) break;          // experimental number
        }
        if (_highestWriteLevel<_maxlevel)
            log->info("Will be outputting 3D grid data up to level " + QString::number(_highestWriteLevel) +
                      ", i.e. " + QString::number(cumulativeCells) + " cells.");
    }

    // Add neighbors to the tree structure (but only if required for the search method)

    if (_search == Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->addneighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortneighbors();
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::subdivide(TreeNode* node)
{
    // If level is below or at minlevel, there is always subdivision, and the subdivision is "regular"
    int level = node->level();
    if (level <= _minlevel)
    {
        node->createchildren(_tree.size());
        _tree.insert(_tree.end(), node->children().begin(), node->children().end());
    }

    // if level is below maxlevel, there may be subdivision depending on various stopping criteria
    else if (level < _maxlevel)
    {
        // construct an appropriate density calculator to estimate properties for stopping criteria and division
        TreeNodeDensityCalculator* calc;
        if (_useDmibForSubdivide)
        {
            // use the DustMassInBox interface
            calc = new TreeNodeBoxDensityCalculator(_dmib, node);
        }
        else
        {
            // sample the density in the cell
            TreeNodeSampleDensityCalculator* sampleCalc =
                    new TreeNodeSampleDensityCalculator(_random, _Nrandom, _dd, node);
            _parallel->call(sampleCalc, _Nrandom);
            calc = sampleCalc;
        }

        // if no stopping criteria are enabled, we keep subdividing indefinitely
        bool needDivision = (_maxOpticalDepth == 0 && _maxMassFraction == 0 && _maxDensDispFraction == 0);

        // otherwise, we subdivide if at least one stopping criterion is not satisfied

        // check mass fraction
        if (!needDivision && _maxMassFraction > 0)
        {
            double massfraction = calc->mass() / _totalmass;
            if (massfraction >= _maxMassFraction) needDivision = true;
        }

        // check optical depth
        if (!needDivision && _maxOpticalDepth > 0)
        {
            double opticaldepth = calc->opticalDepth();
            if (opticaldepth >= _maxOpticalDepth) needDivision = true;
        }

        // check density dispersion fraction
        if (!needDivision && _maxDensDispFraction > 0)
        {
            double densdispfraction = calc->densityDispersion();
            if (densdispfraction >= _maxDensDispFraction) needDivision = true;
        }

        if (needDivision)
        {
            // there is subdivision, possibly using calculated properties such as barycenter
            node->createchildren(_tree.size(), calc);
            _tree.insert(_tree.end(), node->children().begin(), node->children().end());
        }

        delete calc;
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMinLevel(int value)
{
    _minlevel = value;
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::minLevel() const
{
    return _minlevel;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxLevel(int value)
{
    _maxlevel = value;
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::maxLevel() const
{
    return _maxlevel;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setSearchMethod(TreeDustGrid::SearchMethod value)
{
    _search = value;
}

//////////////////////////////////////////////////////////////////////

TreeDustGrid::SearchMethod TreeDustGrid::searchMethod() const
{
    return _search;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setSampleCount(int value)
{
    _Nrandom = value;
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::sampleCount() const
{
    return _Nrandom;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxOpticalDepth(double value)
{
    _maxOpticalDepth = value;
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::maxOpticalDepth() const
{
    return _maxOpticalDepth;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxMassFraction(double value)
{
    _maxMassFraction = value;
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::maxMassFraction() const
{
    return _maxMassFraction;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxDensDispFraction(double value)
{
    _maxDensDispFraction = value;
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::maxDensDispFraction() const
{
    return _maxDensDispFraction;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxTempVolFraction(double value)
{
    _maxTempVolFraction = value;
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::maxTempVolFraction() const
{
    return _maxTempVolFraction;
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::setMaxTempGradVolFraction(double value)
{
    _maxTempGradVolFraction = value;
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::maxTempGradVolFraction() const
{
    return _maxTempGradVolFraction;
}

////////////////////////////////////////////////////////////////////

double TreeDustGrid::volume(int m) const
{
    if (m<0 || m>numCells())
        throw FATALERROR("Invalid cell number: " + QString::number(m));
    TreeNode* node = getnode(m);
    return node->xwidth() * node->ywidth() * node->zwidth();
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::numCells() const
{
    return _idv.size();
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::whichcell(Position bfr) const
{
    const TreeNode* node = root()->whichnode(bfr);
    return node ? cellnumber(node) : -1;
}

//////////////////////////////////////////////////////////////////////

Position TreeDustGrid::centralPositionInCell(int m) const
{
    return Position(getnode(m)->extent().center());
}

//////////////////////////////////////////////////////////////////////

Position TreeDustGrid::randomPositionInCell(int m) const
{
    return _random->position(getnode(m)->extent());
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::path(DustGridPath* path) const
{
    // Initialize the path
    path->clear();

    // If the photon package starts outside the dust grid, move it into the first grid cell that it will pass
    Position bfr = path->moveInside(extent(), _eps);

    // Get the node containing the current location;
    // if the position is not inside the grid, return an empty path
    const TreeNode* node = root()->whichnode(bfr);
    if (!node) return path->clear();

    // Start the loop over nodes/path segments until we leave the grid.
    // Use a different code segment depending on the search method.
    double x,y,z;
    bfr.cartesian(x,y,z);
    double kx,ky,kz;
    path->direction().cartesian(kx,ky,kz);

    // ----------- Top-down -----------

    if (_search == TopDown)
    {
        while (node)
        {
            double xnext = (kx<0.0) ? node->xmin() : node->xmax();
            double ynext = (ky<0.0) ? node->ymin() : node->ymax();
            double znext = (kz<0.0) ? node->zmin() : node->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            double ds;
            if (dsx<=dsy && dsx<=dsz) ds = dsx;
            else if (dsy<=dsx && dsy<=dsz) ds = dsy;
            else ds = dsz;
            path->addSegment(cellnumber(node), ds);
            x += (ds+_eps)*kx;
            y += (ds+_eps)*ky;
            z += (ds+_eps)*kz;

            // always search from the root node down
            const TreeNode* oldnode = node;
            node = root()->whichnode(Vec(x,y,z));

            // if we're stuck in the same node...
            if (node==oldnode)
            {
                // try to escape by advancing the position to the next representable coordinates
                find<Log>()->warning("Photon package seems stuck in dust cell "
                                     + QString::number(node->id()) + " -- escaping");
                x = nextafter(x, (kx<0.0) ? -DBL_MAX : DBL_MAX);
                y = nextafter(y, (ky<0.0) ? -DBL_MAX : DBL_MAX);
                z = nextafter(z, (kz<0.0) ? -DBL_MAX : DBL_MAX);
                node = root()->whichnode(Vec(x,y,z));

                // if that didn't work, terminate the path
                if (node==oldnode)
                {
                    find<Log>()->warning("Photon package is stuck in dust cell "
                                         + QString::number(node->id()) + " -- terminating this path");
                    break;
                }
            }
        }
    }

    // ----------- Neighbor -----------

    else if (_search == Neighbor)
    {
        while (node)
        {
            double xnext = (kx<0.0) ? node->xmin() : node->xmax();
            double ynext = (ky<0.0) ? node->ymin() : node->ymax();
            double znext = (kz<0.0) ? node->zmin() : node->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            double ds;
            TreeNode::Wall wall;
            if (dsx<=dsy && dsx<=dsz)
            {
                ds = dsx;
                wall = (kx<0.0) ? TreeNode::BACK : TreeNode::FRONT;
            }
            else if (dsy<=dsx && dsy<=dsz)
            {
                ds = dsy;
                wall = (ky<0.0) ? TreeNode::LEFT : TreeNode::RIGHT;
            }
            else
            {
                ds = dsz;
                wall = (kz<0.0) ? TreeNode::BOTTOM : TreeNode::TOP;
            }
            path->addSegment(cellnumber(node), ds);
            x += (ds+_eps)*kx;
            y += (ds+_eps)*ky;
            z += (ds+_eps)*kz;

            // attempt to find the new node among the neighbors of the current node;
            // this should not fail unless the new location is outside the grid,
            // however on rare occasions it fails due to rounding errors (e.g. in a corner),
            // thus we use top-down search as a fall-back
            const TreeNode* oldnode = node;
            node = node->whichnode(wall, Vec(x,y,z));
            if (!node) node = root()->whichnode(Vec(x,y,z));

            // if we're stuck in the same node...
            if (node==oldnode)
            {
                // try to escape by advancing the position to the next representable coordinates
                find<Log>()->warning("Photon package seems stuck in dust cell "
                                     + QString::number(node->id()) + " -- escaping");
                x = nextafter(x, (kx<0.0) ? -DBL_MAX : DBL_MAX);
                y = nextafter(y, (ky<0.0) ? -DBL_MAX : DBL_MAX);
                z = nextafter(z, (kz<0.0) ? -DBL_MAX : DBL_MAX);
                node = root()->whichnode(Vec(x,y,z));

                // if that didn't work, terminate the path
                if (node==oldnode)
                {
                    find<Log>()->warning("Photon package is stuck in dust cell "
                                         + QString::number(node->id()) + " -- terminating this path");
                    break;
                }
            }
        }
    }

    // ----------- Bookkeeping -----------

    // !! This code section relies on the fact that an octtree node is used !!

    else if (_search == Bookkeeping)
    {
        int l = node->id();  // node index in the tree
        while (true)
        {
            double xnext = (kx<0.0) ? _tree[l]->xmin() : _tree[l]->xmax();
            double ynext = (ky<0.0) ? _tree[l]->ymin() : _tree[l]->ymax();
            double znext = (kz<0.0) ? _tree[l]->zmin() : _tree[l]->zmax();
            double dsx = (fabs(kx)>1e-15) ? (xnext-x)/kx : DBL_MAX;
            double dsy = (fabs(ky)>1e-15) ? (ynext-y)/ky : DBL_MAX;
            double dsz = (fabs(kz)>1e-15) ? (znext-z)/kz : DBL_MAX;

            // First option: the x-wall is hit first. After moving
            // towards the boundary, we have to find the next cell. First we
            // check whether the node is on the right or left side of his
            // father node. If the movement is towards positive x (i.e. if
            // kx>0) we move up in the tree until we find a node on the left
            // side. The next cell will then be the corresponding right node
            // (if it is a leaf) or one of its children. If we have to move
            // up until we hit the root node, this means our path has ended.

            if (dsx<=dsy && dsx<=dsz)
            {
                path->addSegment(_cellnumberv[l], dsx);
                x = xnext;
                y += ky*dsx;
                z += kz*dsx;
                while (true)
                {
                    int oct = ((l - 1) % 8) + 1;
                    bool place = (kx<0.0) ? (oct % 2 == 1) : (oct % 2 == 0);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (kx<0.0) ? -1 : 1;
                while (_cellnumberv[l] == -1)
                {
                    double yM = _tree[l]->child(0)->ymax();
                    double zM = _tree[l]->child(0)->zmax();
                    if (kx<0.0)
                    {
                        if (y<=yM)
                            l = (z<=zM) ? _tree[l]->child(1)->id() : _tree[l]->child(5)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(3)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (y<=yM)
                            l = (z<=zM) ? _tree[l]->child(0)->id() : _tree[l]->child(4)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(2)->id() : _tree[l]->child(6)->id();
                    }
                }
            }

            // Repeat the same exercise, but now the y-wall is hit first...

            else if (dsy<dsx && dsy<=dsz)
            {
                path->addSegment(_cellnumberv[l], dsy);
                x += kx*dsy;
                y  = ynext;
                z += kz*dsy;
                while (true)
                {
                    bool place = (ky<0.0) ? ((l-1) % 4 < 2) : ((l-1) % 4 > 1);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (ky<0.0) ? -2 : 2;
                while (_cellnumberv[l] == -1)
                {
                    double xM = _tree[l]->child(0)->xmax();
                    double zM = _tree[l]->child(0)->zmax();
                    if (ky<0.0)
                    {
                        if (x<=xM)
                            l = (z<=zM) ? _tree[l]->child(2)->id() : _tree[l]->child(6)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(3)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (x<=xM)
                            l = (z<=zM) ? _tree[l]->child(0)->id() : _tree[l]->child(4)->id();
                        else
                            l = (z<=zM) ? _tree[l]->child(1)->id() : _tree[l]->child(5)->id();
                    }
                }
            }

            // Finally, repeat the same exercise, but now the z-wall is hit first...

            else if (dsz< dsx && dsz< dsy)
            {
                path->addSegment(_cellnumberv[l], dsz);
                x += kx*dsz;
                y += ky*dsz;
                z  = znext;
                while (true)
                {
                    int oct = ((l-1) % 8) + 1;
                    bool place = (kz<0.0) ? (oct < 5) : (oct > 4);
                    if (!place) break;
                    l = _tree[l]->father()->id();
                    if (l == 0) return;
                }
                l += (kz<0.0) ? -4 : 4;
                while (_cellnumberv[l] == -1)
                {
                    double xM = _tree[l]->child(0)->xmax();
                    double yM = _tree[l]->child(0)->ymax();
                    if (kz<0.0)
                    {
                        if (x<=xM)
                            l = (y<=yM) ? _tree[l]->child(4)->id() : _tree[l]->child(6)->id();
                        else
                            l = (y<=yM) ? _tree[l]->child(5)->id() : _tree[l]->child(7)->id();
                    }
                    else
                    {
                        if (x<=xM)
                            l = (y<=yM) ? _tree[l]->child(0)->id() : _tree[l]->child(2)->id();
                        else
                            l = (y<=yM) ? _tree[l]->child(1)->id() : _tree[l]->child(3)->id();
                    }
                }
            }
        }
    }

    // ------------------------------
}

//////////////////////////////////////////////////////////////////////

QList<SimulationItem*> TreeDustGrid::interfaceCandidates(const type_info &interfaceTypeInfo)
{
    if (interfaceTypeInfo == typeid(DustGridDensityInterface) && !_dmib)
        return QList<SimulationItem*>();
    return BoxDustGrid::interfaceCandidates(interfaceTypeInfo);
}

//////////////////////////////////////////////////////////////////////

double TreeDustGrid::density(int h, int m) const
{
    TreeNode* node = getnode(m);
    return _dmib->massInBox(h, node->extent()) / node->volume();
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::subdivideTemperatureRecursive(TreeNode* node, double vol, double temp, double tempGrad)
{
    // If level is below or at minlevel, there is always subdivision, and the subdivision is "regular"
    int level = node->level();
    int cellNr = cellnumber(node);

    if (vol == 0 && tempGrad == 0)
    {
        vol = _ds->volume(cellNr);
        temp = _ds->temperature(cellNr);
    }

    // if level is below maxlevel, there may be subdivision depending on various stopping criteria
    if (level < _maxlevel && node->ynchildless())
    {
        // default to no subdivision, unless there is a stopping criteria
        bool needDivision = false;

        // check temperature * volume fraction
        if (!needDivision && _maxTempVolFraction > 0)
        {
            double tempVolfraction = temp*vol/_totalvolume;
            if (tempVolfraction >= _maxTempVolFraction) needDivision = true;
        }

        // check tempGrad * volume fraction
        if (!needDivision && _maxTempGradVolFraction > 0)
        {
            // The temperature gradient starts by calculating the average temperature for
            // each wall (at the neighbor's side). This takes into account the volume of
            // the neighboring cells (a neighbor which subdivision level is one more than the
            // current cell gets a weight of 1/4, since there are 4 of those neighbors).
            // The temperature gradient is then defined as the average (over all 6 walls)
            // of the absolute difference between the average wall temperature and the
            // current cell temperature. So:
            // tempGrad = (SUM_wall abs(T_wall-T_cell)) / 6
            if (tempGrad == 0)
            {
                // Loop over all 6 walls
                for (int wallInt = 0; wallInt < 6; wallInt++)
                {
                    TreeNode::Wall wall = static_cast<TreeNode::Wall>(wallInt);
                    // The mean wall temperature is acquired by multiplying the temperature of
                    // each neighbor with it's weight (1/4^{wall_level-cell_level}), and summing
                    // over all these weighted neighbor temperatures.
                    double meanWallT = 0;
                    vector<TreeNode*> neighbors = node->getneighbors(wall);
                    for (unsigned int i = 0 ; i < neighbors.size() ; i++)
                    {
                        int neighborcellNr = cellnumber(neighbors[i]);
                        meanWallT += _ds->temperature(neighborcellNr)/pow(4, neighbors[i]->level() - node->level());
                    }
                    tempGrad += abs(meanWallT - temp);
                }
                tempGrad /= 6; // 6 walls: divide by 6 to average
            }
            double tempGradVolFraction = tempGrad * vol / _totalvolume;
            if (tempGradVolFraction >= _maxTempGradVolFraction) needDivision = true;
        }

        if (needDivision)
        {
            // there is subdivision. Note that this does use regular subdivision
            // So this needs further attention when using barycentric trees...
            node->createchildren(_tree.size());
            _tree.insert(_tree.end(), node->children().begin(), node->children().end());
            // Recursion
            vector<TreeNode*> children = node->children();
            int childSize = children.size(); // 8 for octree, 2 for bintree
            for (int i = 0 ; i < childSize; i++)
            {
                // Volume is extensive, temperature is intensive.
                // Temperature gradient depends on how the temperature varies, but is
                // assumed to be intensive here (by definition, since the tempGrad is just a heuristic).
                subdivideTemperatureRecursive(children[i], vol/childSize, temp, tempGrad);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::dynamicGrid()
{
    Log* log = find<Log>();
    log->info("Dynamic grid: subdividing the tree dust grid based on radiation criteria...");
    _ds->calculateTemperature(); // Calculates the temperature in each dust cell (saved in _ds)
    // Loop over all leaf nodes
    for (unsigned int i = 0 ; i < _idv.size() ; i++)
    {
        TreeNode* node = _tree[_idv[i]];
        // Do recursive subdivision of the node
        // Use of recursion is done since a node already calculates its volume, temperature and
        // temperature gradient, and these are passed to each child (if subdivided).
        // These can not be calculated from the child itself, since the child is not (yet) part
        // of the dust grid.
        subdivideTemperatureRecursive(node);
    }
    // -- Pretty much copy paste from setupSelfBefore()... --
    _Nnodes = _tree.size();
    // Reconstruct _cellnumberv and _idv
    int m = 0;
    _cellnumberv.resize(_Nnodes,-1);
    _idv.clear();
    for (int l=0; l<_Nnodes; l++)
    {
        if (_tree[l]->ynchildless())
        {
            _idv.push_back(l);
            _cellnumberv[l] = m;
            m++;
        }
    }
    int Ncells = _idv.size();

    // Log the number of cells

    log->info("Reconstruction of the tree finished.");
    log->info("  Total number of nodes: " + QString::number(_Nnodes));
    log->info("  Total number of leaves: " + QString::number(Ncells));
    vector<int> countv(_maxlevel+1);
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = _tree[_idv[m]];
        int level = node->level();
        countv[level]++;
    }
    log->info("  Number of leaf cells of each level:");
    for (int level=0; level<=_maxlevel; level++)
        log->info("    Level " + QString::number(level) + ": " + QString::number(countv[level]) + " cells");

    // Determine the number of levels to be included in 3D grid output (if such output is requested)

    if (writeGrid())
    {
        int cumulativeCells = 0;
        for (_highestWriteLevel=0; _highestWriteLevel<=_maxlevel; _highestWriteLevel++)
        {
            cumulativeCells += countv[_highestWriteLevel];
            if (cumulativeCells > 1500) break;          // experimental number
        }
        if (_highestWriteLevel<_maxlevel)
            log->info("Will be outputting 3D grid data up to level " + QString::number(_highestWriteLevel) +
                      ", i.e. " + QString::number(cumulativeCells) + " cells.");
    }

    // Add neighbors to the tree structure (but only if required for the search method)

    if (_search == Neighbor)
    {
        log->info("Adding neighbors to the tree nodes...");
        for (int l=0; l<_Nnodes; l++) _tree[l]->deleteallneighbors(); // First start by resetting all neighbor lists
        for (int l=0; l<_Nnodes; l++) _tree[l]->addneighbors();
        for (int l=0; l<_Nnodes; l++) _tree[l]->sortneighbors();
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xy(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_xmin, _ymin, _xmax, _ymax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getnode(m);
        if (fabs(node->zmin()) < 1e-8*extent().zwidth())
        {
            outfile->writeRectangle(node->xmin(), node->ymin(), node->xmax(), node->ymax());
        }
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xz(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_xmin, _zmin, _xmax, _zmax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getnode(m);
        if (fabs(node->ymin()) < 1e-8*extent().ywidth())
        {
            outfile->writeRectangle(node->xmin(), node->zmin(), node->xmax(), node->zmax());
        }
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_yz(DustGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(_ymin, _zmin, _ymax, _zmax);
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getnode(m);
        if (fabs(node->xmin()) < 1e-8*extent().xwidth())
        {
            outfile->writeRectangle(node->ymin(), node->zmin(), node->ymax(), node->zmax());
        }
    }
}

//////////////////////////////////////////////////////////////////////

void TreeDustGrid::write_xyz(DustGridPlotFile* outfile) const
{
    // Output all leaf cells up to a certain level
    int Ncells = numCells();
    for (int m=0; m<Ncells; m++)
    {
        TreeNode* node = getnode(m);
        if (node->level() <= _highestWriteLevel)
            outfile->writeCube(node->xmin(), node->ymin(), node->zmin(), node->xmax(), node->ymax(), node->zmax());
    }
}

//////////////////////////////////////////////////////////////////////

TreeNode* TreeDustGrid::root() const
{
    return _tree[0];
}

//////////////////////////////////////////////////////////////////////

TreeNode* TreeDustGrid::getnode(int m) const
{
    return _tree[_idv[m]];
}

//////////////////////////////////////////////////////////////////////

int TreeDustGrid::cellnumber(const TreeNode* node) const
{
    return _cellnumberv[node->id()];
}

//////////////////////////////////////////////////////////////////////
