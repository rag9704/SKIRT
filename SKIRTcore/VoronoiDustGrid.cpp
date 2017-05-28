/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include <cmath>
#include "DustDistribution.hpp"
#include "DustGridPlotFile.hpp"
#include "DustParticleInterface.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "Units.hpp"
#include "VoronoiDustGrid.hpp"
#include "VoronoiMesh.hpp"
#include "VoronoiMeshFile.hpp"
#include "VoronoiMeshInterface.hpp"
#include "PanMonteCarloSimulation.hpp"
#include "container.hh"
#include "PanDustSystem.hpp"
#include "TextOutFile.hpp"
#include "DustMix.hpp"
#include "NR.hpp"
#include "WavelengthGrid.hpp"
#include "OligoDustSystem.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////

VoronoiDustGrid::VoronoiDustGrid()
    : _numParticles(0), _distribution(DustDensity), _meshfile(0), _random(0), _mesh(0), _meshOwned(true)
{
}

//////////////////////////////////////////////////////////////////////

VoronoiDustGrid::~VoronoiDustGrid()
{
    if (_meshOwned) delete _mesh;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setupSelfBefore()
{
    BoxDustGrid::setupSelfBefore();

    // Cache the random number generator
    _random = find<Random>();

    Log* log = find<Log>();

    // Prepackage phase (for setting up the dynamic grid)
    _totalNumParticles = _numParticles;
    if(_tempDistFraction+_tempGradFraction+_tempMassFraction > 0)
    {
        // Use a fraction of the total particles for the prepackage phase
        _numParticles = ceil(_preGridPointFraction*_totalNumParticles);
        log->info("Prepackage phase: Using " + QString::number(_numParticles) + " of the total "+
                  QString::number(_totalNumParticles)+" grid points to calculate a temperature distribution.");
    }

    // Determine an appropriate set of particles and construct the Voronoi mesh
    switch (_distribution)
    {
    case Uniform:
        {
            if (_numParticles < 10) throw FATALERROR("The number of particles should be at least 10");
            vector<Vec> rv(_numParticles);
            for (int m=0; m<_numParticles; m++)
            {
                rv[m] = _random->position(extent());
            }
            log->info("Computing Voronoi tesselation for " + QString::number(_numParticles)
                      + " uniformly distributed random particles...");
            _mesh = new VoronoiMesh(rv, extent(), log, _relaxationSteps);
            break;
        }
    case CentralPeak:
        {
            if (_numParticles < 10) throw FATALERROR("The number of particles should be at least 10");
            const int a = 1000;                     // steepness of the peak; the central 1/a portion is NOT covered
            const double rscale = extent().rmax().norm();
            vector<Vec> rv(_numParticles);
            for (int m=1; m<_numParticles; m++)     // skip first particle so that it remains (0,0,0)
            {
                while (true)
                {
                    double r = rscale * pow(1./a, _random->uniform());   // random distribution according to 1/x
                    Direction k = _random->direction();
                    Position p = Position(r,k);
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            log->info("Computing Voronoi tesselation for " + QString::number(_numParticles)
                      + " random particles distributed in a central peak...");
            _mesh = new VoronoiMesh(rv, extent(), log, _relaxationSteps);
            break;
        }
    case DustDensity:
        {
            if (_numParticles < 10) throw FATALERROR("The number of particles should be at least 10");
            DustDistribution* dd = find<DustDistribution>();
            vector<Vec> rv(_numParticles);
            for (int m=0; m<_numParticles; m++)
            {
                while (true)
                {
                    Position p = dd->generatePosition();
                    if (extent().contains(p))       // discard any points outside of the domain
                    {
                        rv[m] = p;
                        break;
                    }
                }
            }
            log->info("Computing Voronoi tesselation for " + QString::number(_numParticles)
                      + " random particles distributed according to dust density...");
            _mesh = new VoronoiMesh(rv, extent(), log, _relaxationSteps);
            break;
        }
    case DustTesselation:
        {
            VoronoiMeshInterface* vmi = find<DustDistribution>()->interface<VoronoiMeshInterface>();
            if (!vmi) throw FATALERROR("Can't retrieve Voronoi mesh from this dust distribution");
            _mesh = vmi->mesh();
            _meshOwned = false;
            log->info("Using Voronoi tesselation from dust distribution with " + QString::number(_mesh->Ncells())
                      + " particles...");
            break;
        }
    case SPHParticles:
        {
            DustParticleInterface* dpi = find<DustDistribution>()->interface<DustParticleInterface>();
            if (!dpi) throw FATALERROR("Can't retrieve particle locations from this dust distribution");
            log->info("Computing Voronoi tesselation for " + QString::number(dpi->numParticles())
                      + " dust distribution particles...");
            _mesh = new VoronoiMesh(dpi, extent(), log, _relaxationSteps);
            break;
        }
    case File:
        {
            if (!_meshfile) throw FATALERROR("File containing particle locations is not defined");
            log->info("Computing Voronoi tesselation for particles loaded from file " + _meshfile->filename() + "...");
            _mesh = new VoronoiMesh(_meshfile, QList<int>(), extent(), log, _relaxationSteps);
            break;
        }
    default:
        throw FATALERROR("Unknown distribution type");
    }

    int Ncells = _mesh->Ncells();

    // Log statistics on the cell neighbors
    double avgNeighbors;
    int minNeighbors, maxNeighbors;
    _mesh->neighborStatistics(avgNeighbors, minNeighbors, maxNeighbors);
    log->info("Computed Voronoi tesselation with " + QString::number(Ncells) + " cells:");
    log->info("  Average number of neighbors per cell: " + QString::number(avgNeighbors,'f',1));
    log->info("  Minimum number of neighbors per cell: " + QString::number(minNeighbors));
    log->info("  Maximum number of neighbors per cell: " + QString::number(maxNeighbors));

    // Log statistics on the block lists
    int nblocks = _mesh->Nblocks();
    double avgRefsPerBlock;
    int minRefsPerBlock, maxRefsPerBlock;
    _mesh->blockStatistics(avgRefsPerBlock, minRefsPerBlock, maxRefsPerBlock);
    log->info("Created grid to accelerate which-cell operations:");
    log->info("  Number of cells                  : " + QString::number(Ncells));
    log->info("  Number of blocks                 : " + QString::number(nblocks*nblocks*nblocks) +
              " (" + QString::number(nblocks) + " in each dimension)");
    log->info("  Average number of cells per block: " + QString::number(avgRefsPerBlock,'f',1));
    log->info("  Minimum number of cells per block: " + QString::number(minRefsPerBlock));
    log->info("  Maximum number of cells per block: " + QString::number(maxRefsPerBlock));

    // Log statistics on the search trees
    double avgRefsPerTree;
    int nTrees, minRefsPerTree, maxRefsPerTree;
    _mesh->treeStatistics(nTrees, avgRefsPerTree, minRefsPerTree, maxRefsPerTree);
    log->info("Created search trees to accelerate which-cell operations:");
    log->info("  Number of trees                  : " + QString::number(nTrees) +
              " (" + QString::number(100.*nTrees/(nblocks*nblocks*nblocks),'f',1) + "% of blocks)");
    log->info("  Average number of cells per tree : " + QString::number(avgRefsPerTree,'f',1));
    log->info("  Minimum number of cells per tree : " + QString::number(minRefsPerTree));
    log->info("  Maximum number of cells per tree : " + QString::number(maxRefsPerTree));

    // If requested, output the plot files (we have to reconstruct the Voronoi tesselation...)
    if (writeGrid())
    {
        setWriteGrid(false);    // keep the base class from overwriting our plot files

        // create the plot files
        DustGridPlotFile plotxy(this, "ds_gridxy");
        DustGridPlotFile plotxz(this, "ds_gridxz");
        DustGridPlotFile plotyz(this, "ds_gridyz");
        DustGridPlotFile plotxyz(this, "ds_gridxyz");

        // load all particles in a Voro container
        int nb = max(3, min(1000, static_cast<int>(pow(Ncells/5.,1./3.)) ));
        voro::container con(xmin(), xmax(), ymin(), ymax(), zmin(), zmax(), nb, nb, nb, false,false,false, 8);
        for (int m=0; m<Ncells; m++)
        {
            Vec r = _mesh->particlePosition(m);
            con.put(m, r.x(),r.y(),r.z());
        }

        // loop over all Voro cells
        voro::c_loop_all loop(con);
        if (loop.start()) do
        {
            // Compute the cell
            voro::voronoicell fullcell;
            con.compute_cell(fullcell, loop);

            // Get the edges of the cell
            double x,y,z;
            loop.pos(x,y,z);
            vector<double> coords;
            fullcell.vertices(x,y,z, coords);
            vector<int> indices;
            fullcell.face_vertices(indices);

            // Write the edges of the cell to the plot files
            Box bounds = _mesh->extent(loop.pid());
            if (bounds.zmin()<=0 && bounds.zmax()>=0) plotxy.writePolyhedron(coords, indices);
            if (bounds.ymin()<=0 && bounds.ymax()>=0) plotxz.writePolyhedron(coords, indices);
            if (bounds.xmin()<=0 && bounds.xmax()>=0) plotyz.writePolyhedron(coords, indices);
            if (loop.pid() <= 1000) plotxyz.writePolyhedron(coords, indices);
        }
        while (loop.inc());
    }
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setupSelfAfter()
{
    if(!find<WavelengthGrid>()->issampledrange()) // OLIGO
    {
        if(_tempDistFraction > 0 && !find<OligoDustSystem>()->writeMeanIntensity())
            throw FATALERROR("Set writeMeanIntensity to true when drawing from a mean intensity distribution");
    }
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setNumParticles(int value)
{
    _numParticles = value;
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::numParticles() const
{
    return _numParticles;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setDistribution(VoronoiDustGrid::Distribution value)
{
    _distribution = value;
}

//////////////////////////////////////////////////////////////////////

VoronoiDustGrid::Distribution VoronoiDustGrid::distribution() const
{
    return _distribution;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setRelaxationSteps(int value)
{
    _relaxationSteps = value;
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::relaxationSteps() const
{
    return _relaxationSteps;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setPreGridPointFraction(double value)
{
    _preGridPointFraction = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::preGridPointFraction() const
{
    return _preGridPointFraction;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setTempDistFraction(double value)
{
    _tempDistFraction = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::tempDistFraction() const
{
    return _tempDistFraction;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setTempGradFraction(double value)
{
    _tempGradFraction = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::tempGradFraction() const
{
    return _tempGradFraction;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setTempMassFraction(double value)
{
    _tempMassFraction = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::tempMassFraction() const
{
    return _tempMassFraction;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setTempImportance(double value)
{
    _tempImportance = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::tempImportance() const
{
    return _tempImportance;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setMassImportance(double value)
{
    _massImportance = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::massImportance() const
{
    return _massImportance;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setTempGradImportance(double value)
{
    _tempGradImportance = value;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::tempGradImportance() const
{
    return _tempGradImportance;
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::setVoronoiMeshFile(VoronoiMeshFile* value)
{
    if (_meshfile) delete _meshfile;
    _meshfile = value;
    if (_meshfile) _meshfile->setParent(this);
}

//////////////////////////////////////////////////////////////////////

VoronoiMeshFile* VoronoiDustGrid::voronoiMeshFile() const
{
    return _meshfile;
}

//////////////////////////////////////////////////////////////////////

double VoronoiDustGrid::volume(int m) const
{
    return _mesh->volume(m);
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::numCells() const
{
    return _mesh->Ncells();
}

//////////////////////////////////////////////////////////////////////

int VoronoiDustGrid::whichcell(Position bfr) const
{
    return _mesh->cellIndex(bfr);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiDustGrid::centralPositionInCell(int m) const
{
    return _mesh->centralPosition(m);
}

//////////////////////////////////////////////////////////////////////

Position VoronoiDustGrid::randomPositionInCell(int m) const
{
    return _mesh->randomPosition(_random, m);
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::path(DustGridPath* path) const
{
    _mesh->path(path);
}

//////////////////////////////////////////////////////////////////////

void VoronoiDustGrid::drawFromTemperatureDistribution()
{
    // If we don't draw from a dynamic distribution, return
    // Commented out so that if we use prepackages, we still draw from a dust distribution (and not Uniform)
    //if(_tempDistFraction+_tempGradFraction+_tempMassFraction == 0)
    //    return;
    Log* log = find<Log>();
    log->info("Drawing extra voronoi points from a dynamic distribution.");
    DustSystem* ds = find<DustSystem>(); // Cache the dust system

    // Construct cdf of V*T (volume since we want to sample more from bigger cells, if temperature is equal)
    // Start with a V*T vector
    Array VTv; // Volume * Temperature for every cell (_numParticles is the old grid size)
    VTv.resize(_numParticles);
    ds->calculateTemperature();
    for(int m=0; m<_numParticles; m++)
    {
        VTv[m] = ds->temperature(m)*ds->volume(m);
    }

    // Temperature gradient
    Array VgradTv; // V*gradT vector
    VgradTv.resize(_numParticles);
    if (_tempGradFraction > 0)
    {
        // Loop over all particles
        for (int m=0; m<_numParticles; m++)
        {
            double totalgrad = 0; // The sum of all individual gradients
            vector<int> neighborID = _mesh->getNeighbors(m);
            int nrNeighbors = neighborID.size();
            Position cellpos = _mesh->particlePosition(m);
            for (int mni=0; mni<nrNeighbors; mni++)
            {
                int mn = neighborID[mni];
                if (mn >= 0) {
                    // gradient with neighbor
                    totalgrad += abs(ds->temperature(mn) - ds->temperature(m)) /
                                 (_mesh->particlePosition(mn) - cellpos).norm();
                }
            }
            VgradTv[m] = ds->volume(m)*totalgrad/nrNeighbors; // Mean gradient
        }
    }
    // Temperature times density (use mass instead of density to draw more from large cells)
    Array TMv; // temperature * mass vector
    TMv.resize(_numParticles);
    if (_tempMassFraction > 0)
    {
        // Loop over all particles
        for (int m=0; m<_numParticles; m++)
        {
            double tempgrad = 0;
            if (_tempGradImportance > 0)
            {
                vector<int> neighborID = _mesh->getNeighbors(m);
                int nrNeighbors = neighborID.size();
                Position cellpos = _mesh->particlePosition(m);
                for (int mni=0; mni < nrNeighbors; mni++)
                {
                    int mn = neighborID[mni];
                    if (mn >= 0)
                    {
                        tempgrad += abs(ds->temperature(mn) - ds->temperature(m)) /
                                    (_mesh->particlePosition(mn) - cellpos).norm();
                    }
                }
                tempgrad /= nrNeighbors;
            }
            TMv[m] = pow(ds->temperature(m), _tempImportance) * pow(ds->density(m), _massImportance) *
                     pow(tempgrad, _tempGradImportance) * ds->volume(m);
        }
    }

    // normalized cdf (will be of length _numParticles+1)
    Array VTcumv; // temperature
    NR::cdf(VTcumv, VTv);
    Array VgradTcumv; // temperature gradient
    NR::cdf(VgradTcumv, VgradTv);
    Array TMcumv; // temperature times density
    NR::cdf(TMcumv, TMv);

    // Write new points to file
    DustGridPlotFile plotPrePoints(this, "ds_pregridpoints");
    DustGridPlotFile plotNewPoints(this, "ds_newgridpoints");
    plotPrePoints.writeLine("#X\tY\tZ\tT\tV\tamount of sampled points\tamount of sampled points per volume");
    plotNewPoints.writeLine("#X\tY\tZ");

    // Save the old points (usually drawn from the dust distribution)
    int oldNumParticles = _numParticles;
    _numParticles = _totalNumParticles; // Increase amount of particles to the total specified in the ski file
    vector<Vec> rv(_numParticles); // Position vector for all the new grid points

    // Sample grid points according to the original distribution (usually dust distribution)
    int dustNumParticles = floor(_numParticles*(1-_tempDistFraction-_tempGradFraction-_tempMassFraction));
    if (dustNumParticles > 0)
    {
        // Assume that we sample along the dust density, even though initial drawing can be different
        DustDistribution* dd = find<DustDistribution>();
        for (int m=0; m<dustNumParticles; m++)
        {
            while (true)
            {
                Position p = dd->generatePosition();
                if (extent().contains(p))       // discard any points outside of the domain
                {
                    rv[m] = p;
                    plotNewPoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z());
                    break;
                }
            }
        }
    }

    Array nrSampledPoints(oldNumParticles); // How many new points sampled in old cell? (write to file)
    // Now pick new points according to the temperature distribution
    int tempDistNumParticles = floor(_numParticles*_tempDistFraction); // Amount of particles drawn from temperature
    int breakIteration = dustNumParticles+tempDistNumParticles; // End of iteration (gets updated)
    for (int m=dustNumParticles; m<breakIteration; m++)
    {
        int cellidx = NR::locate_clip(VTcumv, _random->uniform()); // Determine cell (where we generate new point)
        rv[m] = _mesh->randomPosition(_random, cellidx); // Generate random location in cell
        plotNewPoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z());
        nrSampledPoints[cellidx] = nrSampledPoints[cellidx] + 1;
    }
    // Now pick new points according to the temperature gradient distribution
    int tempGradNumParticles = floor(_numParticles*_tempGradFraction);
    for (int m=breakIteration; m<breakIteration+tempGradNumParticles; m++)
    {
        int cellidx = NR::locate_clip(VgradTcumv, _random->uniform()); // Determine cell (where we generate new point)
        rv[m] = _mesh->randomPosition(_random, cellidx); // Generate random location in cell
        plotNewPoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z());
        nrSampledPoints[cellidx] = nrSampledPoints[cellidx] + 1;
    }
    breakIteration += tempGradNumParticles;
    // Finally, pick according to the temperature times density distribution
    for (int m=breakIteration; m<_numParticles; m++)
    {
        int cellidx = NR::locate_clip(TMcumv, _random->uniform()); // Determine cell (where we generate new point)
        rv[m] = _mesh->randomPosition(_random, cellidx); // Generate random location in cell
        plotNewPoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z());
        nrSampledPoints[cellidx] = nrSampledPoints[cellidx] + 1;
    }
    // Write out interesting properties about the old grid points
    for (int m=0; m<oldNumParticles; m++)
    {
        Array values(4); // T, V and amount of sampled particles in cell
        values[0] = VTv[m]/_mesh->volume(m);
        values[1] = _mesh->volume(m);
        values[2] = nrSampledPoints[m];
        values[3] = nrSampledPoints[m]/_mesh->volume(m);
        Position pos = _mesh->particlePosition(m);
        plotPrePoints.writePoint(pos.x(), pos.y(), pos.z(), values); // Write xyz and T, V and nrsampledpoints
    }
    // With the new particle positions, generate a new voronoi mesh
    log->info("Computing Voronoi tesselation for " + QString::number(dustNumParticles)
                  + " random particles from a dust distribution, "
                  + QString::number(tempDistNumParticles)
                  + " particles distributed according to a temperature distribution, "
                  + QString::number(tempGradNumParticles)
                  + " particles distributed according to a temperature gradient distribution, and "
                  + QString::number(_numParticles-dustNumParticles-tempDistNumParticles-tempGradNumParticles)
                  + " particles distributed according to a temperature times density distribution, "
                  + "with temperature importance " + QString::number(_tempImportance)
                  +  " and mass importance " + QString::number(_massImportance) + ".");
    delete _mesh; // Delete old mesh
    _mesh = new VoronoiMesh(rv, extent(), log, _relaxationSteps);
}

//////////////////////////////////////////////////////////////////////
