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
    if(_tempDistFraction+_tempGradFraction > 0)
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
    if(_tempDistFraction+_tempGradFraction == 0) // If we don't draw from a dynamic distribution, return
        return;
    Log* log = find<Log>();
    log->info("Drawing extra voronoi points from a temperature distribution.");
    DustSystem* ds = find<DustSystem>(); // Cache the dust system

    // Construct cdf of V*T (volume since we want to sample more from bigger cells, if temperature is equal)
    // Start with a V*T vector
    Array VTv; // Volume * Temperature for every cell (_numParticles is the old grid size)
    VTv.resize(_numParticles);

    // Check if we are dealing with a panchromatic or oligochromatic simulation
    bool pan = find<WavelengthGrid>()->issampledrange();
    // Panchromatic: draw from temperature (fill VTv)
    if(pan)
    {
        ds->sumResults(); // To call switchScheme (else there is an error)
        for(int m=0; m<_numParticles; m++)
        {
            // indicative temperature = average population equilibrium temperature weighed by population mass fraction
            const Array& Jv = ds->meanintensityv(m);
            // average over dust components
            double sumRho_h = 0;
            double sumRhoT_h = 0;
            for (int h=0; h<ds->Ncomp(); h++)
            {
                double rho_h = ds->density(m,h);
                if (rho_h>0.0)
                {
                    // average over dust populations within component
                    double sumMu_c = 0;
                    double sumMuT_c = 0;
                    for (int c=0; c<ds->mix(h)->Npop(); c++)
                    {
                        double mu_c = ds->mix(h)->mu(c);
                        double T_c = ds->mix(h)->equilibrium(Jv,c);
                        sumMu_c += mu_c;
                        sumMuT_c += mu_c * T_c;
                    }
                    double T_h = sumMuT_c / sumMu_c;

                    sumRho_h += rho_h;
                    sumRhoT_h += rho_h * T_h;
                }
            }
            double T = sumRhoT_h / sumRho_h;
            VTv[m] = _mesh->volume(m) * T; // V*T vector: volume * temperature
        }
    }
    else // Oligochromatic: use mean intensity instead, for the first wavelength
    {
        ds->meanintensityv(0);
        for(int m=0; m<_numParticles; m++)
        {
            const Array& Jv = ds->meanintensityv(m);
            VTv[m] = _mesh->volume(m) * Jv[0]; // V*T vector: volume * mean intensity (first wavelength)
        }
    }
    // Now for the (normalized) cdf
    Array VTcumv; // Initialize cdf (which will be of length _numParticles+1
    NR::cdf(VTcumv, VTv);

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
            for (int mn=0; mn<nrNeighbors; mn++)
            {
                totalgrad += abs((VTv[mn]/_mesh->volume(mn)) - (VTv[m]/_mesh->volume(m))); // temperature difference
            }
            VgradTv[m] = totalgrad/nrNeighbors; // Mean gradient
        }
    }

    // normalized cdf
    Array VgradTcumv; // Initialize cdf (which will be of length _numParticles+1
    NR::cdf(VgradTcumv, VgradTv);

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
    int dustNumParticles = floor(_numParticles*(1-_tempDistFraction-_tempGradFraction));
    if (dustNumParticles > 0)
    {
        if (_distribution != DustDensity)
            throw FATALERROR("Dynamic grid currently only supports Dust Density distributions!");
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
    for (int m=dustNumParticles; m<dustNumParticles+tempDistNumParticles; m++)
    {
        int cellidx = NR::locate_clip(VTcumv, _random->uniform()); // Determine cell (where we generate new point)
        rv[m] = _mesh->randomPosition(_random, cellidx); // Generate random location in cell
        plotNewPoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z());
        nrSampledPoints[cellidx] = nrSampledPoints[cellidx] + 1;
    }
    // Now pick new points according to the temperature gradient distribution
    for (int m=dustNumParticles+tempDistNumParticles; m<_numParticles; m++)
    {
        int cellidx = NR::locate_clip(VgradTcumv, _random->uniform()); // Determine cell (where we generate new point)
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
        plotPrePoints.writePoint(rv[m].x(), rv[m].y(), rv[m].z(), values); // Write xyz and T, V and nrsampledpoints
    }

    // With the new particle positions, generate a new voronoi mesh
    log->info("Computing Voronoi tesselation for " + QString::number(dustNumParticles)
                  + " random particles from a dust distribution, "
                  + QString::number(tempDistNumParticles)
                  + " particles distributed according to a temperature distribution, and "
                  + QString::number(_numParticles-dustNumParticles-tempDistNumParticles)
                  + " particles distributed according to a temperature gradient distribution.");
    delete _mesh; // Delete old mesh
    _mesh = new VoronoiMesh(rv, extent(), log, _relaxationSteps);
}

//////////////////////////////////////////////////////////////////////
