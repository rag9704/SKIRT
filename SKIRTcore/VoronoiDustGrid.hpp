/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIDUSTGRID_HPP
#define VORONOIDUSTGRID_HPP

#include "BoxDustGrid.hpp"
#include "Random.hpp"
class DustDistribution;
class VoronoiMesh;
class VoronoiMeshFile;

//////////////////////////////////////////////////////////////////////

/** VoronoiDustGrid is a concrete subclass of the BoxDustGrid class. It is used for representing a
    three-dimensional dust grid based on a Voronoi tesselation of the cuboid containing
    substantially all of the dust. See the VoronoiMesh class for more information on Voronoi
    tesselations. The class offers several options for determining the locations of the particles
    generating the Voronoi tesselation. A specified number of particles can be distributed randomly
    over the domain, either uniformly or with the same overall density distribution as the dust.
    Alternatively, the locations can be copied from the particles in an SPH dust distribution. This
    class uses the Voro++ code written by Chris H. Rycroft (LBL / UC Berkeley) to generate output
    files for plotting the Voronoi grid. */
class VoronoiDustGrid : public BoxDustGrid
{
    Q_OBJECT
    Q_CLASSINFO("Title", "a Voronoi dust grid")

    Q_CLASSINFO("Property", "numParticles")
    Q_CLASSINFO("Title", "the number of random particles (or cells in the grid)")
    Q_CLASSINFO("MinValue", "10")
    Q_CLASSINFO("Default", "500")

    Q_CLASSINFO("Property", "distribution")
    Q_CLASSINFO("Title", "the probablity distribution for the particles")
    Q_CLASSINFO("Uniform", "uniform")
    Q_CLASSINFO("CentralPeak", "with a steep central peak")
    Q_CLASSINFO("DustDensity", "same as dust density")
    Q_CLASSINFO("DustTesselation", "copy Voronoi tesselation in dust distribution")
    Q_CLASSINFO("SPHParticles", "exact locations of SPH particles in dust distribution")
    Q_CLASSINFO("File", "particle locations in a specified data file")
    Q_CLASSINFO("Default", "DustDensity")
    Q_CLASSINFO("TrueIf", "File")

    Q_CLASSINFO("Property", "relaxationSteps")
    Q_CLASSINFO("Title", "the amount of times the voronoi cells will relax")
    Q_CLASSINFO("MinValue", "0")
    Q_CLASSINFO("Default", "0")

    Q_CLASSINFO("Property", "preGridPointFraction")
    Q_CLASSINFO("Title", "the percentage of grid particles used for the original grid (when using dynamic grids)")
    Q_CLASSINFO("MinValue", "0.0001")
    Q_CLASSINFO("MaxValue", "1")
    Q_CLASSINFO("Default", "0.8")

    Q_CLASSINFO("Property", "tempDistFraction")
    Q_CLASSINFO("Title", "the percentage of grid particles drawn from a temperature distribution")
    Q_CLASSINFO("MinValue", "0")
    Q_CLASSINFO("MaxValue", "1")
    Q_CLASSINFO("Default", "0")

    Q_CLASSINFO("Property", "tempGradFraction")
    Q_CLASSINFO("Title", "the percentage of grid particles drawn from a temperature gradient distribution")
    Q_CLASSINFO("MinValue", "0")
    Q_CLASSINFO("MaxValue", "1")
    Q_CLASSINFO("Default", "0")

    Q_CLASSINFO("Property", "tempDensityFraction")
    Q_CLASSINFO("Title", "the percentage of grid particles drawn from a temperature times density distribution")
    Q_CLASSINFO("MinValue", "0")
    Q_CLASSINFO("MaxValue", "1")
    Q_CLASSINFO("Default", "0")

    Q_CLASSINFO("Property", "densityImportance")
    Q_CLASSINFO("Title", "the importance of density when using tempDensityFraction")
    Q_CLASSINFO("MinValue", "0")
    Q_CLASSINFO("MaxValue", "100")
    Q_CLASSINFO("Default", "1")

    Q_CLASSINFO("Property", "voronoiMeshFile")
    Q_CLASSINFO("Title", "the Voronoi mesh data file")
    Q_CLASSINFO("Optional", "true")
    Q_CLASSINFO("Default", "VoronoiMeshAsciiFile")
    Q_CLASSINFO("RelevantIf", "distribution")

    //============= Construction - Setup - Destruction =============

public:
    /** The default constructor. */
    Q_INVOKABLE VoronoiDustGrid();

    /** The destructor cleans up all memory allocated in this class. */
    ~VoronoiDustGrid();

protected:
    /** This function verifies that the attributes have been appropriately set, selects the
        requested particles for generating the Voronoi tesselation, and finally constructs the
        Voronoi tesselation through an instance of the VoronoiMesh class. If requested, it also
        outputs files that can be used for plotting the grid. To generate the random particles
        according to the simulation's dust density distribution, the setup function partitions the
        domain using a three-dimensional cuboidal cell structure, called the foam. The distribution
        of the foam cells is determined automatically from the dust density distribution. The foam
        allows to efficiently generate random points drawn from this probability distribution. Once
        the particles have been generated, the foam is discarded. */
    void setupSelfBefore();

    /** This function makes sure that we only use the meanintensity when writemeanintensity is set to
        true. */
    void setupSelfAfter();

    //======== Setters & Getters for Discoverable Attributes =======

public:
    /** Sets the number of random particles (or cells in the grid). */
    Q_INVOKABLE void setNumParticles(int value);

    /** Returns number of random particles (or cells in the grid). */
    Q_INVOKABLE int numParticles() const;

    /** The enumeration type indicating the probability distribution used for generating the random
        particles. */
    Q_ENUMS(Distribution)
    enum Distribution { Uniform, CentralPeak, DustDensity, DustTesselation, SPHParticles, File };

    /** Sets the enumeration value indicating the probability distribution used for generating the
        random particles. */
    Q_INVOKABLE void setDistribution(Distribution value);

    /** Returns the enumeration value indicating the probability distribution used for generating
        the random particles. */
    Q_INVOKABLE Distribution distribution() const;

    /** Sets the amount of times Lloyd's algorithm will be used to relaxate the voronoi grid */
    Q_INVOKABLE void setRelaxationSteps(int value);

    /** Returns the amound of times Lloyd's algorithm will be used to relaxate the voronoi grid */
    Q_INVOKABLE int relaxationSteps() const;

    /** Sets the fraction of grid particles drawn from a temperature distribution */
    Q_INVOKABLE void setPreGridPointFraction(double value);

    /** Returns the fraction of grid particles drawn from a temperature distribution */
    Q_INVOKABLE double preGridPointFraction() const;

    /** Sets the fraction of grid particles drawn from a temperature distribution */
    Q_INVOKABLE void setTempDistFraction(double value);

    /** Returns the fraction of grid particles drawn from a temperature distribution */
    Q_INVOKABLE double tempDistFraction() const;

    /** Sets the fraction of grid particles drawn from a temperature gradient distribution */
    Q_INVOKABLE void setTempGradFraction(double value);

    /** Returns the fraction of grid particles drawn from a temperature gradient distribution */
    Q_INVOKABLE double tempGradFraction() const;

    /** Sets the fraction of grid particles drawn from a temperature times density distribution */
    Q_INVOKABLE void setTempDensityFraction(double value);

    /** Returns the fraction of grid particles drawn from a temperature times density distribution */
    Q_INVOKABLE double tempDensityFraction() const;

    /** Sets the importance of density when using a temperature times density distribution */
    Q_INVOKABLE void setDensityImportance(double value);

    /** Returns the importance of density when using a temperature times density distribution */
    Q_INVOKABLE double densityImportance() const;

    /** Sets the file containing the Voronoi particle locations in case \em distribution has the
        value \em File. */
    Q_INVOKABLE void setVoronoiMeshFile(VoronoiMeshFile* value);

    /** Returns the file containing the Voronoi particle locations in case \em distribution has the
        value \em File. */
    Q_INVOKABLE VoronoiMeshFile* voronoiMeshFile() const;

    //======================== Other Functions =======================

public:
    /** This function returns the volume of the dust cell with cell number \f$m\f$. */
    double volume(int m) const;

    /** This function returns the number of cells in the dust grid. */
    int numCells() const;

    /** This function returns the number of the dust cell that contains the position
        \f${\bf{r}}\f$. See the VoronoiMesh class for more information. */
    int whichcell(Position bfr) const;

    /** This function returns the central location of the dust cell with cell number \f$m\f$. In
        this class the function returns the centroid of the Voronoi cell. */
    Position centralPositionInCell(int m) const;

    /** This function returns a random location from the dust cell with cell number \f$m\f$. */
    Position randomPositionInCell(int m) const;

    /** This function calculates a path through the grid. The DustGridPath object passed as an
        argument specifies the starting position \f${\bf{r}}\f$ and the direction \f${\bf{k}}\f$
        for the path. The data on the calculated path are added back into the same object. See the
        VoronoiMesh class for more information. */
    void path(DustGridPath* path) const;

    /** This function draws extra grid points from a temperature distribution */
    void drawFromTemperatureDistribution();

    //======================== Data Members ========================

private:
    // discoverable attributes (in addition to extent which is stored in inherited Box)
    int _numParticles;
    Distribution _distribution;
    int _relaxationSteps;
    double _preGridPointFraction;
    double _tempDistFraction;
    double _tempGradFraction;
    double _tempDensityFraction;
    double _densityImportance;
    VoronoiMeshFile* _meshfile;

    // data members initialized during setup
    Random* _random;
    VoronoiMesh* _mesh;
    bool _meshOwned;
    int _totalNumParticles; // Total amount of particles (current can be lower during prepackage phase)
};

//////////////////////////////////////////////////////////////////////

#endif // VORONOIDUSTGRID_HPP
