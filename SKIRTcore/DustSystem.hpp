/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTSYSTEM_HPP
#define DUSTSYSTEM_HPP

#include <mutex>
#include <vector>
#include "Array.hpp"
#include "Position.hpp"
#include "SimulationItem.hpp"
#include "Table.hpp"

class DustDistribution;
class DustGrid;
class DustGridDensityInterface;
class DustMix;
class PhotonPackage;
class ProcessAssigner;

//////////////////////////////////////////////////////////////////////

/** DustSystem is an abstract class for representing complete dust systems, including both a dust
    distribution (i.e. a complete description of the spatial distribution and optical properties of
    the dust) and a grid on which this distribution is discretized. There are specialized
    subclasses for use with oligochromatic and panchromatic simulations respectively. A DustSystem
    object contains a vector of dust cells, each of which contain all useful information on the
    dust within that particular piece of the configuration space. Furthermore, a DustSystem object
    contains pointers a DustDistribution object and to a DustGrid object). A subclass may of course
    maintain additional information depending on its needs. */
class DustSystem : public SimulationItem
{
    Q_OBJECT
    Q_CLASSINFO("Title", "a dust system")

    Q_CLASSINFO("Property", "dustDistribution")
    Q_CLASSINFO("Title", "the dust distribution")
    Q_CLASSINFO("Default", "CompDustDistribution")

    Q_CLASSINFO("Property", "dustGrid")
    Q_CLASSINFO("Title", "the dust grid")
    Q_CLASSINFO("Default", "OctTreeDustGrid")

    Q_CLASSINFO("Property", "sampleCount")
    Q_CLASSINFO("Title", "the number of random density samples for determining cell mass")
    Q_CLASSINFO("MinValue", "10")
    Q_CLASSINFO("MaxValue", "1000")
    Q_CLASSINFO("Default", "100")

    Q_CLASSINFO("Property", "writeConvergence")
    Q_CLASSINFO("Title", "output a data file with convergence checks on the dust system")
    Q_CLASSINFO("Default", "yes")

    Q_CLASSINFO("Property", "writeDensity")
    Q_CLASSINFO("Title", "output FITS files displaying the dust density distribution")
    Q_CLASSINFO("Default", "yes")

    Q_CLASSINFO("Property", "writeDepthMap")
    Q_CLASSINFO("Title", "output FITS file with a V-band optical depth map seen from the center")
    Q_CLASSINFO("Default", "no")

    Q_CLASSINFO("Property", "writeQuality")
    Q_CLASSINFO("Title", "calculate and output quality metrics for the dust grid")
    Q_CLASSINFO("Default", "no")
    Q_CLASSINFO("Silent", "true")

    Q_CLASSINFO("Property", "writeCellProperties")
    Q_CLASSINFO("Title", "output a data file with relevant properties for all dust cells")
    Q_CLASSINFO("Default", "no")
    Q_CLASSINFO("Silent", "true")

    Q_CLASSINFO("Property", "writeCellsCrossed")
    Q_CLASSINFO("Title", "output statistics on the number of cells crossed per path")
    Q_CLASSINFO("Default", "no")
    Q_CLASSINFO("Silent", "true")

    //============= Construction - Setup - Destruction =============

protected:
    /** The default constructor; it is protected since this is an abstract class. */
    DustSystem();

    /** This function verifies that all attribute values have been appropriately set. */
    void setupSelfBefore();

    /** This function performs setup for the dust system, which includes several tasks. First, the
        function verifies that either all dust mixes in the dust system support polarization, or
        none of them do. The next task consists of calculating and storing the volume and the dust
        density (for every dust component) of all the cells. To calculate the volume of a given
        dust cell, we just call the corresponding function of the dust grid. To calculate and set
        the density corresponding to given dust component in a given cell, a number of random
        positions are generated within the cell (see sampleCount()). The density in the cell is
        calculated as the mean of the density values (found using a call to the corresponding
        function of the dust distribution) in these points. The calculation of both volume and
        density is parallellized. Finally, the function optionally invokes various writeXXX()
        functions depending on the state of the corresponding write flags. */
    void setupSelfAfter();

private:
    /** This function serves as the parallelization body for calculating the volume of each cell.
        */
    void setVolumeBody(size_t m);

    /** This function serves as the parallelization body for setting the density value of each cell
        through the DustGridDensityInterface interface, if available. */
    void setGridDensityBody(size_t m);

    /** This function serves as the parallelization body for setting the density value of each cell
        by taking random density sample. */
    void setSampleDensityBody(size_t m);

    /** This function is used to assemble the container that stores the densities of all dust cells
        for each dust component. If multiprocessing is enabled, the calculation of these densities
        can be performed in parallel by the different processes, depending on the type of
        ProcessAssigner that is used for this dust system. When a ProcessAssigner subclass is used
        which distributes the calculation of the dust cell densities amongst the different parallel
        processes, each process contains only the densities for a particular set of dust cells (but
        for each dust component). Therefore, this function uses the PeerToPeerCommunicator object
        to broadcast the densities of the dust cells assigned to a particular process to all other
        processes and storing them in the appropriate place in the container. This is implemented
        as an element-wise summation of this container across all processes, where the density for
        a particular dust cell and dust component will be added to a series of zeros, coming from
        the processes that were not assigned to that dust cell. The end result will be an assembly
        of the densities, where each process stores the density over the entire dust grid. */
    void assemble();

    /** This function writes out a simple text file, named <tt>prefix_ds_convergence.dat</tt>,
        providing a convergence check on the dust system. The function calculates the total dust
        mass, the face-on surface density and the edge-on surface density by directly integrating
        over the dust grid, and compares these values with the expected "theoretical" values found
        by just calling the corresponding functions of the dust distribution. The results are
        written to the file and can be studied to see whether the chosen dust grid is adequate for
        the chosen dust distribution. */
    void writeconvergence() const;

    /** This function writes out FITS files with the theoretical dust density and the
        grid-discretized dust density in the coordinate planes. Each of these maps contains 1024 x
        1024 pixels, and covers as a field of view the total extension of the grid. The number of
        data files written depends on the dimension of the dust system's geometry: for spherical
        symmetry only the intersection with the xy plane is written, for axial symmetry the
        intersections with the xy and xz planes are written, and for general geometries all three
        intersections are written. The difference between the theoretical dust density maps (named
        <tt>prefix_ds_trhoXX.fits</tt>) and the grid-discretized dust density maps (named
        <tt>prefix_ds_grhoXX.fits</tt>) is the following: the theoretical dust density is the total
        dust density of the dust distribution, i.e.\ the actual dust density that would correspond
        to an infinitely fine dust grid. The grid-discretized dust density maps on the other hand
        give the value of the dust as read from the finite-resolution dust grid. A comparison of
        both sets of maps can reveal whether the chosen dust grid is suitable (in the ideal case,
        there would be no difference between both sets of maps). */
    void writedensity() const;

    /** This function writes out a FITS file named <tt>prefix_ds_tau.fits</tt> with an all-sky
        V-band optical depth map as seen from the coordinate origin. The map has 1600 x 800 pixels
        and uses the Mollweide projection to project the complete sky onto a proportional 2:1
        ellipse. The values outside of the ellipse are set to zero. The direction
        \f$\bf{k}=(\theta,\phi)\f$ corresponding to the pixel in the map with horizontal and
        vertical indices \f$(i,j)\f$ can be found through the inverse Mollweide projection, which
        in this case can be written as follows: \f[ x=(i+\frac{1}{2})/N_\mathrm{pixels,x} \f] \f[
        y=(j+\frac{1}{2})/N_\mathrm{pixels,y} \f] \f[ \alpha=\arcsin(2y-1) \f] \f[
        \theta=\arccos(\frac{2\alpha+\sin 2\alpha}{\pi}) \f] \f[ \phi=\frac{\pi(2x-1)}{\cos\alpha}
        \f] */
    void writedepthmap() const;

    /** This function writes out a simple text file, named <tt>prefix_ds_quality.dat</tt>,
        providing some basic quality metrics for the dust grid. The first metric consists of the
        mean value and the standard deviation for the absolute difference \f$|\rho_g-\rho_t|\f$
        between the theoretical and grid density in a large number of randomly chosen points,
        uniformly distributed over the dust grid volume. The second metric consists of the mean
        value and the standard deviation for the difference \f$|\tau_g-\tau_t|\f$ between the
        theoretical and grid optical depth, calculated for a large number of line segments with
        random end points uniformly distributed over the dust grid volume. */
    void writequality() const;

    /** This function writes out a text data file, named <tt>prefix_ds_cellprops.dat</tt>, which
        contains a line for each cell in the dust grid. Each line contains four columns
        representing the following dust cell properties: volume, density, mass fraction and optical
        depth. The first two values are listed in the corresponding output units; the last two are
        dimensionless quantities. */
    void writecellproperties() const;

    //======== Setters & Getters for Discoverable Attributes =======

public:
    /** Sets the dust distribution. */
    Q_INVOKABLE void setDustDistribution(DustDistribution* value);

    /** Returns the dust distribution. */
    Q_INVOKABLE DustDistribution* dustDistribution() const;

    /** Sets the dust grid. */
    Q_INVOKABLE void setDustGrid(DustGrid* value);

    /** Returns the dust grid. See also grid(). */
    Q_INVOKABLE DustGrid* dustGrid() const;

    /** Sets the number of random positions on which the density is sampled for each cell in the
        dust grid. The default value is 100 samples per cell. */
    Q_INVOKABLE void setSampleCount(int value);

    /** Returns the number of random positions on which the density is sampled for each cell in the
        dust grid. */
    Q_INVOKABLE int sampleCount() const;

    /** Sets the flag that indicates whether or not to output a data file with convergence checks
        on the dust system. The default value is true. */
    Q_INVOKABLE void setWriteConvergence(bool value);

    /** Returns the flag that indicates whether or not to output a data file with convergence
        checks on the dust system. */
    Q_INVOKABLE bool writeConvergence() const;

    /** Sets the flag that indicates whether or not to output FITS files displaying the dust
        density distribution for further analysis. The default value is true. */
    Q_INVOKABLE void setWriteDensity(bool value);

    /** Returns the flag that indicates whether or not to output FITS files displaying the dust
        density distribution for further analysis. */
    Q_INVOKABLE bool writeDensity() const;

    /** Sets the flag that indicates whether or not to output a FITS file with an all-sky V-band
        optical depth map as seen from the center. The default value is false. */
    Q_INVOKABLE void setWriteDepthMap(bool value);

    /** Returns the flag that indicates whether or not to output a FITS file with an all-sky V-band
        optical depth map as seen from the center. */
    Q_INVOKABLE bool writeDepthMap() const;

    /** Sets the flag that indicates whether or not to calculate and output quality metrics for the
        dust grid. The default value is false. */
    Q_INVOKABLE void setWriteQuality(bool value);

    /** Returns the flag that indicates whether or not to calculate and output quality metrics for
        the dust grid. */
    Q_INVOKABLE bool writeQuality() const;

    /** Sets the flag that indicates whether or not to output a data file with relevant properties
        for all dust cells, including volume, density, mass fraction and optical depth. The default
        value is false. */
    Q_INVOKABLE void setWriteCellProperties(bool value);

    /** Returns the flag that indicates whether or not to output a data file with relevant
        properties for all dust cells. */
    Q_INVOKABLE bool writeCellProperties() const;

    /** Sets the flag that indicates whether or not to output a text file with statistics on the
        number of dust grid cells crossed per path calculated through the grid. The default value
        is false. */
    Q_INVOKABLE void setWriteCellsCrossed(bool value);

    /** Returns the flag that indicates whether or not to output a text file with statistics on the
        number of dust grid cells crossed per path calculated through the grid. */
    Q_INVOKABLE bool writeCellsCrossed() const;

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the dust system, which depends on the (lack of)
        symmetry in the geometry of its distribution. A value of 1 means spherical symmetry, 2
        means axial symmetry and 3 means none of these symmetries. */
    int dimension() const;

    /** This function returns the number of dust cells. */
    int Ncells() const;

    /** This function returns the number of dust components. */
    int Ncomp() const;

    /** This function returns true if the dust mixes in this dust system support polarization;
        false otherwise. During setup it is verified that either all dust mixes support
        polarization, or none of them do. */
    bool polarization() const;

    /** This function returns a pointer to the dust mixture corresponding to the \f$h\f$'th dust
        component. */
    DustMix* mix(int h) const;

    /** This function returns a pointer to a dust mixture that is selected randomly among the dust
        mixes of the dust components in the dust system. If we have just a single dust component,
        this is simple. If there are multiple dust components, the relative probability of
        selecting the dust mix for a particular dust component \f$h\f$ is given by
        \f$\kappa_\ell^{\text{sca}}(h)\,\rho({\bf{r}},h)\f$ at wavelength index \f$\ell\f$. */
    DustMix* randomMixForPosition(Position bfr, int ell) const;

    /** This function returns the number of the dust cell that contains the position
        \f${\boldsymbol{r}}\f$. The function just passes the call to corresponding function of the
        dust grid. */
    int whichcell(Position bfr) const;

    /** This function returns a random location in the dust cell with cell number \f$m\f$. The
        function just passes the call to corresponding function of the dust grid. */
    Position randomPositionInCell(int m) const;

    /** This function returns the volume of the dust cell with cell number \f$m\f$. */
    double volume(int m) const;

    /** This function returns the dust density corresponding to dust component \f$h\f$ of the dust
        cell with cell number \f$m\f$. If \f$m=-1\f$, i.e. if the cell number corresponds to a
        non-existing cell outside the grid, the value zero is returned. */
    double density(int m, int h) const;

    /** This function returns the total dust density of the dust cell with cell number \f$m\f$. If
        \f$m=-1\f$, i.e. if the cell number corresponds to a non-existing cell outside the grid,
        the value zero is returned. */
    double density(int m) const;

    /** This function calculates the optical depth
        \f$\tau_{\ell,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}})\f$ at wavelength index
        \f$\ell\f$ along a path through the dust system starting at the position
        \f${\boldsymbol{r}}\f$ into the direction \f${\boldsymbol{k}}\f$, where \f$\ell\f$,
        \f${\boldsymbol{r}}\f$ and \f${\boldsymbol{k}}\f$ are obtained from the specified
        PhotonPackage object, and it stores the resulting details back into the photon package
        object. The hard work is done by calling the DustGrid::path() function which stores the
        geometrical information on the path through the dust grid into the photon package: the cell
        numbers \f$m\f$ of the cells that are crossed by the path, the pathlength \f$(\Delta
        s)_m\f$ covered in that particular cell and a total path length counter \f$s_m\f$ that
        gives the total path length covered between the starting point \f${\boldsymbol{r}}\f$ and
        the boundary of the cell. With this information given, the calculation of the optical depth
        is rather straightforward: it is calculated as \f[
        \tau_{\ell,{\text{path}}}({\boldsymbol{r}},{\boldsymbol{k}}) = \sum_m (\Delta s)_m \sum_h
        \kappa_{\ell,h}^{\text{ext}}\, \rho_m, \f] where \f$\kappa_{\ell,h}^{\text{abs}}\f$ is the
        extinction coefficient corresponding to the \f$h\f$'th dust component at wavelength index
        \f$\ell\f$ and \f$\rho_{m,h}\f$ the dust density in the cell with cell number \f$m\f$
        corresponding to the \f$h\f$'th dust component. The function also stores the details on the
        calculation of the optical depth in the photon package, specifically it stores the optical
        depth covered within the \f$m\f$'th dust cell, \f[ (\Delta\tau_\ell)_m = (\Delta s)_m
        \sum_h \kappa_{\ell,h}^{\text{ext}}\, \rho_m, \f] and the total optical depth
        \f$\tau_{\ell,m}\f$ covered between the starting point \f${\boldsymbol{r}}\f$ and the
        boundary of the cell. */
    void fillOpticalDepth(PhotonPackage* pp);

    /** This function returns the optical depth
        \f$\tau_{\ell,{\text{d}}}({\boldsymbol{r}},{\boldsymbol{k}})\f$ at wavelength index
        \f$\ell\f$ along a path through the dust system starting at the position
        \f${\boldsymbol{r}}\f$ into the direction \f${\boldsymbol{k}}\f$ for a distance \f$d\f$,
        where \f$\ell\f$, \f${\boldsymbol{r}}\f$ and \f${\boldsymbol{k}}\f$ are obtained from the
        specified PhotonPackage object. The function first determines the photon package's path
        through the dust grid, storing the geometric information about the path segments through
        each cell into the photon package, and then calculates the optical depth at the specified
        distance. The calculation proceeds as described for the fillOpticalDepth() function; the
        differences being that the path length is limited to the specified distance, and that this
        function does not store the optical depth information back into the PhotonPackage object.
        */
    double opticaldepth(PhotonPackage* pp, double distance);

    /** If the writeCellsCrossed attribute is true, this function writes out a data file (named
        <tt>prefix_ds_crossed.dat</tt>) with statistics on the number of dust grid cells crossed
        per path calculated through the grid. The first column on each line specifies a particular
        number of cells crossed; the second column indicates the number of paths that crossed this
        precise number of cells. In effect this provides a histogram for the distribution of the
        path length (measured in the number of cells crossed). This virtual function can be
        overridden in a subclass to write out additional results of the simulation stored in the
        dust system. In that case, the overriding function must also call the implementation in
        this base class. The filenameSuffix gets added to the end of the filename. */
    virtual void write(QString filenameSuffix = "") const;

    /** This pure virtual function must be implemented in each subclass to indicate whether dust
        emission is turned on for this dust system. The function returns true if dust emission is
        turned on, and false otherwise. It is provided in this base class because it is invoked
        from the general MonteCarloSimulation class. */
    virtual bool dustemission() const = 0;

    /** This function synchronizes the results of the absorption by calling the sync() function on the
        absorption tables. **/
    virtual void sumResults();

    /** This pure virtual function must be implemented in each subclass to indicate whether the
        absorption rates in each cell need to be stored for this dust system. This is needed if
        dust emission is turned on, but it can also be chosen if the user wants to study the mean
        intensity of the radiation field. It is provided in this base class because it is invoked
        from the general MonteCarloSimulation class. */
    virtual bool storeabsorptionrates() const = 0;

    /** This pure virtual function must be implemented in each subclass to simulate absorption of
        of a monochromatic luminosity package in the specified dust cell. The function should be
        called only if dustemission() returns true. It is provided in this base class because it is
        referenced from the general MonteCarloSimulation class (although it is actually invoked
        only for panchromatic simulations). */
    virtual void absorb(int m, int ell, double DeltaL, bool ynstellar) = 0;

    /** This pure virtual function returns the absorbed luminosity \f$L_{\ell,m}\f$ at wavelength index
        \f$\ell\f$ in the dust cell with cell number \f$m\f$. */
    virtual double Labs(int m, int ell) const = 0;

    /** This function returns a vector with the mean radiation field \f$J_{\ell,m}\f$ at all
        wavelength indices in the dust cell with cell number \f$m\f$. It is calculated as \f[
        J_{\ell,m} = \frac{ L_{\ell,m}^{\text{abs}} }{ 4\pi\, V_m\, (\Delta\lambda)_\ell \sum_h
        \kappa_{\ell,h}^{\text{abs}}\, \rho_{m,h} } \f] with \f$L_{\ell,m}^{\text{abs}}\f$ the
        absorbed luminosity, \f$\kappa_{\ell,h}^{\text{abs}}\f$ the absorption coefficient
        corresponding to the \f$h\f$'th dust component, \f$\rho_{m,h}\f$ the dust density
        corresponding to the \f$h\f$'th dust component, and \f$V_m\f$ the volume of the cell. */
    Array meanintensityv(int m) const;

    /** This function reinitialises the grid, used when building a dynamic grid */
    virtual void reinitialiseGrid() = 0;

    //======================== Data Members ========================

protected:
    // data members to be set before setup is invoked
    DustDistribution* _dd;
    DustGrid* _grid;
    DustGridDensityInterface* _gdi;
    int _Nrandom;
    bool _writeConvergence;
    bool _writeDensity;
    bool _writeDepthMap;
    bool _writeQuality;
    bool _writeCellProperties;
    bool _writeCellsCrossed;

    // data members initialized during setup
    ProcessAssigner* _setupAssigner;  // determines which dust cells are assigned to this process
                                      // for various calculations during setup
    int _Ncomp;
    int _Ncells;
    Array _volumev;     // volume for each cell (indexed on m)
    Table<2> _rhovv;    // density for each cell and each dust component (indexed on m,h)
    std::vector<qint64> _crossed;
    std::mutex _crossedMutex;
};

//////////////////////////////////////////////////////////////////////

#endif // DUSTSYSTEM_HPP
