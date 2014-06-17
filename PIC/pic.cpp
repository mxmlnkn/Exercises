// del pic.exe && g++ pic.cpp -o pic.exe -Wall -fopenmp && pic.exe
/******************************************************************************
 * ToDo:                                                                      *
 *   - use cuda                                                               *
 *   - implement weighting                                                    *
 *   - make it work for more than one cell (only force of particles in same cell)
 *   - weighting ~ m ~ 1/a ~ dt_min_needed ?!                                 *
 *   - copy relativistic formula in here                                      *
 *   - copy exact (relativistic) pusher from picongpu to Verlet               *
 *   - using more particles see them taking on maxwell distribution (fit)     *
 *   - make work for 3D                                                       *
 *   - seed with maxwell temperature                                          *
 *   - ensure momentum conservation for reflecting boundaries in one cell     *
 *     by also giving the simulation box a velocity !                         *
 *   - look why momentum rises, but kinetic energy is constant ...            *
 *       P^2 = (p1+p2+p3+..)^2 !~ Ekin ~ p1^2+p2^2+p3^2+.. (difference 2p1p2) *
 *   - use also protons => will that attenuate the error for energy           *
 *     conservation by itself?                                                *
 *   - when calculating forces of neighbors, sum up beginning with smallest   *
 *     force to minimize numeric roundoff errors                              *
 *   - use particle lists per cell -> implement cell change in CheckBoundary  *
 *   - F(r), V(r) output for all forces and for horizintal and diagonal dist. *
 *   - cellwise lists to make more cells only linearly more demanding         *
 *   - output electrons and ions separatedly                                  *
 *   - calculate force in 3 different steps (clear, calc, apply)              *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <omp.h>
#include <iostream>
#include <fstream>

#define DEBUG 1

using namespace std;

#include "Vector.cpp"
#include "Parameters.cpp"

class teeStream{
    public:
        ofstream fileStream;
        teeStream( const char * filename ) {
            fileStream.open( filename, std::ofstream::out | std::ofstream::app );
        }
        ~teeStream(void) {
            fileStream.close();
        }
        template <class T> teeStream& operator<< (T val) {
            fileStream << val;
            cout << val;
            fileStream.flush();
            return *this;
        }
        teeStream& operator<< (ostream& (*pfun)(ostream&)) {
            pfun(fileStream);
            pfun(cout);
            return *this;
        }
};

teeStream tout("output.txt");

const uint16_t STEPS_TO_REMEMBER = 2;
struct Particle {
    // previous steps needed for example for relativistic formulas because of time retardation
    Vec r,p;
    Vec r_prev[STEPS_TO_REMEMBER],p_prev[STEPS_TO_REMEMBER];
    Vec dr,dp;
    double m,q;
};

inline Vec ForceActingOnParticle1( const Particle & particle1, const Particle & particle2, uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    // if not in same cell, then force is 0
    if (!FULL_N_SQUARED_FORCE)
        // may not round to 3.0, but 3.00...1, but should still round to one and the same 3.00...1, so that they are comparable. If not, the C-standard and floor
        if ( floor( particle1.r.x / CELL_SIZE_X ) != floor( particle2.r.x / CELL_SIZE_X )
          or floor( particle1.r.y / CELL_SIZE_Y ) != floor( particle2.r.y / CELL_SIZE_Y )
          or floor( particle1.r.z / CELL_SIZE_Z ) != floor( particle2.r.z / CELL_SIZE_Z ) ) {
            if ( NUMBER_OF_CELLS_X * NUMBER_OF_CELLS_Y * NUMBER_OF_CELLS_Z == 1 )
                tout << "ASSERT ERROR: Found particle outside of cell when calculating force, even though there is only 1 cell.\n";
            return Vec(0,0,0);
          }

    /* Force felt by Particle 1, because it moves in the field of Particle 2 *
     * If sgn(q1) = sgn(q2 ) and x1 = 0, x2 > 0 => force should point to     *
     * negative x-axis direction. r12 = x1-x2 < 0, rest is positive => true  */
    /* shape types:
     * 00 - point-point                                                      *
     * 01 - ball-ball (CIC radialsymmetric equivalent)                       *
     * 02 - TSC equivalent (ball with linearly falling charge distr.)        *
     * 03 -                                                                  *
     * 04 - 4th order                                                        *
     * ..                                                                    *
     * 11 - forces between two rectangles                                    *
     * 12 - correct rectangular TSC shape                                    *
     * ..                                                                    *
     * 99 - sphere-sphere                                                    */

    const double cloudRadius = CELL_SIZE_X / 2.0;
    const double alpha       = particle1.q * particle2.q / ( 4.*M_PI*EPS0 );
    Vec F(0,0,0);
    
    // Loop over periodically continued simulation box i.e. Supercells, because this simulation only supports one supercell at a time with not so clearly separated cells
    const uint16_t supercells_to_consider_X = int( ceil( CELL_SIZE_X * CONSIDERATION_RATIO / ( CELL_SIZE_X * CELL_SIZE_X ) ) + 0.1 );
    const uint16_t supercells_to_consider_Y = int( ceil( CELL_SIZE_Y * CONSIDERATION_RATIO / ( CELL_SIZE_Y * CELL_SIZE_Y ) ) + 0.1 );
    const uint16_t supercells_to_consider_Z = int( ceil( CELL_SIZE_Z * CONSIDERATION_RATIO / ( CELL_SIZE_Z * CELL_SIZE_Z ) ) + 0.1 );
    const uint16_t supercells_to_consider   = max( supercells_to_consider_X, max( supercells_to_consider_Y, supercells_to_consider_Z ) );
    
    // if maxXYZ == 0, then the for loop is rendered basically useless
    const uint16_t maxX = (SIMDIM >= 1) ? supercells_to_consider : 0;
    const uint16_t maxY = (SIMDIM >= 2) ? supercells_to_consider : 0;
    const uint16_t maxZ = (SIMDIM >= 3) ? supercells_to_consider : 0;
    
    for (int iPeriodZ = -maxZ; iPeriodZ <= maxZ; iPeriodZ++)
    for (int iPeriodY = -maxY; iPeriodY <= maxY; iPeriodY++)
    for (int iPeriodX = -maxX; iPeriodX <= maxX; iPeriodX++) {
        // this is the only line changing in the loops 
        const Vec r12  = particle1.r - ( particle2.r + Vec( iPeriodX * CELL_SIZE_X, iPeriodY * CELL_SIZE_Y, iPeriodZ * CELL_SIZE_Z) );
        const double r = r12.norm();
        Vec F12 = alpha * r12/r;
        
        if( r < CELL_SIZE_MIN * CONSIDERATION_RATIO ) {
            if (colshape == 99) {
                if (r < 2.*cloudRadius)
                    F12 = F12 / ( 4.*cloudRadius*cloudRadius );
                else
                    F12 = F12 / ( r*r );
            } else if (colshape == 0) {
                F12 = F12 / ( r*r );
            } else if( colshape == 1) {
                if (r < 2.*cloudRadius) {
                    const double radii_ratio = r/cloudRadius;
                    // will result in F12 = F12/(4*cloudRadius^2) for radii_ratio = 2 :)
                    F12 = F12 * double( (1./32.)*pow(radii_ratio,4) - (9./16.)*radii_ratio*radii_ratio + radii_ratio ) / ( cloudRadius*cloudRadius );
                } else
                    F12 = F12 / ( r*r );
            }
        } else
            F12 *= 0;
        
        F += F12;
    }
    
    return F;
}

inline double Energy( const Particle & particle1, const Particle & particle2, uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    const double cloudRadius = CELL_SIZE_X / 2.0;
    const Vec r12            = particle1.r - particle2.r;
    const double r           = r12.norm();
    // for two positive charges the potential energy is supposed to increase for smaller distances, meaning the particles slow down
    double V = particle1.q * particle2.q / ( 4.*M_PI*EPS0 );

    if (colshape == 99) {
        if (r < 2*cloudRadius)
            // r=2*cloudRadius => V = 1/(2*cloudRadius) => continuous
            V *= 1./cloudRadius - r/(4.*cloudRadius*cloudRadius);
        else
            V = V / r;
    } else if (colshape == 0)
        V = V / r;
    else if( colshape == 1) {
        const double radii_ratio = r / cloudRadius;
        if (r < 2*cloudRadius)
            // r=2*cloudRadius => radii_ratio=2 => V = 1/(2*cloudRadius) (correct :) )
            V *= (192. - 80. * pow(radii_ratio,2) + 30 * pow(radii_ratio,3) - pow(radii_ratio,5) ) / (160.*cloudRadius);
        else
            V = V / r;
    }
    return V;
}

/* checks if particles has left the simulation space. Should be called every  *
 * time the position has been changed. Therefore would be nice to implement   *
 * this in all positions change routines, but then Vector.hpp wouldn't be     *
 * general anymore                                                            *
 * If some particles has left, then put it back to the wall in that dimension *
 * and clear the momentum of that direction like a infinitely small elastic   *
 * coefficient                                                                */
bool CheckBoundaryConditions( Particle & particle ) {
    uint16_t outOfBounds = 0;
    // Cycle through all (max. 3) dimensions
    for (uint16_t dim = 0; dim < SIMDIM; dim++) {
        if ( (particle.r[dim] < 0) or (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) )
            outOfBounds++;

        if (BOUNDARY_CONDITION == 0) {
            /* Periodic Boundary Conditions */
            if (particle.r[dim] < 0) {
                particle.r[dim] += NUMBER_OF_CELLS[dim] * CELL_SIZE[dim];
            } else if (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
                particle.r[dim] -= NUMBER_OF_CELLS[dim] * CELL_SIZE[dim]; //not safe for two cells at one push ! => courant criterium !
            }
        }
        
        else if (BOUNDARY_CONDITION == 1) {
            /* Reflecting Boundary Condition */
            if (particle.r[dim] < 0) {
                particle.p[dim] *= -1;
                particle.r[dim] = 0 + ( 0 - particle.r[dim] );
            } else if ( particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
                particle.p[dim] *= -1;
                particle.r[dim] = NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] - (particle.r[dim] - NUMBER_OF_CELLS[dim] * CELL_SIZE[dim]);
            }
        }
        
        else if (BOUNDARY_CONDITION == 2) {
            /* Adhering Boundary Condition */
            if (particle.r[dim] < 0) {
                particle.r[dim] = 0;
                particle.p[dim] = 0;
            } else if (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
                particle.r[dim] = 0;
                particle.p[dim] = 0;
            }
        }
    }
    return outOfBounds > 0;
}

void Verlet( const uint32_t NUMBER_OF_PARTICLES, Particle particle[], const double dt, bool initialize = false ) {
    /* Algorithm:                                   *
     *     p += ForceActingOnThisParticle * dt/2.0; *
     *     r += p / m * dt/2.0;                     *
     *     r += p / m * dt/2.0;                     *
     *     p += ForceActingOnThisParticle * dt/2.0; */

    // Clear dp
    #pragma omp parallel for
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        particle[i].dp = Vec(0,0,0);

	// Calculate dp
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        //Calculate Net Force on Particle i
        for (uint32_t j = 0; j < i; j++) {
            // this part relativistically ...
            const Vec dp12 = ForceActingOnParticle1( particle[i], particle[j] ) * DELTA_T;
            particle[i].dp += dp12;
            particle[j].dp -= dp12;    // only possible if (i,j) and (j,i) aren't both calculated => j <= i condition, and j != i, because that is self collision
        }

    // Apply Momentum Change and push
    #pragma omp parallel for
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
        // r_prev[highest] is the oldest position => shift r_prev[i] to r_prev[i+1]
        for (uint16_t j = STEPS_TO_REMEMBER-1; j > 0; j--) {
            particle[i].r_prev[j] = particle[i].r_prev[j-1];
            particle[i].p_prev[j] = particle[i].p_prev[j-1];
        }
        particle[i].r_prev[0] = particle[i].r;
        particle[i].p_prev[0] = particle[i].p;
        if (initialize) {
            // normally p += dp and r += dp/m*dt, but this is the initialization
            particle[i].p += particle[i].dp / 2.0;
        } else {
            particle[i].p += particle[i].dp / 2.0;
            particle[i].r += particle[i].p / particle[i].m * dt;
            CheckBoundaryConditions(particle[i]);
        }
    }

    if (initialize)
        return;

    // Clear dp
    #pragma omp parallel for
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        particle[i].dp = Vec(0,0,0);

	// Calculate dp again
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        //Calculate Net Force on Particle i
        for (uint32_t j = 0; j < i; j++) {
            if (i==j) continue;
            // this part relativistically ...
            const Vec dp12 = ForceActingOnParticle1( particle[i], particle[j] ) * DELTA_T;
            particle[i].dp += dp12;
            particle[j].dp -= dp12;    // only possible if (i,j) and (j,i) aren't both calculated => j <= i condition, and j != i, because that is self collision
        }

    // Only Apply new Momentum Change
    #pragma omp parallel for
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
        particle[i].p += particle[i].dp / 2.0;
        //particle[i].p *= 0.9999;  // simulate stopping (make dependent of time step...)
    }
    // no shift of previous positions here, because this is only a step for calculating the current momentum !

    return;
}

double CalcTotalKineticEnergy( const unsigned int particleCount, Particle particle[] ) {
    double T = 0;
    for (unsigned int i = 0; i < particleCount; i++)
        T += particle[i].p.norm2() / (2.0 * particle[i].m);
    return T;
}

double CalcTotalPotentialEnergy( const unsigned int particleCount, Particle particle[], uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    double E = 0;
    #pragma omp parallel for
    for (unsigned int i = 0; i < particleCount; i++) {
        double V = 0;
        for (unsigned int j = 0; j < i; j++) {
            // only if in same cell
            if ( floor( particle[i].r.x / CELL_SIZE_X ) == floor( particle[j].r.x / CELL_SIZE_X )
             and floor( particle[i].r.y / CELL_SIZE_Y ) == floor( particle[j].r.y / CELL_SIZE_Y )
             and floor( particle[i].r.z / CELL_SIZE_Z ) == floor( particle[j].r.z / CELL_SIZE_Z ) )
                V += Energy( particle[i], particle[j], colshape );
        }
        #pragma omp atomic
        E += V;
    }
    return E;
}

double CalcTotalAngularMomentum( const unsigned int particleCount, Particle particle[] ) {
    Vec L;
    for (unsigned int i = 0; i < particleCount; i++) {
        L.x += particle[i].r.y * particle[i].p.z - particle[i].r.z * particle[i].p.y;
        L.y += particle[i].r.z * particle[i].p.x - particle[i].r.x * particle[i].p.z;
        L.z += particle[i].r.x * particle[i].p.y - particle[i].r.y * particle[i].p.x;
    }
    return L.norm();
}

double CalcTotalMomentum( const unsigned int particleCount, Particle particle[] ) {
    Vec P;
    for (unsigned int i = 0; i < particleCount; i++) {
        P.x += particle[i].p.x;
        P.y += particle[i].p.y;
        P.z += particle[i].p.z;
    }
    return P.norm();
}


template<typename T_Datatype>
class DataBinPlot {
    T_Datatype min, max;
    uint16_t Nbins;
    T_Datatype * data;
    const char * filename;
    bool writtenToFile;

public:
    DataBinPlot( T_Datatype min, T_Datatype max, uint16_t Nbins, const char * filename) : min(min), max(max), Nbins(Nbins), filename(filename), writtenToFile(false) {
        data = (T_Datatype*)malloc( sizeof(T_Datatype)*Nbins );
        clearData();
    }

    int addData( T_Datatype dataentry ) {
        // 0 <= dataentry < (max-min)/Nbins will be in 1st bin
        uint16_t ibin = floor( (dataentry - min) / ( (max-min)/Nbins ) );
        if (ibin >= Nbins)
            return 0;
        else {
            data[ibin]++;
            return ibin+1;
        }
    }

    void clearData( void ) {
        for (uint16_t i = 0; i < Nbins; i++)
            data[i] = 0;
    }

    void writeToFile( void ) {
        FILE * bindata;
        if (!writtenToFile) {
            bindata = fopen (filename,"w");
            assert(bindata != NULL);

            // comments and bin positions
            fprintf(bindata,"# 1st line is bin position. E.g. Bin Position 0.5,1.5,2.5,... mans that the first bin contains counts from [0.5,1.5)\n");
            for (uint16_t i = 0; i < Nbins; i++)
                fprintf(bindata, "%e\t", min + (i+0.5) * (max-min)/Nbins );
            fprintf(bindata, "\n");

            writtenToFile = true;
        }
        else {
            bindata = fopen (filename,"a");
            assert(bindata != NULL);
        }

        for (uint16_t i = 0; i < Nbins; i++)
            fprintf(bindata, "%e\t", data[i] );
        fprintf(bindata, "\n");

        fclose(bindata);
    }

    ~DataBinPlot( void ) {
        free( data );
    }
};

Vec getRelativeDirections(const uint32_t ex) {
        Vec tmp(0,0,0);

        switch (ex % 3) {
        case 1:
            tmp.x = 1;
            break;
        case 2:
            tmp.x = -1;
            break;
        }

        switch (ex / 3 % 3) {
        case 1: /*BOTTOM*/
            tmp.y = 1;
            break;
        case 2: /*TOP*/
            tmp.y = -1;
            break;
        }

        switch (ex / 3 / 3)  {
        case 1: /*BACK*/
            tmp.z = 1;
            break;
        case 2: /*FRONT*/
            tmp.z = -1;
            break;
        }

        return tmp;
    }

double CalcTotalKineticEnergyIons( const unsigned int particleCount, Particle particle[] ) {
    double T = 0;
    for (unsigned int i = 0; i < particleCount; i++)
        if (particle[i].q > 0)
            T += particle[i].p.norm2() / (2.0 * particle[i].m);
    return T;
}

double CalcTotalKineticEnergyElectrons( const unsigned int particleCount, Particle particle[] ) {
    double T = 0;
    for (unsigned int i = 0; i < particleCount; i++)
        if (particle[i].q < 0)
            T += particle[i].p.norm2() / (2.0 * particle[i].m);
    return T;
}

int main(void) {
    /* for (uint16_t i = 0; i <= 27; i++) {
        Vec tmp = getRelativeDirections( i );
        printf( "direction:%u => (%i,%i,%i)\n", i, int(tmp.x), int(tmp.y), int(tmp.z) );
    }
    return 0; */

    tout << "MUE0                 : " << MUE0                       << "\n";
    tout << "EPS0                 : " << EPS0                       << "\n";
    tout << "SPEED_OF_LIGHT       : " << SPEED_OF_LIGHT             << "\n";
    tout << "CELL_SIZE_SI         : " << CELL_SIZE_SI               << "\n";
    tout << "CELL_SIZE            : " << CELL_SIZE_SI / UNIT_LENGTH << "\n";
    tout << "NUMBER_OF_CELLS_X    : " << NUMBER_OF_CELLS_X          << "\n";
    tout << "NUMBER_OF_CELLS_Y    : " << NUMBER_OF_CELLS_Y          << "\n";
    tout << "NUMBER_OF_CELLS_Z    : " << NUMBER_OF_CELLS_Z          << "\n";
    tout << "DELTA_T_SI           : " << DELTA_T_SI                 << "\n";
    tout << "UNIT_ENERGY          : " << UNIT_ENERGY                << "\n";
    tout << "UNIT_MOMENTUM        : " << UNIT_MOMENTUM              << "\n";
    tout << "UNIT_ANGULAR_MOMENTUM: " << UNIT_ANGULAR_MOMENTUM      << "\n";
    tout << "ELECTRON_MASS        : " << ELECTRON_MASS              << "\n";
    tout << "\n";

    DataBinPlot<double> binCellEnergies( .08, .24, 100, "CellEnergies.dat" ); // 50 bins ranging from 3 keV to 6 keV per Cell

    if (RANDOM_SEED == 0)
        srand(time(NULL));
    else
        srand(RANDOM_SEED);

    // Energy distribution for randomly seeded cells

    uint32_t energiesOutsideRange = 0;
    uint16_t shape = 0;
    const uint32_t numberOfRuns = 20*100;
    for (uint32_t run = 0; run < 3*numberOfRuns; run++) {
        Particle electrons[NUMBER_OF_PARTICLES];
        for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
            electrons[i].m   = ELECTRON_MASS;
            if (i >= NUMBER_OF_PARTICLES / 2.)
                electrons[i].q = ELECTRON_CHARGE;
            else
                electrons[i].q =-ELECTRON_CHARGE;
            electrons[i].r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
            electrons[i].r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
            electrons[i].r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
            electrons[i].p.x = 0;
            electrons[i].p.y = 0;
            electrons[i].p.z = 0;
        }

        double E = CalcTotalPotentialEnergy( NUMBER_OF_PARTICLES, electrons, shape ) / 
                   ( NUMBER_OF_PARTICLES_PER_CELL*NUMBER_OF_CELLS_X*NUMBER_OF_CELLS_Y*NUMBER_OF_CELLS_Z );
        //tout << "E:" << E*UNIT_ENERGY*UNITCONV_Joule_to_keV << "\n";

        bool inRange = binCellEnergies.addData( E*UNIT_ENERGY*UNITCONV_Joule_to_keV );
        if (!inRange) energiesOutsideRange++;

        // File will contain one energy pos and three bincounts in the order: point, CIC, Ring
        if ((run+1) % numberOfRuns == 0) {
            if (shape == 1)
                shape = 99;
            else
                shape = 1;
            tout << "Number of Runs: " << numberOfRuns << ", Energies not in specified Range: " << energiesOutsideRange << "\n";
            energiesOutsideRange = 0;
            binCellEnergies.writeToFile();
            binCellEnergies.clearData();
        }
    }

    // Find the lowest Energy configuration: adjust Verlet with p *= 0.9999

    // Open log file
    FILE * simdata = NULL;
    simdata = fopen ("simdata.dat","w");
    fprintf( simdata, "# Horizontal: particles, Vertical: x,y,x,y,...\n" );
    assert(simdata != NULL);

    FILE * stats = NULL;
    stats = fopen ("stats.dat","w");
    fprintf( stats, "# t\tE/keV\tV/keV\tL/(m*kg*m/s)\tP/(kg*m/s)\tT/keV\n" );
    assert(stats != NULL);

    Particle electrons[NUMBER_OF_PARTICLES];
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
        electrons[i].m   = ELECTRON_MASS;
        if (i >= NUMBER_OF_PARTICLES / 2.)
            electrons[i].q = ELECTRON_CHARGE;
        else
            electrons[i].q =-ELECTRON_CHARGE;
        electrons[i].r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
        electrons[i].r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
        electrons[i].r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
        electrons[i].p.x = 0;
        electrons[i].p.y = 0;
        electrons[i].p.z = 0;
    }

    tout << "Stats fprintf interval: " << PRINTF_SIMDATA_INTERVAL << " meaning every " << PRINTF_SIMDATA_INTERVAL * DELTA_T_SI << " s\n";

    // Initialize Verlet Integration (Leapfrog): shift momentum half a time step
    Verlet( NUMBER_OF_PARTICLES, electrons, DELTA_T, true);
    for (uint32_t currentStep = 0; currentStep < NUMBER_OF_STEPS; currentStep++) {
        Verlet( NUMBER_OF_PARTICLES, electrons, DELTA_T );

        if ((currentStep % PRINTF_INTERVAL == 0) or (currentStep % PRINT_INTERVAL == 0)) {
            const double Te = CalcTotalKineticEnergyElectrons  ( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double Ti = CalcTotalKineticEnergyIons( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double T  = CalcTotalKineticEnergy    ( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double V  = CalcTotalPotentialEnergy  ( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double L  = CalcTotalAngularMomentum  ( NUMBER_OF_PARTICLES, electrons ) * UNIT_ANGULAR_MOMENTUM;
            const double P  = CalcTotalMomentum         ( NUMBER_OF_PARTICLES, electrons ) * UNIT_MOMENTUM;
            const double E  = T + V;
            if (currentStep % PRINTF_INTERVAL == 0) 
                fprintf( stats, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", currentStep * DELTA_T_SI, E, V, L, P, T, Te, Ti);
            if (currentStep % PRINTF_SIMDATA_INTERVAL == 0) {
                for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
                    fprintf( simdata, "%e\t", electrons[i].r.x );
                fprintf( simdata, "\n" );
                for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
                    fprintf( simdata, "%e\t", electrons[i].r.y );
                fprintf( simdata, "\n" );
            }
            if (currentStep % PRINT_INTERVAL == 0)
                tout << "time: " << currentStep * DELTA_T_SI << " (" << 100*(double)currentStep/(double)NUMBER_OF_STEPS << "%), E: " << E << ", V: " << V << ", L: " << L << ", P:" << P << "\n";
        }
    }

    fclose(stats);
    fclose(simdata);

    return 0;
}
