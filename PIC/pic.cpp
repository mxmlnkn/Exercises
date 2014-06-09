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
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <omp.h>
#include <iostream>

#define DEBUG 1

using namespace std;

#include "Vector.cpp"
#include "Parameters.cpp"

const uint16_t STEPS_TO_REMEMBER = 2;
struct Particle {
    // previous steps needed for example for relativistic formulas because of time retardation
    Vec r,p;
    Vec r_prev[STEPS_TO_REMEMBER],p_prev[STEPS_TO_REMEMBER];
    Vec dr,dp;
    double m,q;
};

inline Vec ForceActingOnParticle1( const Particle & particle1, const Particle & particle2, uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
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
    
    // Loop about surrounding cells where the particle is also because of periodic boundary conditions 
    // Try to find analytic formula for this ??? Like in solid state physics !
    assert( SIMDIM == 2 );
    for (int iPeriodY = -NUMBER_OF_NEIGHBORS; iPeriodY <= NUMBER_OF_NEIGHBORS; iPeriodY++)
    for (int iPeriodX = -NUMBER_OF_NEIGHBORS; iPeriodX <= NUMBER_OF_NEIGHBORS; iPeriodX++) {
        
        const Vec r12  = particle1.r - ( particle2.r + Vec( iPeriodX * CELL_SIZE_X, iPeriodY * CELL_SIZE_Y, 0) );
        const double r = r12.norm();
        Vec F12 = alpha * r12/r;
        
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
        
        F += F12;
    }
    
    return F;
}

inline double Energy( const Particle & particle1, const Particle & particle2, uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    const double cellWidth = CELL_SIZE_X;
    const Vec r12  = particle1.r - particle2.r;
    const double r = r12.norm();
    // for two positive charges the potential energy is supposed to increase for smaller distances, meaning the particles slow down
    double V = particle1.q * particle2.q / ( 4.*M_PI*EPS0 );

    if (colshape == 99) {
        if (r < cellWidth)
            // r=cellWidth => V = 1/cellWidth => continuous
            V *= 2./cellWidth - r/(cellWidth*cellWidth);
    } else if (colshape == 0)
        V = V / r;
    else if( colshape == 1) {
        const double radii_ratio = r / ( cellWidth/2.0 );
        if (r < cellWidth)
            // r=cellWidth => radii_ratio=2 => V = 1/cellWidth (correct :) )
            V *= (192. - 80. * pow(radii_ratio,2) + 30 * pow(radii_ratio,3) - pow(radii_ratio,5) ) / (160.*cellWidth/2.0);
    }
    return V;
}

/* checks if particles has left the simulation space. Should be called every  *
 * time the position has been changed. Therefore would be nice to implement   *
 * this in all positions change routines, but then Vector.hpp wouldn't be     *
 * general anymore                                                            *
 * If some particles has left, then but it back to the wall in that dimension *
 * and clear the momentum of that direction like a infinitely small elastic   *
 * coefficient                                                                */
bool CheckBoundaryConditions( Particle & particle ) {
    bool outOfBounds = false;
    // Cycle through all (max. 3) dimensions
    for (uint16_t dim = 0; dim < SIMDIM; dim++) {
        if ( (particle.r[dim] < 0) or (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) )
            outOfBounds = true;

        #if BOUNDARY_CONDITION == 0
        /* Periodic Boundary Conditions */
        if (particle.r[dim] < 0) {
            particle.r[dim] += NUMBER_OF_CELLS[dim] * CELL_SIZE[dim];
        } else if (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
            particle.r[dim] -= NUMBER_OF_CELLS[dim] * CELL_SIZE[dim]; //not safe for two cells at one push ! => courant criterium !
        }
        
        #elif BOUNDARY_CONDITION == 1
        /* Reflecting Boundary Condition */
        // This will make an error, if a particle is out of bounds in more than one dimension !
        if (particle.r[dim] < 0) {
            double t1 = ( 0 - particle.r[dim] ) / ( particle.p[dim] / particle.m );
            particle.p[dim] *= -1;
            particle.r[dim] = 0 + (DELTA_T - t1) * ( particle.p[dim] / particle.m );
        } else if (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
            double t1 = ( particle.r[dim] - NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) / ( particle.p[dim] / particle.m );
            particle.p[dim] *= -1;
            particle.r[dim] = NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] + (DELTA_T - t1) * ( particle.p[dim] / particle.m );
        }
        
        #elif BOUNDARY_CONDITION == 2
        /* Adhering Boundary Condition */
        if (particle.r[dim] < 0) {
            particle.r[dim] = 0;
            particle.p[dim] = 0;
        } else if (particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
            particle.r[dim] = 0;
            particle.p[dim] = 0;
        }
        
        #endif
    }
    return outOfBounds;
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
            if (i==j) continue;
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

double CalcTotalEnergy( const unsigned int particleCount, struct Particle particle[], uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    double E = 0;
    for (unsigned int i = 0; i < particleCount; i++) {
        double p = particle[i].p.norm();
        double T = p*p / (2.0 * particle[i].m);
        double V = 0;
        for (unsigned int j = 0; j < i; j++)
            V += Energy( particle[i], particle[j], colshape );
        E += T + V;
    }
    return E;
}

double CalcTotalPotentialEnergy( const unsigned int particleCount, struct Particle particle[], uint16_t colshape = DEFAULT_PARTICLE_SHAPE ) {
    double E = 0;
    for (unsigned int i = 0; i < particleCount; i++) {
        double V = 0;
        for (unsigned int j = 0; j < i; j++)
            V += Energy( particle[i], particle[j], colshape );
        E += V;
    }
    return E;
}

double CalcTotalAngularMomentum( const unsigned int particleCount, struct Particle particle[] ) {
    Vec L;
    for (unsigned int i = 0; i < particleCount; i++) {
        L.x += particle[i].r.y * particle[i].p.z - particle[i].r.z * particle[i].p.y;
        L.y += particle[i].r.z * particle[i].p.x - particle[i].r.x * particle[i].p.z;
        L.z += particle[i].r.x * particle[i].p.y - particle[i].r.y * particle[i].p.x;
    }
    return L.norm();
}

double CalcTotalMomentum( const unsigned int particleCount, struct Particle particle[] ) {
    Vec P;
    for (unsigned int i = 0; i < particleCount; i++) {
        P.x += particle[i].p.x;
        P.y += particle[i].p.y;
        P.z += particle[i].p.z;
    }
    return P.norm();
}


struct return_data {
    double theta_max, theta_max_calc, rmin, rmin_calc;
};


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

int main(void) {
    printf("MUE0                 : %e\n",MUE0);
    printf("EPS0                 : %e\n",EPS0);
    printf("SPEED_OF_LIGHT       : %e\n",SPEED_OF_LIGHT);
    printf("CELL_SIZE_SI         : %e\n",CELL_SIZE_SI);
    printf("CELL_SIZE            : %e\n",CELL_SIZE_SI / UNIT_LENGTH);
    printf("NUMBER_OF_CELLS_X    : %d\n",NUMBER_OF_CELLS_X);
    printf("DELTA_T_SI           : %e\n",DELTA_T_SI);
    printf("UNIT_ENERGY          : %e\n",UNIT_ENERGY);
    printf("UNIT_MOMENTUM        : %e\n",UNIT_MOMENTUM);
    printf("UNIT_ANGULAR_MOMENTUM: %e\n",UNIT_ANGULAR_MOMENTUM);
    printf("ELECTRON_MASS        : %e\n",ELECTRON_MASS);
    printf("\n");

    DataBinPlot<double> binCellEnergies( 2., 6., 100, "CellEnergies.dat" ); // 50 bins ranging from 3 keV to 6 keV

    srand(time(NULL));

    // Energy distribution for randomly seeded cells

    uint32_t energiesOutsideRange = 0;
    uint16_t shape = 0;
    const uint32_t numberOfRuns = 20*100;
    for (uint32_t run = 0; run < 3*numberOfRuns; run++) {
        Particle electrons[NUMBER_OF_PARTICLES];
        for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
            electrons[i].m   = ELECTRON_MASS;
            electrons[i].q   = ELECTRON_CHARGE;
            electrons[i].r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
            electrons[i].r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
            electrons[i].r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
            electrons[i].p.x = 0;
            electrons[i].p.y = 0;
            electrons[i].p.z = 0;
        }

        double E = CalcTotalEnergy( NUMBER_OF_PARTICLES, electrons, shape );
        //printf("%e\n", E*UNIT_ENERGY*UNITCONV_Joule_to_keV);

        //printf( "%e\t\n", E*UNIT_ENERGY*UNITCONV_Joule_to_keV );
        bool inRange = binCellEnergies.addData( E*UNIT_ENERGY*UNITCONV_Joule_to_keV );
        if (!inRange) energiesOutsideRange++;

        // File will contain one energy pos and three bincounts in the order: point, CIC, Ring
        if ((run+1) % numberOfRuns == 0) {
            if (shape == 1)
                shape = 99;
            else
                shape = 1;
            printf( "Number of Runs: %d, Energies not in specified Range: %d\n", numberOfRuns, energiesOutsideRange );
            energiesOutsideRange = 0;
            binCellEnergies.writeToFile();
            binCellEnergies.clearData();
        }
    }

    // Find the lowest Energy configuration

    // Open log file
    FILE * simdata = NULL;
    simdata = fopen ("simdata.dat","w");
    fprintf( simdata, "# Horizontal: particles, Vertical: x,y,x,y,...\n" );
    assert(simdata != NULL);

    FILE * stats = NULL;
    stats = fopen ("stats.dat","w");
    fprintf( stats, "# t\tE/keV\tV/keV\tL/(m*kg*m/s)\tP/(kg*m/s)\n" );
    assert(stats != NULL);

    Particle electrons[NUMBER_OF_PARTICLES];
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
        electrons[i].m   = ELECTRON_MASS;
        electrons[i].q   = ELECTRON_CHARGE;
        electrons[i].r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
        electrons[i].r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
        electrons[i].r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
        electrons[i].p.x = 0;
        electrons[i].p.y = 0;
        electrons[i].p.z = 0;
    }

    const uint32_t printInterval   = 2000;
    const uint32_t fprintInterval  = 200;

    // Initialize Verlet Integration (Leapfrog): shift momentum half a time step
    Verlet( NUMBER_OF_PARTICLES, electrons, DELTA_T, true);
    for (uint32_t currentStep = 0; currentStep < NUMBER_OF_STEPS; currentStep++) {
        Verlet( NUMBER_OF_PARTICLES, electrons, DELTA_T );

        if ((currentStep % fprintInterval == 0) or (currentStep % printInterval == 0)) {
            const double E = CalcTotalEnergy         ( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double V = CalcTotalPotentialEnergy( NUMBER_OF_PARTICLES, electrons ) * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double L = CalcTotalAngularMomentum( NUMBER_OF_PARTICLES, electrons ) * UNIT_ANGULAR_MOMENTUM;
            const double P = CalcTotalMomentum       ( NUMBER_OF_PARTICLES, electrons ) * UNIT_MOMENTUM;
            if (currentStep % fprintInterval == 0) {
                fprintf( stats, "%e\t%e\t%e\t%e\t%e\n", currentStep * DELTA_T_SI, E, V, L, P);
                for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
                    fprintf( simdata, "%e\t", electrons[i].r.x );
                fprintf( simdata, "\n" );
                for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
                    fprintf( simdata, "%e\t", electrons[i].r.y );
                fprintf( simdata, "\n" );
            }
            if (currentStep % printInterval == 0)
                printf( "time: %e (%3.0f%%), E: %e, V: %e, L: %e, P:%e\n", currentStep * DELTA_T_SI, 100*(double)currentStep/(double)NUMBER_OF_STEPS, E, V, L, P);
        }
    }

    fclose(stats);
    fclose(simdata);

    return 0;
}
