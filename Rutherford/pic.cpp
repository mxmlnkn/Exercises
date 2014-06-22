// del pic.exe && g++ pic.cpp -o pic.exe -Wall -Wextra -Wimplicit -pedantic-errors -pedantic -Wcomment -Wconversion -Wuninitialized -std=c++0x -O3 && start pic.exe
/******************************************************************************
 * ToDo:                                                                      *
 *   - make work for center of mass (i.e. for similar particle masses)        *
 *   - substitute verlet with analytical integrated trajectory (but only in   *
 *     one time step!)                                                        *
 ******************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <list>
#include <random>

#define M_PI 3.14159265358979323846 // C++11 HASS GODDAMN IT !!!

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
Vec simBoxp(0,0,0);

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
                assert(false);
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

    const double cloudRadius = PARTICLE_RADIUS;
    const double alpha       = particle1.q * particle2.q / ( 4.*M_PI*EPS0 );
    Vec F(0,0,0);

    // Loop over periodically continued simulation box i.e. Supercells, because this simulation only supports one supercell at a time with not so clearly separated cells
    const uint16_t supercells_to_consider_X = int( ceil( CELL_SIZE_X * CONSIDERATION_RATIO / ( NUMBER_OF_CELLS_X * CELL_SIZE_X ) ) + 0.1 );
    const uint16_t supercells_to_consider_Y = int( ceil( CELL_SIZE_Y * CONSIDERATION_RATIO / ( NUMBER_OF_CELLS_Y * CELL_SIZE_Y ) ) + 0.1 );
    const uint16_t supercells_to_consider_Z = int( ceil( CELL_SIZE_Z * CONSIDERATION_RATIO / ( NUMBER_OF_CELLS_Z * CELL_SIZE_Z ) ) + 0.1 );
    uint16_t supercells_to_consider         = max( supercells_to_consider_X, max( supercells_to_consider_Y, supercells_to_consider_Z ) );
    if (PERIODIC_FORCE == false)
        supercells_to_consider = 0;

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

        if( r < CELL_SIZE_MIN * CONSIDERATION_RATIO or !PERIODIC_FORCE ) {
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
    const double cloudRadius = PARTICLE_RADIUS;
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

class SimulationBox {
public:
    Vec p;
    SimulationBox();
    ~SimulationBox();
    double E0, Te0, Ti0;
    Vec L0;
    uint64_t numberOfAppliedBoundaryConditions;
    Vec P0;
    list<Particle> electrons[ NUMBER_OF_CELLS_X ][ NUMBER_OF_CELLS_Y ][ NUMBER_OF_CELLS_Z ];
    list<Particle>      ions[ NUMBER_OF_CELLS_X ][ NUMBER_OF_CELLS_Y ][ NUMBER_OF_CELLS_Z ];

    FILE * simdata_eons;
    FILE * simdata_ions;
    void DumpData( void );

    double CalcTotalPotentialEnergy( void );
    double CalcTotalKineticEnergyElectrons( void );
    double CalcTotalKineticEnergyIons( void );
    double CalcTotalKineticEnergy( void );
    double CalcTotalEnergy( void );
    Vec CalcTotalMomentum( void );
    Vec MaxMomentum( void );
    Vec CalcTotalAngularMomentum( void );

    void PushLocation( const double dt );
    void ClearDp( void );
    void CalcDp ( const double dt );
    void ApplyDp( void );
    void Verlet ( const double dt, bool initialize = false );
    bool CheckBoundaryConditions( Particle & particle );
};

SimulationBox::SimulationBox() : p(Vec(0,0,0)), numberOfAppliedBoundaryConditions(0) {

    // Open log file
    this->simdata_eons = fopen ("Electrons.dat","w");
    fprintf( simdata_eons, "# Horizontal: particles, Vertical: x,y,x,y,...\n" );
    assert(simdata_eons != NULL);
    this->simdata_ions = fopen ("Ions.dat","w");
    fprintf( simdata_ions, "# Horizontal: particles, Vertical: x,y,x,y,...\n" );
    assert(simdata_ions != NULL);
/*
    FILE * stats = NULL;
    stats = fopen ("stats.dat","w");
    fprintf( stats, "# t\tE/keV\tV/keV\tL/(m*kg*m/s)\tP/(kg*m/s)\tT/keV\tTe/keV\tTi/keV(E-E0)/keV\n" );
    assert(stats != NULL);*/
}

SimulationBox::~SimulationBox(void)     {
    fclose(simdata_eons);
    fclose(simdata_ions);
}

void SimulationBox::DumpData( void ) {
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            fprintf( simdata_eons, "%e\t", i->r.x + ix );
    }
    fprintf( simdata_eons, "\n" );

    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            fprintf( simdata_eons, "%e\t", i->r.y + iy );
    }
    fprintf( simdata_eons, "\n" );

    // Ions
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & ions = this->ions[ix][iy][iz];
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            fprintf( simdata_ions, "%e\t", i->r.x + ix );
    }
    fprintf( simdata_ions, "\n" );

    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & ions = this->ions[ix][iy][iz];
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            fprintf( simdata_ions, "%e\t", i->r.y + iy );
    }
    fprintf( simdata_ions, "\n" );
}

SimulationBox simBox;

void SimulationBox::ClearDp( void ) {
    #pragma omp parallel for
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            i->dp = Vec(0,0,0);
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            i->dp = Vec(0,0,0);
    }
}

/* Loop Template
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        // electron-electron and electron-ion pairing
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            // electron-electron pairing
            for (list<Particle>::iterator j = eons.begin(); j != i; j++) // j != i is only equivalent to j < i, if both loops are incremented in the same manner, which should be true if the list don't change (do CheckBoundary somewhere else). It it's not true... fuck STL...
                V += 0.5 * Energy( *i, *j );
            // electron-ion pairing
            for (list<Particle>::iterator j = ions.begin(); j != ions.end(); j++) // Because different list, go through completely
                V += Energy( *i, *j );
        }
        // ion-ion pairing
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++) {
            for (list<Particle>::iterator j = ions.begin(); j != j; j++)
                V += 0.5 * Energy( *i, *j );
        }
    }
*/


void SimulationBox::CalcDp( const double dt ) {
    this->ClearDp();
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        // electron-electron and electron-ion pairing
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            // electron-electron pairing
            for (list<Particle>::iterator j = eons.begin(); j != i; j++) { // exclude at least self-pairing. j != i is only equivalent to j < i, if both loops are incremented in the same manner, which should be true if the list don't change (do CheckBoundary somewhere else). It it's not true... fuck STL...
                const Vec dp12 = ForceActingOnParticle1( *i, *j ) * dt;
                i->dp += dp12;
                j->dp -= dp12;    // only possible if (i,j) and (j,i) aren't both calculated => j <= i condition, and j != i, because that is self collision
            }
            // electron-ion pairing
            for (list<Particle>::iterator j = ions.begin(); j != ions.end(); j++) { // Because different list, go through completely
                const Vec dp12 = ForceActingOnParticle1( *i, *j ) * dt;
                i->dp += dp12;
                j->dp -= dp12;
            }
        }
        // ion-ion pairing
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            for (list<Particle>::iterator j = ions.begin(); j != i; j++) {
                const Vec dp12 = ForceActingOnParticle1( *i, *j ) * dt;
                i->dp += dp12;
                j->dp -= dp12;
            }
    }
}

/* checks if particles has left the simulation space. Should be called every  *
 * time the position has been changed. Therefore would be nice to implement   *
 * this in all positions change routines, but then Vector.hpp wouldn't be     *
 * general anymore                                                            *
 * If some particles has left, then put it back to the wall in that dimension *
 * and clear the momentum of that direction like a infinitely small elastic   *
 * coefficient                                                                */
bool SimulationBox::CheckBoundaryConditions( Particle & particle ) {
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

        #define AWESOME_REFLEXION 1
        else if (BOUNDARY_CONDITION == 1) {
            /* Reflecting Boundary Condition */
            if (particle.r[dim] < 0) {
                #if AWESOME_REFLEXION == 1
                    //const double v  = particle.p[dim] / particle.m;
                    //const double t1 = ( particle.r[dim] - 0 ) / v;
                    //particle.r[dim] = 0 + (DELTA_T - t1) * (-v);   // seems to result in the best energy conservation ever seen ... That's fucking awesome Oo !!!. This basically symmetrizes the reflection. Meaning the length of the incoming 'ray' is the same as the outgoing... But this means, that in extreme cases where x = -epsilon, t1->0 => x = Delta_t * |v|. Problem is, that this should be the same v used to calculate x = -epsilon = x(i-1) + v*dt, therefore moving the particle v*2*dt in total ... I don't really get why this works, but it works... could by a nice discovery... the symmetry is awesome, so ... also this just seems compatible to the Verlet method
                    particle.r[dim] = particle.r[dim] - DELTA_T * particle.p[dim] / particle.m; // so instead of correctly being reflected it's like the particle never moved -> the step is being inversed :S... This can actually explain the small errors. Those occur when a = (a -c) + c !=a. Although this method will conserve energy
                #else
                    particle.r[dim] = 0 + ( 0 - particle.r[dim] ); // seems to be the correct way, but is not energy conserving, the total time the particle needs for a certain distance (e.g. let it reflect 200x inside the box) will have an error ... Compare that error with the error that would result from higher kinetic energy !!!
                #endif
                particle.p[dim] *= -1;
                #pragma omp atomic
                simBox.p[dim] -= 2.*particle.p[dim];
            } else if ( particle.r[dim] > NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) {
                #if AWESOME_REFLEXION == 1
                    //const double v  = particle.p[dim] / particle.m;
                    //const double t1 = ( particle.r[dim] - NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] ) / v;
                    //particle.r[dim] = NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] - (DELTA_T - t1) * v;
                    particle.r[dim] = particle.r[dim] - DELTA_T * particle.p[dim] / particle.m; // inverse last Verlet-step/push
                #else
                    particle.r[dim] = NUMBER_OF_CELLS[dim] * CELL_SIZE[dim] - (particle.r[dim] - NUMBER_OF_CELLS[dim] * CELL_SIZE[dim]);
                #endif
                particle.p[dim] *= -1;
                #pragma omp atomic
                simBox.p[dim] -= 2.*particle.p[dim];
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
    this->numberOfAppliedBoundaryConditions += outOfBounds;
    return outOfBounds > 0;
}

Vec Velocity( Vec p, double m ) {
    // for the non-relativistic limit the result will be p/m. For ultrarelativistic case result will be c*e_p
    return p / m;
    //return p / sqrt( m*m + p.norm2() / ( SPEED_OF_LIGHT * SPEED_OF_LIGHT ) );
}

void SimulationBox::ApplyDp( void ) {
    #pragma omp parallel for
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            i->p += i->dp;
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            i->p += i->dp;
    }
}

void SimulationBox::PushLocation( const double dt ) {
    #pragma omp parallel for
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            i->r += Velocity( i->p, i->m ) * dt;
            this->CheckBoundaryConditions( *i );
        }
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++) {
            i->r += Velocity( i->p, i->m ) * dt;
            this->CheckBoundaryConditions( *i );
        }
    }
}


void SimulationBox::Verlet( const double dt, bool initialize ) {
    /* Algorithm:                                   *
     *     p += ForceActingOnThisParticle * dt/2.0; *
     *     r += p / m * dt/2.0;                     *
     *     r += p / m * dt/2.0;                     *
     *     p += ForceActingOnThisParticle * dt/2.0; */

    this->CalcDp( 0.5 * dt );
    this->ApplyDp();

    // First Verlet step only moves particle momentum half a step and doesn't need the other steps
    if (initialize)
        return;

    this->PushLocation( dt ); // could be combined in applyDp too save time when iterating over the lists :S

    /* Reimplement this:
        // r_prev[highest] is the oldest position => shift r_prev[i] to r_prev[i+1]
        for (uint16_t j = STEPS_TO_REMEMBER-1; j > 0; j--) {
            i->r_prev[j] = i->r_prev[j-1];
            i->p_prev[j] = i->p_prev[j-1];
        }
        i->r_prev[0] = i->r;
        i->p_prev[0] = i->p;
    */

    this->CalcDp( 0.5 * dt );
    this->ApplyDp();

    // Apply simBox Momentum -> Galilei Transformation (also true for LT?)
    /* This wouldn't behave nicely for one particle for example
    #pragma omp parallel for
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++)
        particle[i].dp = simBox.p / NUMBER_OF_PARTICLES;
    simBox.p = Vec(0,0,0); */

    return;
}

double SimulationBox::CalcTotalKineticEnergy( void ) {
    return this->CalcTotalKineticEnergyIons() + this->CalcTotalKineticEnergyElectrons();
}

double SimulationBox::CalcTotalEnergy( void ) {
    return this->CalcTotalKineticEnergy() + this->CalcTotalPotentialEnergy();
}

double SimulationBox::CalcTotalPotentialEnergy( void ) {
    double E = 0;
    // for every cell calculate potential of particles in the same cell. If we want the full potential over all cells the loop will have to be more complex !!!
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        // electron-electron and electron-ion potential
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            double V = 0;
            for (list<Particle>::iterator j = eons.begin(); j != i; j++)  // i < j, to exclude double calculations between particle 1-2 and 2-1 would have better performance, but list::iterator is not capable of operator<
                if ( i != j ) // exclude at least self-collisions
                    V += Energy( *i, *j );
            for (list<Particle>::iterator j = ions.begin(); j != ions.end(); j++)  // because the inner list ist not the same as the outer we don't use i < j, but the fool loop
                V += Energy( *i, *j );
            E += V;
        }
        // ion-ion potential
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++) {
            double V = 0;
            for (list<Particle>::iterator j = ions.begin(); j != i; j++)
                if ( i != j ) // exclude self-collisions
                    V += Energy( *i, *j );
            E += V;
        }
    }
    return E;
}

Vec SimulationBox::CalcTotalAngularMomentum( void ) {
    Vec L(0,0,0);
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        // electrons
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            L.x += i->r.y * i->p.z - i->r.z * i->p.y;
            L.y += i->r.z * i->p.x - i->r.x * i->p.z;
            L.z += i->r.x * i->p.y - i->r.y * i->p.x;
        }
        // ions
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++) {
            L.x += i->r.y * i->p.z - i->r.z * i->p.y;
            L.y += i->r.z * i->p.x - i->r.x * i->p.z;
            L.z += i->r.x * i->p.y - i->r.y * i->p.x;
        }
    }
    return L;
}

Vec SimulationBox::CalcTotalMomentum( void ) {
    Vec P(0,0,0);
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        list<Particle> & ions = this->ions     [ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            P += i->p;
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            P += i->p;
    }
    return P;
}

Vec SimulationBox::MaxMomentum( void ) {
    Vec P(0,0,0);
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & ions = this->ions[ix][iy][iz];
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++) {
            P.x = max( P.x, abs( i->p.x ) );
            P.y = max( P.y, abs( i->p.y ) );
            P.z = max( P.z, abs( i->p.z ) );
        }
        list<Particle> & eons = this->electrons[ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++) {
            P.x = max( P.x, abs( i->p.x ) );
            P.y = max( P.y, abs( i->p.y ) );
            P.z = max( P.z, abs( i->p.z ) );
        }
    }
    return P;
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
        uint16_t ibin = uint16_t( (dataentry - min) / ( (max-min)/Nbins ) );
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

double SimulationBox::CalcTotalKineticEnergyIons( void ) {
    double T = 0;
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & ions = this->ions[ix][iy][iz];
        for (list<Particle>::iterator i = ions.begin(); i != ions.end(); i++)
            T += i->p.norm2() / (2.0 * i->m);
    }
    return T;
}

double SimulationBox::CalcTotalKineticEnergyElectrons( void ) {
    double T = 0;
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
        list<Particle> & eons = this->electrons[ix][iy][iz];
        for (list<Particle>::iterator i = eons.begin(); i != eons.end(); i++)
            T += i->p.norm2() / (2.0 * i->m);
    }
    return T;
}

int main(void) {
    assert( ((NUMBER_OF_PARTICLES_PER_CELL % 2) == 0) or SPECIES != 3 );

    tout << "MUE0                 : " << MUE0                       << "\n";
    tout << "EPS0                 : " << EPS0                       << "\n";
    tout << "SPEED_OF_LIGHT       : " << SPEED_OF_LIGHT             << "\n";
    tout << "CELL_SIZE_SI         : " << CELL_SIZE_SI               << "\n";
    tout << "CELL_SIZE            : " << CELL_SIZE_SI / UNIT_LENGTH << "\n";
    tout << "PARTICLE_RADIUS      : " << PARTICLE_RADIUS            << "\n";
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

    // Open log file
    FILE * stats = NULL;
    stats = fopen ("stats.dat","w");
    fprintf( stats, "# t\tE/keV\tV/keV\tL/(m*kg*m/s)\tP/(kg*m/s)\tT/keV\tTe/keV\tTi/keV(E-E0)/keV\n" );
    assert(stats != NULL);

    Particle particle;
    particle.m    = ELECTRON_MASS;
    particle.q    = ELECTRON_CHARGE;

    // coulomb force factor: |F| = alpha/r*r, V = alpha/r
    double alpha  = ELECTRON_CHARGE*ELECTRON_CHARGE / ( 4*M_PI*EPS0 );
    // collision parameter (parallel distance)
    double s      = 0.1*CELL_SIZE_X;//250e-15;
    // initial momentum in infinity in x direction, here: T = 10eV = p*p/2m => p = sqrt(2mT)
    double p      = sqrt( 2* ELECTRON_MASS* 0.06 * UNITCONV_keV_to_Joule / UNIT_ENERGY );
    // E(-inf) = E_kin(-inf)
    double E      = p*p / ( 2*ELECTRON_MASS );
    // angular momentum for t=-inf. L = r.x * p.y - r.y * p.x = 0 - p*s
    double L      =-p*s;
    // dimensionless energy. just some term which we need more than once
    double Ew     = sqrt( 2*E / ELECTRON_MASS ) * L/alpha;
    /* minimum distance given by v_radial(!) = 0 => E = E_pot(t) +    *
     * E_kin_radial(t) + E_kin_angular(t) = alpha/r_min + 0 +         *
     * + L*L/2*m*r_min^2 = E_pot(-inf) + E_kin(-inf) = 0 + p*p/2*m    *
     *  => r_min =  alpha/2E *[ 1 +- sqrt( 1 + 4*E*L*L/2*m*alpha^2 )] *
     *           =  alpha*m/p^2 * [ 1 + sqrt( 1 + (p*L/m*alpha)^2 ) ] *
     *           =: alpha*m/p^2 * [ 1 + sqrt( 1 + Ew*Ew ) ]           *
     *  Ew = p*L/m*alpha = sqrt(2*E/m) L/alpha
     * alpha/(p*p/2m)). This formula will be the same if Ew = 0.      */
    double rmin   = alpha/(2*E) * ( 1 + sqrt( 1 + Ew*Ew ) );
    tout << "E(t=-inf):" << E * UNIT_ENERGY * UNITCONV_Joule_to_keV << ", L(t=-inf):" << L * UNIT_ANGULAR_MOMENTUM << "\n";

    // initial radius for particle at simulation start
    double r0     = 2*rmin;//rmin * pow(10,r0Expo);
    assert( r0 > rmin ); // else asin later on will fail
    /* initial momentum calculated simply from total energy, as there *
     * is no need here to distinguish between radial and angular      *
     * momentum. We did this only because we needed especially the    *
     * radial momentum to be set 0                                    *
     * momentum conservation: E(t) = alpha/r(t) + p(t)^2/2m           *
     *  => p = sqrt[ ( E - alpha/r0 - L*L/2*m*r0 ) * 2m ]             */
    double p0   = sqrt(2*ELECTRON_MASS) * sqrt( E - alpha/r0 );
    /* Needed for x,y EXPLAIN DERIVATION !!! */
    double phi0 = M_PI + asin( 1 / sqrt(Ew*Ew+1) )
                       - asin( ( Ew*L/( sqrt(2*ELECTRON_MASS*E)*r0 ) + 1 ) / sqrt(Ew*Ew+1) );

    /* r.x * 10 => deviation of theta / 10 approx. Shouldn't be because       *
     * of wrongly calculated p vs. p_inf, because p_inf is correctly          *
     * calculated from total energy. Therefore it should be a geometrical     *
     * problem, because we want to know the asymptotes of the trajectories.   *
     * Is also almost independent of Integration method for 1e5 or more steps *
     * p may be correct, but the direction of p could still be error prone !  */
    const double x =-r0*cos(M_PI - phi0);
    const double y = r0*sin(M_PI - phi0);
    /* r.x * p.y - r.y * p.x = L and p.x^2 + p.y^2 = p0^2             *
     * x*sqrt(p0^2 - px^2) - y*px = L =>                              *
     *      ( p0^2 - px^2 ) * x^2 = L^2 + px^2*y^2 + 2*L*px*y <=>     *
     *  0 = px^2 ( x^2 + y^2 ) + 2*L*px*y + L^2 - p0^2*x^2     =>     *
     * px =-L*y/r0^2 +- sqrt[(L*y/r0^2)^2 - L^2/r0^2 + p0^2*x^2/r0^2] *
     *    =-L/r0^2 * { y + sqrt[ y^2 - r0^2 + p0^2 x^2 r0^2/ L^2 ] }  *
     *    =-L/r0^2 * { y + sqrt[ -1 + p0^2 r0^2/ L^2 ] *x }           */
    particle.p.x = -L/(r0*r0) * ( y - x * sqrt( -1 + pow( r0*p0/L ,2) ) );
    particle.p.y = sqrt(p0*p0 - particle.p.x * particle.p.x);
    particle.p.z = 0;
    particle.r.x = x + 0.50 * CELL_SIZE_X;
    particle.r.y = y + 0.25 * CELL_SIZE_X;
    particle.r.z = 0;
    simBox.electrons[0][0][0].push_back( particle );
    // total scattering, so the angle between the incoming and outgoing trajectory asymptote
    double theta_inf = 2*asin( 1/sqrt( 1 + Ew*Ew ) ); 

    double initial_angle = acos( particle.p * Vec(1,0,0) / particle.p.norm() ) *180./M_PI;
    printf( "phi0:%E, p0:%E, p_inf:%E, <(p0,x):%E\n", (M_PI-phi0)/M_PI*180, p0,p,initial_angle );
    printf("     ( %e )         ( %e )\n  ", particle.r.x, particle.p.x );
    printf("r0 = ( %e )  , p0 = ( %e )\n  ", particle.r.y, particle.p.y );
    printf("     ( %e )         ( %e )\n\n", particle.r.z, particle.p.z );
    
    // the much much heavier ion which is approximately at rest
    particle.m   = ION_MASS;
    particle.q   = ION_CHARGE;
    particle.r.x = 0.50 * CELL_SIZE_X;
    particle.r.y = 0.25 * CELL_SIZE_Y;
    particle.r.z = 0;
    particle.p.x = 0;
    particle.p.y = 0;
    particle.p.z = 0;
    simBox.ions[0][0][0].push_back( particle );

    // Initial analysis and output
    simBox.Te0 = simBox.CalcTotalKineticEnergyElectrons();
    simBox.Ti0 = simBox.CalcTotalKineticEnergyIons();
    simBox.L0  = simBox.CalcTotalAngularMomentum();
    simBox.E0  = simBox.Te0 + simBox.Ti0 + simBox.CalcTotalPotentialEnergy();
    simBox.P0  = simBox.CalcTotalMomentum();
    tout << "Initial total Momentum P0.x : " << simBox.P0.x * UNIT_MOMENTUM << " kg m/s\n";
    tout << "Initial total Ang.Mom L0.z  : " << simBox.L0.x * UNIT_ANGULAR_MOMENTUM << " kg m^2/s\n";
    tout << "Initial total   Energy E0   : " << simBox.E0   * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Initial kinetic Energy Te0  : " << simBox.Te0  * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Initial kinetic Energy Ti0  : " << simBox.Ti0  * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Note that <E_kin_e> = Te0/N_Particles = " << simBox.Te0 / (NUMBER_OF_PARTICLES/2.0) * UNIT_ENERGY * UNITCONV_Joule_to_keV * 1000 << "eV = f/2 kT, where f=" << SIMDIM << "\n";
    tout << "Note that <E_kin_i> = Ti0/N_Particles = " << simBox.Ti0 / (NUMBER_OF_PARTICLES/2.0) * UNIT_ENERGY * UNITCONV_Joule_to_keV * 1000 << "eV = f/2 kT, where f=" << SIMDIM << "\n";

    tout << "Stats fprintf interval: " << PRINTF_SIMDATA_INTERVAL << " meaning every " << PRINTF_SIMDATA_INTERVAL * DELTA_T_SI << " s\n";

    //open log file
    FILE * log = NULL;
    log = fopen ("traj.dat","w");
    fprintf(log,"#t\trx\try\trz\tpx\tpy\tpz\tr\ttheta\t(Et-E)/E\t(Lt-L)/L"
                "\ttheta analyt.\t(th-th_anal/th_anal\n");
    assert(log != NULL);
    
	//some constants which only depend on the initial conditions
    Particle & p1 = simBox.electrons[0][0][0].front();
    Particle & p2 = simBox.ions     [0][0][0].front();
    Vec dir0      = p1.p / p1.p.norm();      //to remember initial direction of movement
	double theta0 = asin( ( Ew*L/( sqrt(2*p1.m*E) * r0  ) + 1 ) / sqrt( Ew*Ew+1 ) ) 
       - asin( 1/sqrt( Ew*Ew+1 ) ) -(M_PI - acos( (p1.r-p2.r)*Vec(1,0,0)/(p1.r-p2.r).norm() ) );
    tout << "L:" << L * UNIT_ANGULAR_MOMENTUM << " , E:" << E << " , epsilon:" << theta0 << "\n";
    
    // initialize some values we will calculate from the simulation
    double rmin_num = r0;
    double phi_min; //holds angle of position nearest to scattering center
    bool returning = false;
    
    // Initialize Verlet Integration (Leapfrog): shift momentum half a time step
    // this half step makes many calculated deviations wrong, because we prepare the initial state so carefully !!
    //simBox.Verlet( DELTA_T, true );
    // Main Simulation Loop
    for (uint32_t currentStep = 0; currentStep < NUMBER_OF_STEPS; currentStep++) {
        simBox.Verlet( DELTA_T );
        

        //calc theta, r
        Particle & p1 = simBox.electrons[0][0][0].front();
        Particle & p2 = simBox.ions     [0][0][0].front();
        double rt     = ( p1.r - p2.r ).norm();
        double Et     = simBox.CalcTotalEnergy();
        double Lt     = simBox.CalcTotalAngularMomentum().z;
        double thetat = acos( (p1.r-p2.r)/(p1.r-p2.r).norm() * Vec(1,0,0)  );
        if (rt < rmin_num) {
            rmin_num = rt;
            phi_min = thetat;
        } else if (!returning) {
            returning = true;
            tout << "Going back home now!\n"; //but calculate phi0 beforehand
            tout << "phi_min:" << phi_min*180/M_PI << "\n";
            tout << "rmin num.: " << rmin_num << " <-> " << rmin << " :rmin anal. => dev:" << (rmin_num-rmin)/rmin << "\n";
        }

        //calc deviation to analytical theta
        double theta      = acos( (p1.r-p2.r) / (p1.r-p2.r).norm() * Vec(1,0,0) );
        double phit       = asin( (Ew*L/(sqrt(2*p1.m*E)*rt)+1) / sqrt(Ew*Ew+1) );
        double theta_calc = M_PI + theta0 + theta_inf/2.0 - phit ;
        if (returning)
            theta_calc = theta0 + theta_inf/2.0 + phit;
            
        //Debug Output
        if ( currentStep % PRINTF_INTERVAL == 0 || rt > r0 ) {
            fprintf(log, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", currentStep*DELTA_T_SI, p1.r.x * UNIT_LENGTH, p1.r.y * UNIT_LENGTH, p1.r.z * UNIT_LENGTH, p1.p.x * UNIT_LENGTH, p1.p.y * UNIT_LENGTH, p1.p.z * UNIT_LENGTH, rt * UNIT_LENGTH, theta, (Et-simBox.E0)/simBox.E0, (Lt-simBox.L0.z)/simBox.L0.z, theta_calc, (theta-theta_calc)/theta_calc );
        }
        if ( currentStep % PRINT_INTERVAL == 0 ) {
            tout << "r" << rt << ", theta: num:" << theta*180/M_PI << " deg, anal:" << theta_calc*180/M_PI << "\n";
            tout << "    ( " << p1.r.x << " )        ( " << p1.p.x << " )\n";
            tout << "r = ( " << p1.r.y << " )  , p = ( " << p1.p.y << " )\n";
            tout << "    ( " << p1.r.z << " )        ( " << p1.p.z << " )\n";
        }
        
        /*
        struct return_data data = simulate_scattering( 2, particle, dt, "traj.dat", steps/200, steps );

        double dev_theta = (data.theta_max - data.theta_max_calc) / data.theta_max_calc;
        fprintf( stat, "%e\t%e\t%e\t%e\\n", particle[1].r.x, dt, steps, dev_theta );
        printf( "i:%i, j:%i => x:%e, dt:%e => devTheta:%e\n", i,j, particle[1].r.x, dt, dev_theta );
        */
        if ((currentStep % PRINTF_INTERVAL == 0) or (currentStep % PRINT_INTERVAL == 0)) {
            const Vec Pvec  = simBox.CalcTotalMomentum              () * UNIT_MOMENTUM + simBox.p * UNIT_MOMENTUM;
            const double Te = simBox.CalcTotalKineticEnergyElectrons() * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double Ti = simBox.CalcTotalKineticEnergyIons     () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double T  = simBox.CalcTotalKineticEnergy         () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double V  = simBox.CalcTotalPotentialEnergy       () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const Vec Pmax  = simBox.MaxMomentum                    () * UNIT_MOMENTUM;
            const double E  = T + V;
            const double P  = Pvec.x + Pvec.y + Pvec.z - (simBox.P0.x + simBox.P0.y + simBox.P0.z);
            const double dE = E - simBox.E0 * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            if (currentStep % PRINTF_INTERVAL == 0)
                fprintf( stats, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", currentStep * DELTA_T_SI, E, V, L, P, T, Te, Ti, dE);
            if (currentStep % PRINTF_SIMDATA_INTERVAL == 0)
                simBox.DumpData();
            if (currentStep % PRINT_INTERVAL == 0)
                tout << "[" << 100*(double)currentStep/(double)NUMBER_OF_STEPS << "%] E: " << E << ", V: " << V << ", <P.x>:" << Pvec.x/NUMBER_OF_PARTICLES << ", (P.x-P0)/Pmax.x:" << ( Pvec.x - simBox.P0.x * UNIT_MOMENTUM ) / Pmax.x /*max(max(Pmax.x,Pmax.y),Pmax.z)*/ << "\n";
        }
        
        if (simBox.numberOfAppliedBoundaryConditions > 0)
            break;
    }

    fclose(stats);
    fclose(log);

    return 0;
}

