// del pic.exe && g++ pic.cpp -o pic.exe -Wall -Wextra -Wimplicit -pedantic-errors -pedantic -Wcomment -Wconversion -std=c++0x && pic.exe
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

#define M_PI 3.14159265358979323846

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

class SimulationBox {
public:
    Vec p;
    SimulationBox();
    ~SimulationBox();
    double E0, Te0, Ti0;
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
    Vec CalcTotalMomentum( void );
    Vec MaxMomentum( void );
    double CalcTotalAngularMomentum( void );
    
    void PushLocation( const double dt );
    void ClearDp( void );
    void CalcDp ( const double dt );
    void ApplyDp( void );
    void Verlet ( const double dt, bool initialize = false );
    bool CheckBoundaryConditions( Particle & particle );
};

SimulationBox::SimulationBox() : p(Vec(0,0,0)) {

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
    #if AWESOME_REFLEXION == 1
        if (outOfBounds > 1 and BOUNDARY_CONDITION == 1)
            tout << "Double Reflexion at one of the corners or edges => Trajectory non-physical!\n";
    #endif
    return outOfBounds > 0;
}

Vec Velocity( Vec p, double m ) {
    // for the non-relativistic limit the result will be p/m
    return p / sqrt( m*m + p.norm2() / ( SPEED_OF_LIGHT * SPEED_OF_LIGHT ) );
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

double SimulationBox::CalcTotalAngularMomentum( void ) {
    Vec L(0,0,0);
    /*for (unsigned int i = 0; i < particleCount; i++) {
        L.x += particle[i].r.y * particle[i].p.z - particle[i].r.z * particle[i].p.y;
        L.y += particle[i].r.z * particle[i].p.x - particle[i].r.x * particle[i].p.z;
        L.z += particle[i].r.x * particle[i].p.y - particle[i].r.y * particle[i].p.x;
    }*/
    return L.norm();
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

    list<int> L;
    L.push_back(0);              // Insert a new element at the end
    L.push_front(0);             // Insert a new element at the beginning
    L.insert(++L.begin(),2);     // Insert "2" before position of first argument
                                 // (Place before second argument)
    L.push_back(5);
    L.push_back(6);

    list<int>::iterator i;

    for(i=L.begin(); i != L.end(); i++)
        cout << *i << " ";
    cout << endl;
    int k = 0;
    for(i=L.begin(); i != L.end(); i++) {
        if (k++ == 2)
            L.erase(i++);
        cout << *i << " ";
    }
    cout << endl;
    for(i=L.begin(); i != L.end(); i++)
        cout << *i << " ";
    cout << endl;

  const int nrolls=10000;  // number of experiments
  const int nstars=100;    // maximum number of stars to distribute

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0,2.0);

  int p[10]={};

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    if ((number>=0.0)&&(number<10.0)) ++p[int(number)];
  }

  std::cout << "normal_distribution (5.0,2.0):" << std::endl;

  for (int i=0; i<10; ++i) {
    std::cout << i << "-" << (i+1) << ": ";
    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
  }


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
            if ( SPECIES == 3 && i >= NUMBER_OF_PARTICLES / 2.)
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

        /*double E = CalcTotalPotentialEnergy( NUMBER_OF_PARTICLES, electrons, shape ) /
                   ( NUMBER_OF_PARTICLES_PER_CELL*NUMBER_OF_CELLS_X*NUMBER_OF_CELLS_Y*NUMBER_OF_CELLS_Z );
        bool inRange = binCellEnergies.addData( E*UNIT_ENERGY*UNITCONV_Joule_to_keV );
        if (!inRange) energiesOutsideRange++; */

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

    // To find the lowest Energy configuration: adjust Verlet with p *= 0.9999

    // Open log file
    FILE * stats = NULL;
    stats = fopen ("stats.dat","w");
    fprintf( stats, "# t\tE/keV\tV/keV\tL/(m*kg*m/s)\tP/(kg*m/s)\tT/keV\tTe/keV\tTi/keV(E-E0)/keV\n" );
    assert(stats != NULL);
    
    for (uint32_t ix = 0; ix < NUMBER_OF_CELLS_X; ix++)
    for (uint32_t iy = 0; iy < NUMBER_OF_CELLS_Y; iy++)
    for (uint32_t iz = 0; iz < NUMBER_OF_CELLS_Z; iz++) {
    list<Particle> & eons = simBox.electrons[ix][iy][iz];
    list<Particle> & ions = simBox.ions     [ix][iy][iz];
        std::default_random_engine rng(RANDOM_SEED);
        const double initialDrift = sqrt(1. - PARTICLE_INIT_DRIFT_GAMMA * PARTICLE_INIT_DRIFT_GAMMA) * SPEED_OF_LIGHT; // g*g = 1-v*v/c*c
        // draw momentum from temperature distribution. Of course this only makes sense non-relativistically
        std::normal_distribution<double> eon_vel_dist( initialDrift, sqrt( ELECTRON_TEMPERATURE / ELECTRON_MASS ) );
        std::normal_distribution<double> ion_vel_dist( initialDrift, sqrt(      ION_TEMPERATURE / ION_MASS      ) );
        for (uint32_t i = 0; i < NUMBER_OF_PARTICLES_PER_CELL / 2.0; i++) {
            Particle particle;
            particle.m   = ELECTRON_MASS;
            particle.q   = ELECTRON_CHARGE;
            particle.r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
            particle.r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
            particle.r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
            particle.p.x = (SIMDIM > 0) * eon_vel_dist(rng) * particle.m;
            particle.p.y = (SIMDIM > 1) * eon_vel_dist(rng) * particle.m;
            particle.p.z = (SIMDIM > 2) * eon_vel_dist(rng) * particle.m;
            eons.push_back( particle );
        }
        for (uint32_t i = 0; i < NUMBER_OF_PARTICLES_PER_CELL / 2.0; i++) {
            Particle particle;
            particle.m   = ION_MASS;
            particle.q   = ION_CHARGE;
            particle.r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
            particle.r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
            particle.r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
            particle.p.x = (SIMDIM > 0) * ion_vel_dist(rng) * particle.m;
            particle.p.y = (SIMDIM > 1) * ion_vel_dist(rng) * particle.m;
            particle.p.z = (SIMDIM > 2) * ion_vel_dist(rng) * particle.m;
            ions.push_back( particle );
        }
    }
    simBox.Te0 = simBox.CalcTotalKineticEnergyElectrons();
    simBox.Ti0 = simBox.CalcTotalKineticEnergyIons();
    simBox.E0  = simBox.Te0 + simBox.Ti0 + simBox.CalcTotalPotentialEnergy();
    simBox.P0  = simBox.CalcTotalMomentum();
    tout << "Initial total   Energy P0.x : " << simBox.P0.x * UNIT_MOMENTUM << " kg m/s\n";
    tout << "Initial total   Energy E0   : " << simBox.E0   * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Initial kinetic Energy Te0  : " << simBox.Te0  * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Initial kinetic Energy Ti0  : " << simBox.Ti0  * UNIT_ENERGY * UNITCONV_Joule_to_keV << "keV\n";
    tout << "Note that <E_kin_e> = Te0/N_Particles = " << simBox.Te0 / (NUMBER_OF_PARTICLES/2.0) * UNIT_ENERGY * UNITCONV_Joule_to_keV * 1000 << "eV = f/2 kT, where f=" << SIMDIM << "\n";
    tout << "Note that <E_kin_i> = Ti0/N_Particles = " << simBox.Ti0 / (NUMBER_OF_PARTICLES/2.0) * UNIT_ENERGY * UNITCONV_Joule_to_keV * 1000 << "eV = f/2 kT, where f=" << SIMDIM << "\n";
    
    tout << "Stats fprintf interval: " << PRINTF_SIMDATA_INTERVAL << " meaning every " << PRINTF_SIMDATA_INTERVAL * DELTA_T_SI << " s\n";

    // Initialize Verlet Integration (Leapfrog): shift momentum half a time step
    simBox.Verlet( DELTA_T, true );
    // Main Simulation Loop
    for (uint32_t currentStep = 0; currentStep < NUMBER_OF_STEPS; currentStep++) {
        simBox.Verlet( DELTA_T );

        if ((currentStep % PRINTF_INTERVAL == 0) or (currentStep % PRINT_INTERVAL == 0)) {
            const Vec Pvec  = simBox.CalcTotalMomentum              () * UNIT_MOMENTUM + simBox.p * UNIT_MOMENTUM;
            const double Te = simBox.CalcTotalKineticEnergyElectrons() * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double Ti = simBox.CalcTotalKineticEnergyIons     () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double T  = simBox.CalcTotalKineticEnergy         () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double V  = simBox.CalcTotalPotentialEnergy       () * UNIT_ENERGY * UNITCONV_Joule_to_keV;
            const double L  = simBox.CalcTotalAngularMomentum       () * UNIT_ANGULAR_MOMENTUM;
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
    }

    fclose(stats);

    return 0;
}

