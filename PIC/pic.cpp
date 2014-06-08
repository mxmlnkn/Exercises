#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <iostream>

#define DEBUG 1

using namespace std;

#include "Vector.cpp"
#include "Parameters.cpp"

struct Particle {
    Vec r;   //m
    Vec p;   //kg m/s
    double m;       //kg
    double q;       //C
};

        Vec force( Particle particle[] ) {
            double r = (particle[1].r - particle[0].r).norm();
            return (particle[0].r - particle[1].r) / r * (-1) * particle[0].q * particle[1].q / ( 4 * M_PI * EPS0 * r*r );
        }

        Vec RK4(Particle particle[], const double dt) {
            /*const double dz = v * dt;
            const double a1 = force(z,v);
            const double a2 = force(z+dz/2., v + a1*dt/2.);
            const double a3 = force(z+dz/2., v + a2*dt/2.);
            const double a4 = force(z+dz   , v + a3*dt   );
            const double a = (a1+ 2*a2+ 2*a3+ a4)/6.0;
            z+=v*dt + a/2*dt*dt;                                    // new location (s=a/2t^2+v*t+z)
            v+=a*dt;                                                // new velocity
            return a;*/
            Vec dr = particle[1].p / particle[1].m * dt;
            Vec F1 = force( particle );
                particle[1].r += dr/2.0;
                particle[1].p += F1*dt/2.0;
            Vec F2 = force( particle );
                particle[1].p -= F1*dt/2.0;
                particle[1].p += F2*dt/2.0;
            Vec F3 = force( particle );
                particle[1].r -= dr/2.0;
                particle[1].r += dr;
                particle[1].p -= F2*dt/2.0;
                particle[1].p += F3*dt;
            Vec F4 = force( particle );
            Vec F  = ( F1 + 2*F2 + 2*F3 + F4 ) / 6.0;
            particle[1].p += F * dt;
            particle[1].r += F / particle[1].m * dt * dt / 2 + particle[1].p / particle[1].m * dt;
            return F;
        }

        Vec Euler(Particle particle[], const double dt) {
            //v += a*dt;
            //s += a/2 t^2 + v*t
            double r = (particle[0].r - particle[1].r).norm();
            Vec F = (particle[0].r - particle[1].r) / r * (-1.0) * particle[0].q * particle[1].q / ( 4 * M_PI * EPS0 * r*r ) ;
            particle[1].p += F * dt;
            particle[1].r += F / particle[1].m * dt * dt / 2 + particle[1].p / particle[1].m * dt;
            return F;
        }

        Vec Own(Particle particle[], const double dt) { //is this Verlet ???
            particle[1].p += force(particle) * dt/2.0 ;
            particle[1].r += particle[1].p / particle[1].m * dt/2.0;
            particle[1].r += particle[1].p / particle[1].m * dt/2.0;
            particle[1].p += force(particle) * dt/2.0 ;
            /*
            //Calculate F
            double r = (particle[0].r - particle[1].r).norm();
            Vec F = (particle[0].r - particle[1].r) / r * (-1.0) * 
                     particle[0].q * particle[1].q / ( 4 * M_PI * EPS0 * r*r );
			#if DEBUG == 1
				if (F != F) {
					printf("r:%e\n",r );
					printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.x, particle[1].p.x, F.x );
					printf("r = ( %e )  , p = ( %e )  , F = ( %e )\n", particle[1].r.y, particle[1].p.y, F.y );
					printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.z, particle[1].p.z, F.z );
				}
				assert(F == F);
            #endif
            
            particle[1].p += F * dt ;*/
            return force(particle);
        }

double CalcTotalEnergy( const unsigned int particleCount, struct Particle particle[] ) {
    double E = 0;
    for (unsigned int i = 0; i < particleCount; i++) {
        double p = particle[i].p.norm();
        double T = p*p / (2.0 * particle[i].m);
        double V = 0;
        for (unsigned int j = 0; j < i; j++)
            V += particle[i].q * particle[j].q / ( 4 * M_PI * EPS0 * (particle[i].r - particle[j].r).norm() );
        E += T + V;
    }
    return E;
}

double CalcAngularMomentum( const unsigned int particleCount, struct Particle particle[] ) {
    Vec L;
    for (unsigned int i = 0; i < particleCount; i++) {
        L.x += particle[i].r.y * particle[i].p.z - particle[i].r.z * particle[i].p.y;
        L.y += particle[i].r.z * particle[i].p.x - particle[i].r.x * particle[i].p.z;
        L.z += particle[i].r.x * particle[i].p.y - particle[i].r.y * particle[i].p.x;
    }
    return sqrt( L.x*L.x + L.y*L.y + L.z*L.z );
}


struct return_data {
    double theta_max, theta_max_calc, rmin, rmin_calc;
};

struct return_data simulate_scattering(const unsigned int particleCount, Particle particle[], const double dt, bool writeToFile = false, bool quiet = true){
    //Calculate some conserved observables
    double E    = CalcTotalEnergy( particleCount, particle );
    double L    = CalcAngularMomentum( particleCount, particle );
    double r0   = (particle[0].r - particle[1].r).norm();
    Vec dir0       = particle[1].p / particle[1].p.norm();      //to remember initial direction of movement
    struct return_data data;
	data.theta_max = 0;
	data.rmin      = r0;
    
	//some constants which only depend on the initial conditions
	double alpha  = particle[0].q * particle[1].q / ( 4 * M_PI * EPS0 );
	double Ew     = sqrt( 2*E/particle[1].m ) * L/alpha;
	double theta0 = asin( ( Ew*L/( sqrt(2*particle[1].m*E) * r0  ) + 1 ) / sqrt( Ew*Ew+1 ) ) - asin( 1/sqrt( Ew*Ew+1 ) );
	
    printf("L:%e, E:%e, epsilon:%E\n", L*UNIT_ANGULAR_MOMENTUM, E*UNIT_ENERGY, theta0);
	
    //open log file
    FILE * log = NULL;
    if (writeToFile) {
        log = fopen ("traj.dat","w");
        fprintf(log,"#t\trx\try\trz\tpx\tpy\tpz\tr\ttheta\t(Et-E)/E\t(Lt-L)/L\ttheta_calc\t(theta-theta_calc)/theta_calc\n");
        assert(log != NULL);
    }

    const double t_end = 2* fabs(particle[1].r.x) / sqrt( 2*E/particle[1].m );
    int print_interval    = t_end / dt / 10;
    int fprint_interval   = t_end / dt / 200;

	//Initialize Verlet Integration (Leapfrog)
	Vec F = force(particle);
	particle[1].p += F*dt/2.0; 
	
    bool returning = false;
    int i = 0;
    do {
        //push particles
        //for (int i = 0; i < particleCount-1; i++) {
            Vec F = RK4( particle, dt );
            //Verlet:
        //}

        //calc theta, r
        double r  = (particle[0].r - particle[1].r).norm();
        //double r2 = r*r;
        double Et = CalcTotalEnergy( particleCount, particle );
        double Lt = CalcAngularMomentum( particleCount, particle );
        if (r < data.rmin) data.rmin = r;
        else if (!returning) {
            returning = true;
            printf("Going back home now!\n"); //but calculate phi0 beforehand
        }

        //calc deviation to analytical theta
        //double theta = acos( (particle[1].p / particle[1].p.norm()) * dir0 )/M_PI*180;
        double theta = acos( (particle[1].r / particle[1].r.norm()) * Vec(1,0,0) )/M_PI*180;
        double theta_calc = M_PI + asin( ( 0*Ew*L/( sqrt(2*particle[1].m*E) * r0 ) + 1 ) / sqrt( Ew*Ew+1 ) )
                                 - asin( (   Ew*L/( sqrt(2*particle[1].m*E) * r  ) + 1 ) / sqrt( Ew*Ew+1 ) );
        if (returning) {
            double theta_inf = 2*asin( 1/sqrt( 1.0 + Ew*Ew ) ) ;
            //theta_inf = 164.1842/180.0*M_PI;	//numerical result... 
            theta_calc = M_PI - theta_calc + theta_inf;// + theta0;
        }
        theta_calc *= 180/M_PI;
        data.theta_max = theta;
        if ( ( i % fprint_interval == 0 || r > r0 ) && writeToFile) {
            fprintf(log, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
            i*dt*UNIT_TIME,
            particle[1].r.x*UNIT_LENGTH,
            particle[1].r.y*UNIT_LENGTH,
            particle[1].r.z*UNIT_LENGTH,
            particle[1].p.x*UNIT_MOMENTUM,
            particle[1].p.y*UNIT_MOMENTUM,
            particle[1].p.z*UNIT_MOMENTUM,
            r*UNIT_LENGTH, theta, (Et-E)/E, (Lt-L)/L, theta_calc, (theta-theta_calc)/theta_calc );
        }
        if ( ( i % print_interval == 0 || r > r0 ) && !quiet) {
            printf("r:%e, theta: num:%e deg, anal:%e, F:%e\n",r*UNIT_LENGTH, theta, theta_calc, F.norm()*UNIT_FORCE );
            printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.x*UNIT_LENGTH, particle[1].p.x*UNIT_MOMENTUM, F.x*UNIT_FORCE );
            printf("r = ( %e )  , p = ( %e )  , F = ( %e )\n", particle[1].r.y*UNIT_LENGTH, particle[1].p.y*UNIT_MOMENTUM, F.y*UNIT_FORCE );
            printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.z*UNIT_LENGTH, particle[1].p.z*UNIT_MOMENTUM, F.z*UNIT_FORCE );
        }

        i++;
    } while( (particle[0].r - particle[1].r).norm() < r0 );
    printf ("%e = r > r0 = %e\n" ,(particle[0].r - particle[1].r).norm()*UNIT_LENGTH ,r0*UNIT_LENGTH);

    double r = (particle[0].r - particle[1].r).norm();
    
    data.rmin_calc = alpha/(2*E)*( 1 + sqrt( 1 + Ew*Ew ) );
    //double theta_max_calc = 2*asin( 1/sqrt( 1+2*E*L*L/(particle[1].m*alpha*alpha) ) ) ;
    double corrector =  1.0 - ( Ew*L/( sqrt(2*particle[1].m*E) * data.rmin_calc  ) + 1 ) /sqrt( Ew*Ew+1 ) ;
    printf("Correction Term: %E deg\n", corrector);
    data.theta_max_calc = 2*asin( ( Ew*L/(sqrt(2*particle[1].m*E)*r) + 1 )/sqrt( Ew*Ew+1 ) ); //2*asin( 1/sqrt( 1+Ew*Ew ) );
    data.theta_max_calc *= 180/M_PI;
    //if (!quiet) {
        printf("rmin num.:     %E <-> %E :rmin anal. => dev:%e\n", data.rmin*UNIT_LENGTH, data.rmin_calc*UNIT_LENGTH, (data.rmin - data.rmin_calc) / data.rmin_calc );
        printf("thetamax num.: %E <-> %E :anal.      => dev:%e\n", data.theta_max, data.theta_max_calc, (data.theta_max - data.theta_max_calc) / data.theta_max_calc );
    //}

    if (writeToFile)
        fclose(log);
    return data;
}

int main(void) {
    const double test = test / 3.; //wth ... this actually works Oo

    printf("MUE0: %e\n",MUE0);
    printf("EPS0: %e\n",EPS0);
    printf("SPEED_OF_LIGHT: %e\n",SPEED_OF_LIGHT);
    printf("CELL_SIZE_SI: %e\n",CELL_SIZE_SI);
    printf("CELL_SIZE: %e\n",CELL_SIZE_SI / UNIT_LENGTH);
    printf("NUMBER_OF_CELLS_X: %d\n",NUMBER_OF_CELLS_X);
    printf("DELTA_T_SI: %e\n",DELTA_T_SI);
    printf("UNIT_ENERGY: %e\n",UNIT_ENERGY);
    printf("ELECTRON_MASS: %e\n",ELECTRON_MASS);
    
    
    srand(time(NULL));
    
    Particle electrons[NUMBER_OF_PARTICLES];
    for (uint32_t i = 0; i < NUMBER_OF_PARTICLES; i++) {
        //Particle to be scattered
        electrons[i].m   = ELECTRON_MASS;
        electrons[i].q   = ELECTRON_CHARGE;
        electrons[i].r.x = (SIMDIM > 0) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_X * CELL_SIZE_X;
        electrons[i].r.y = (SIMDIM > 1) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Y * CELL_SIZE_Y;
        electrons[i].r.z = (SIMDIM > 2) * (double)rand()/RAND_MAX * NUMBER_OF_CELLS_Z * CELL_SIZE_Z;
        electrons[i].p.x = 0;
        electrons[i].p.y = 0;
        electrons[i].p.z = 0;
    }

    
    /* Because of this internally units ARE needed! Calculating with SI the whole time uses up unneeded exponent "space" */
    /*printf("Max Float:%e, Float Min:%e\n", FLT_MAX, FLT_MIN);
    printf("Max Float:%i, Float Min:%i\n", FLT_MAX_EXP, FLT_MIN_EXP);
    float fer = 1.00e-42;   //why the fuck is this working? => "FLT_MAX_EXP: The maximal exponent of a floating point value expressed in base FLT_RADIX; greater exponents are principally possible (up to 16383), but not supported in all math functions."
    float rec = 1.0/fer;
    printf("1.0/%e = %e\n", fer, rec);
    // FURTHER PROBLEM HERE! somehow printf shows the right number, but only if %e instead of %f is used. Is there some kind of buffer overflow happning in both the assignment of "fer" AND the readout in the printf function Oo ???
    return 0; */
    
    double m      = ELECTRON_MASS;
    double q1     = ELECTRON_CHARGE;
    double q2     = ELECTRON_CHARGE;
    double s      = 250e-15 / UNIT_LENGTH;	//collision parameter
    double p_SI   = sqrt( 2.*m*UNIT_MASS* 1.e3 * ELEMENTARY_CHARGE_SI ); //initial momentum in infinity in x direction, here: T=10Mev=>p
    double p      = p_SI / UNIT_MOMENTUM;
    double E      = p*p/(2.*m);
    double L      = p*s;
    double alpha  = q1*q2 / ( 4.*M_PI*EPS0 );
	double Ew     = sqrt(2.*E/m) * L/alpha;
	double rmin   = alpha/(2.*E)*( 1. + sqrt( 1. + Ew*Ew ) );
    printf("E:%E, L:%E\n", E*UNIT_ENERGY, L*UNIT_ANGULAR_MOMENTUM);
    
    //statistics for more than one scattering event
    FILE * stat = fopen( "statistics.dat", "w" );
    fprintf( stat, "#x\tdt\tsteps\tdev theta\n" );

    for (int i = 6; i <= 6; i++) {
        for (int j = 1; j <= 1; j++) {
			double r0   = rmin * pow(10,j); //must be greater than rmin!!!
			double p0   = sqrt(2.*m) * sqrt( E - alpha/r0 );
			assert( r0 >= rmin ); //else asin will fail !
			double phi0 = M_PI + asin( 1. / sqrt(Ew*Ew+1.) )
                               - asin( ( Ew*L/( sqrt(2.*m*E)*r0 ) + 1. ) / sqrt(Ew*Ew+1.) );
			
			//TODO ANGABE VON r0 veraendert streuwinkel!!! shouldn't happen! -> Calculation of vec p wrong?
			
            Particle particle[2];
            //Particle to be scattered
            particle[1].m   = m;
            particle[1].q   = q2;
            /* r.x * 10 => deviation of theta / 10 approx. Shouldn't be because       *
             * of wrongly calculated p vs. p_inf, because p_inf is correctly          *
             * calculated from total energy. Therefore it should be a geometrical     *
             * problem, because we want to know the asymptotes of the trajectories.   *
             * Is also almost independent of Integration method for 1e5 or more steps *
             * p may be correct, but the direction of p could still be error prone !  */
            particle[1].r.x =-r0*cos(M_PI - phi0);
            particle[1].r.y = r0*sin(M_PI - phi0);
            particle[1].p.x = p*s/(r0*r0) * ( particle[1].r.y - particle[1].r.x * sqrt( pow( r0*p0/(p*s) ,2) - 1 ) );   // calculate from initial momentum and collision parameter using momentum and energy conservation
            particle[1].p.y = sqrt(p0*p0 - particle[1].p.x * particle[1].p.x);
            //particle 0 is static scattering center, won't move (!!!)
            particle[0].m   = PROTON_MASS_SI;
            particle[0].q   = q1;
            printf( "phi0:%E, p0:%E, p_inf=%E\n", (M_PI-phi0)/M_PI*180, p0*UNIT_MOMENTUM,p*UNIT_MOMENTUM );
            printf("     ( %e )         ( %e )\n", particle[1].r.x*UNIT_LENGTH, particle[1].p.x*UNIT_MOMENTUM );
            printf("r0 = ( %e )  , p0 = ( %e )\n", particle[1].r.y*UNIT_LENGTH, particle[1].p.y*UNIT_MOMENTUM );
            printf("     ( %e )         ( %e )\n", particle[1].r.z*UNIT_LENGTH, particle[1].p.z*UNIT_MOMENTUM );
            
            E = CalcTotalEnergy(2,particle);
            L = CalcAngularMomentum(2,particle);
            printf( "E:%E, L:%E\n", E*UNIT_ENERGY, L*UNIT_ANGULAR_MOMENTUM );

            //just a rough guess for the endtime to calculate dt. Analytically calculatable?!
            double E    = CalcTotalEnergy( 2, particle );
            const double t_end   = 2* fabs(particle[1].r.x) / sqrt( 2*E/particle[1].m );
            printf("tend:%e\n" , t_end*UNIT_TIME);
            double steps = pow(10,i);
            const double dt = t_end / steps;

            struct return_data data = simulate_scattering( 2, particle, dt, true, false );
            double dev_theta = (data.theta_max - data.theta_max_calc) / data.theta_max_calc;
            fprintf( stat, "%e\t%e\t%e\t%e\\n", particle[1].r.x*UNIT_LENGTH, dt*UNIT_TIME, steps, dev_theta );
            printf( "i:%i, j:%i => x:%e, dt:%e => devTheta:%e\n", i,j, particle[1].r.x*UNIT_LENGTH, dt*UNIT_TIME, dev_theta );
        }
    }

    fclose(stat);
    return 0;
}
