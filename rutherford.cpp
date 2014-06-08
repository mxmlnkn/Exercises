#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

#define DEBUG 1

class Vec {
public:
    double x,y,z;
    Vec(void) : x(0),y(0),z(0) { }
    Vec(const Vec & v) : x(v.x),y(v.y),z(v.z) { }
    Vec(const double x, const double y, const double z) : x(x),y(y),z(z) { }
    Vec& operator+= (const Vec & v);
    Vec operator+ (const Vec & v) const;
    Vec& operator-= (const Vec & v);
    Vec operator- (const Vec & v) const;
    double operator* (const Vec v) const;
    Vec& operator*= (const double a);
    Vec operator* (const double a) const;
    Vec& operator/= (const double a);
    Vec operator/ (const double a) const;
    bool operator== (const Vec v) const;
    bool operator!= (const Vec v) const;
    double norm() const;
};
double Vec::operator* (const Vec v) const {
    return this->x * v.x + this->y * v.y + this->z * v.z;
}

Vec& Vec::operator+= (const Vec & v) {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}
Vec Vec::operator+ (const Vec & v) const {
    Vec res = *this;
    res += v;
    return res;
}

Vec& Vec::operator-= (const Vec & v) {
    return (*this) += v * (-1);
}
Vec Vec::operator- (const Vec & v) const {
    return (*this)+( v*(-1) );
}

Vec& Vec::operator*= (const double a) {
    this->x *= a;
    this->y *= a;
    this->z *= a;
    return *this;
}
Vec Vec::operator* (const double a) const {
    Vec res = *this;
    res *= a;
    return res;
}

Vec& Vec::operator/= (const double a) {
    (*this) *= 1.0/a;
    return *this;
}
Vec Vec::operator/ (const double a) const {
    return (*this)*(1.0/a);
}

bool Vec::operator== (const Vec v) const {
    return ( (this->x == v.x) && (this->y == v.y) && (this->z == v.z) );
}
bool Vec::operator!= (const Vec v) const {
    return !((*this) == v);
}
double Vec::norm() const {
    assert( this->x == this->x );
    assert( this->y == this->y );
    assert( this->z == this->z );
    double res = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    assert( res == res );
    return res;
}


template<typename T>
Vec operator*(T const& scalar, Vec rhs)
{
    // scalar multiplication is commutative: s M = M s
    return rhs *= scalar; // calls rhs.operator*=(scalar);
}


struct Particle {
    Vec r;   //m
    Vec p;   //kg m/s
    double m;       //kg
    double q;       //C
};

#define epsilon0           8.854187817e-12  // C^2/J*m
#define elementary_charge  1.602176565e-19  // C
#define m_e                9.10938291e-31   // kg
#define speed_of_light     2.99792458e8     // m/s
//#define M_PI               3.14159265358979323846
#define m_p                1.6726231e-27    // kg

        Vec force( Particle particle[] ) {
            double r = (particle[1].r - particle[0].r).norm();
            return (particle[0].r - particle[1].r) / r * (-1) * particle[0].q * particle[1].q / ( 4 * M_PI * epsilon0 * r*r );
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
            Vec F = (particle[0].r - particle[1].r) / r * (-1.0) * particle[0].q * particle[1].q / ( 4 * M_PI * epsilon0 * r*r ) ;
            particle[1].p += F * dt;
            particle[1].r += F / particle[1].m * dt * dt / 2 + particle[1].p / particle[1].m * dt;
            return F;
        }

        Vec Verlet(Particle particle[], const double dt) {
            particle[1].p += force(particle) * dt/2.0 ;
            particle[1].r += particle[1].p / particle[1].m * dt/2.0;
            particle[1].r += particle[1].p / particle[1].m * dt/2.0;
            particle[1].p += force(particle) * dt/2.0 ;
            /*
            //Calculate F
            double r = (particle[0].r - particle[1].r).norm();
            Vec F = (particle[0].r - particle[1].r) / r * (-1.0) *
                     particle[0].q * particle[1].q / ( 4 * M_PI * epsilon0 * r*r );
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
            V += particle[i].q * particle[j].q / ( 4 * M_PI * epsilon0 * (particle[i].r - particle[j].r).norm() );
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

struct return_data simulate_scattering(const unsigned int particleCount, Particle particle[], const double dt, const char* filename = "", unsigned int fprint_interval = 0 , unsigned int print_interval = 0) {
    bool writeToFile = strlen(filename) && fprint_interval;
    bool quiet       = !print_interval;

    //Calculate some conserved observables
    double E         = CalcTotalEnergy( particleCount, particle );
    double L         = CalcAngularMomentum( particleCount, particle );
    double r0        = (particle[0].r - particle[1].r).norm();
    Vec dir0         = particle[1].p / particle[1].p.norm();      //to remember initial direction of movement
    double phi_min; //holds angle of position nearest to scattering center
    
	//some constants which only depend on the initial conditions
	double alpha     = particle[0].q * particle[1].q / ( 4 * M_PI * epsilon0 );
	double Ew        = sqrt( 2*E/particle[1].m ) * L/alpha;
	double theta0    = asin( ( Ew*L/( sqrt(2*particle[1].m*E) * r0  ) + 1 ) / sqrt( Ew*Ew+1 ) ) 
       - asin( 1/sqrt( Ew*Ew+1 ) ) -(M_PI - acos( particle[1].r*Vec(1,0,0)/particle[1].r.norm() ) );
    double theta_inf = 2*asin( 1/sqrt( 1.0 + Ew*Ew ) ); //scattering angle for infinity!
    double rmin_calc = alpha/(2*E)*( 1 + sqrt( 1 + Ew*Ew ) );
    if (!quiet)
        printf("L:%e, E:%e, epsilon:%E\n", L,E,theta0);

    //Prepare return packet
    struct return_data data;
        data.rmin           = r0;
        data.rmin_calc      = rmin_calc;
        /* the error on theta_max is fairly large, because the calculated    *
         * theta assumes a particle to come from infinity and leave to it.   *
         * The error on theta_max therefore is the error of the asymptote!   *
         * For asymptote corrected calculated thetas, see the theta which is *
         * written into the log file and there the last entry, which is      *
         * equivalent to theta max.                                          */
        data.theta_max_calc = theta_inf;
    
    //open log file
    FILE * log = NULL;
    if (writeToFile) {
        log = fopen ("traj.dat","w");
        fprintf(log,"#t\trx\try\trz\tpx\tpy\tpz\tr\ttheta\t(Et-E)/E\t(Lt-L)/L"
                    "\ttheta analyt.\t(th-th_anal/th_anal\n");
        assert(log != NULL);
    }

    bool returning = false;
    int i = 0;
    do {
        //push particles
        Vec F = Verlet( particle, dt );

        //calc theta, r
        double r  = (particle[0].r - particle[1].r).norm();
        double Et = CalcTotalEnergy( particleCount, particle );
        double Lt = CalcAngularMomentum( particleCount, particle );
        if (r < data.rmin) {
            data.rmin = r;
            phi_min = acos( particle[1].r * Vec(1,0,0) / particle[1].r.norm() );
        }
        else if (!returning) {
            returning = true;
            printf("phi_min:%E\n",phi_min*180/M_PI);
            printf("Going back home now!\n"); //but calculate phi0 beforehand
        }

        //calc deviation to analytical theta
        double theta      = acos( (particle[1].r / particle[1].r.norm()) * Vec(1,0,0) );
        double theta_inf  = 2*asin( 1/sqrt( 1.0 + Ew*Ew ) ); //scattering angle for infinity!
        double phi_r      = asin( (Ew*L/(sqrt(2*particle[1].m*E)*r)+1) / sqrt(Ew*Ew+1) );
        double theta_calc = M_PI + theta0 + theta_inf/2.0 - phi_r ;
        if (returning)
            theta_calc = theta0 + theta_inf/2.0 + phi_r;
            
        //Debug Output
        if ( ( i % fprint_interval == 0 || r > r0 ) && writeToFile) {
            fprintf(log, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i*dt, particle[1].r.x,
            particle[1].r.y, particle[1].r.z, particle[1].p.x, particle[1].p.y,
            particle[1].p.z, r, theta, (Et-E)/E, (Lt-L)/L, theta_calc, (theta-theta_calc)/theta_calc );
        }
        if ( !quiet && ( i % print_interval == 0 || r > r0 ) ) {
            printf("r:%e, theta: num:%e deg, anal:%e, F:%e\n",r, theta, theta_calc, F.norm() );
            printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.x, particle[1].p.x, F.x );
            printf("r = ( %e )  , p = ( %e )  , F = ( %e )\n", particle[1].r.y, particle[1].p.y, F.y );
            printf("    ( %e )        ( %e )        ( %e )\n", particle[1].r.z, particle[1].p.z, F.z );
        }

        i++;
    } while( (particle[0].r - particle[1].r).norm() < r0 );

    //data.theta_max = acos( (particle[1].p / particle[1].p.norm()) * dir0 );
    data.theta_max = acos( (particle[1].r / particle[1].r.norm()) * Vec(1,0,0) );
    if (!quiet) {
        printf("rmin num.:     %E <-> %E :rmin anal. => dev:%e\n", data.rmin, data.rmin_calc, (data.rmin - data.rmin_calc) / data.rmin_calc );
        printf("thetamax num.: %E <-> %E :anal.      => dev:%e\n", data.theta_max*180/M_PI, data.theta_max_calc*180/M_PI, (data.theta_max - data.theta_max_calc) / data.theta_max_calc );
        
        double thm2 = M_PI - 2*(M_PI - phi_min) + theta0;
        printf("theta_max from epsilon derivation and minimal distance: %E -> dev:%E\n", thm2*180/M_PI, (thm2-data.theta_max)/data.theta_max );
    }

    if (writeToFile)
        fclose(log);
    return data;
}

int main(void) {
    /*srand(time(NULL));
    for (int i=0; i<particleCount; i++) {
        particle[i].p.x = rand()/(float)RAND_MAX;
    }*/

    double m      = m_e;
    double q1     = -elementary_charge, q2 = -elementary_charge;
    double s      = 250e-15;	//collision parameter
    double p      = sqrt( 2*m* 1e3 * elementary_charge ); //initial momentum in infinity in x direction, here: T=10Mev=>p
    double E      = p*p/(2*m);
    double L      = p*s;
    double alpha  = q1*q2 / ( 4*M_PI*epsilon0 );
	double Ew     = sqrt(2*E/m) * L/alpha;
	double rmin   = alpha/(2*E)*( 1 + sqrt( 1 + Ew*Ew ) );
    printf("E:%E, L:%E\n", E,L);

    //statistics for more than one scattering event
    FILE * stat = fopen( "statistics.dat", "w" );
    fprintf( stat, "#x\tdt\tsteps\tdev theta\n" );

    for (int i = 6; i <= 6; i++) {
        for (int j = 1; j <= 1; j++) {
			double r0   = rmin * pow(10,j); //must be greater than rmin!!!
			double p0   = sqrt(2*m) * sqrt( E - alpha/r0 );
			assert( r0 >= rmin ); //else asin will fail !
			double phi0 = M_PI + asin( 1 / sqrt(Ew*Ew+1) )
                               - asin( ( Ew*L/( sqrt(2*m*E)*r0 ) + 1 ) / sqrt(Ew*Ew+1) );

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
            particle[1].r.x = -r0*cos(M_PI - phi0);
            particle[1].r.y =  r0*sin(M_PI - phi0);
            particle[1].p.x = p*s/(r0*r0) * ( particle[1].r.y - particle[1].r.x * sqrt( pow( r0*p0/(p*s) ,2) - 1 ) );
            particle[1].p.y = sqrt(p0*p0 - particle[1].p.x * particle[1].p.x);
            //old pure numerical version. This also works now because of theta0 (asymptote approximation)
              /*particle[1].r.x = -r0;
                particle[1].r.y = s;
                particle[1].p.x = p;
                particle[1].p.y = 0; */
            //particle 0 is static scattering center, won't move (!!!)
            particle[0].m   = m_p;
            particle[0].q   = q1;

            double initial_angle = acos( particle[1].p * Vec(1,0,0) / particle[1].p.norm() ) *180/M_PI;
            printf( "phi0:%E, p0:%E, p_inf:%E, <(p0,x):%E\n", (M_PI-phi0)/M_PI*180, p0,p,initial_angle );
            printf("     ( %e )         ( %e )\n  ", particle[1].r.x, particle[1].p.x );
            printf("r0 = ( %e )  , p0 = ( %e )\n  ", particle[1].r.y, particle[1].p.y );
            printf("     ( %e )         ( %e )\n\n", particle[1].r.z, particle[1].p.z );

            E = CalcTotalEnergy(2,particle);
            L = CalcAngularMomentum(2,particle);
            printf( "E:%E, L:%E\n", E,L );

            //just a rough guess for the endtime to calculate dt. Analytically calculatable?!
            double E    = CalcTotalEnergy( 2, particle );
            const double t_end   = 2* fabs(particle[1].r.x) / sqrt( 2*E/particle[1].m );
            printf("tend:%e\n" , t_end);
            double steps = pow(10,i);
            const double dt = t_end / steps;

            struct return_data data = simulate_scattering( 2, particle, dt, "traj.dat", steps/200, steps );
            
            double dev_theta = (data.theta_max - data.theta_max_calc) / data.theta_max_calc;
            fprintf( stat, "%e\t%e\t%e\t%e\\n", particle[1].r.x, dt, steps, dev_theta );
            printf( "i:%i, j:%i => x:%e, dt:%e => devTheta:%e\n", i,j, particle[1].r.x, dt, dev_theta );
        }
    }

    fclose(stat);
    return 0;
}
