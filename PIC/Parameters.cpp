//componentsConfig.param
const uint16_t SIMDIM                 = 2;
const uint16_t NUMBER_OF_NEIGHBORS    = 0;  //Neighboring Cells to be included for periodic Force. 0 means no periodic summation
const uint32_t NUMBER_OF_STEPS        = 1e7;

//GasConfig.param
const double GAS_DENSITY_SI           = 1.e30; // 1/m^3

//physicalConstants.param
const double ELEMENTARY_CHARGE_SI     = 1.602176565e-19;  // C
const double PROTON_MASS_SI           = 1.6726231e-27;    // kg
const double SPEED_OF_LIGHT_SI        = 2.99792458e8;     // m/s
const double MUE0_SI                  = M_PI * 4.e-7;     // N/A^2 = kg*m/C^2
const double EPS0_SI                  = 1.0/(MUE0_SI*SPEED_OF_LIGHT_SI*SPEED_OF_LIGHT_SI);    // C^2/J*m, 8.854187817e-12

//particleConfig.param
const double ELECTRON_TEMPERATURE_keV = 0.1;
const double ION_TEMPERATURE_keV      = 0.1;
const double ELECTRON_MASS_SI         = 9.109382e-31;         // kg
const double ELECTRON_CHARGE_SI       =-ELEMENTARY_CHARGE_SI; // C
const double ION_MASS_SI              = PROTON_MASS_SI;       // kg
const double ION_CHARGE_SI            = ELEMENTARY_CHARGE_SI; // C
const uint32_t NUMBER_OF_PARTICLES_PER_CELL = 25;   // NUM instead of NUMBER_OF in picongpu is also inconsistent, and there are other longer names, soo ...
const uint16_t DEFAULT_PARTICLE_SHAPE = 1;  //00:point-point, 01:ball-ball (CIC radialsymmetric equivalent), 99:sphere-sphere    

//GridConfig.param
const double DELTA_T_SI               = 1.5e-20;
const uint32_t NUMBER_OF_CELLS_X      = 2;
const uint32_t NUMBER_OF_CELLS_Y      = 2;
const uint32_t NUMBER_OF_CELLS_Z      = 1;
const double CELL_SIZE_SI             = pow( double(NUMBER_OF_PARTICLES_PER_CELL)/GAS_DENSITY_SI, 1./3. );
const double CELL_SIZE_X_SI           = CELL_SIZE_SI;   // (!!!) picongpu naming with width, height and depth seems to be too random => could lead to mixups
const double CELL_SIZE_Y_SI           = CELL_SIZE_SI;
const double CELL_SIZE_Z_SI           = CELL_SIZE_SI;
const uint16_t BOUNDARY_CONDITION     = 0;  //0:periodic, 1:reflecting, 2:adhering

//================================== Units ===================================//
const double UNITCONV_keV_to_Joule    = 1.e3 * ELEMENTARY_CHARGE_SI;
const double UNITCONV_Joule_to_keV    = 1.0 / UNITCONV_keV_to_Joule;

const double UNIT_LENGTH              = CELL_SIZE_SI;
const double UNIT_MASS                = ELECTRON_MASS_SI;
const double UNIT_CHARGE              = ELEMENTARY_CHARGE_SI;
const double UNIT_TIME                = DELTA_T_SI;

//derived units
const double UNIT_SPEED               = UNIT_LENGTH / UNIT_TIME;
const double UNIT_ENERGY              = UNIT_MASS * UNIT_LENGTH*UNIT_LENGTH / (UNIT_TIME*UNIT_TIME); // in general not eV !
const double UNIT_MOMENTUM            = UNIT_MASS * UNIT_SPEED;
const double UNIT_ANGULAR_MOMENTUM    = UNIT_LENGTH * UNIT_MASS * UNIT_SPEED;
const double UNIT_FORCE               = UNIT_MASS * UNIT_LENGTH / (UNIT_TIME*UNIT_TIME);

//physicalConstants.unitless   //In Picongpu the variable names change to Q_EL, M_EL and so on, which I find inconsistent -.- (!!!)
const double ELEMENTARY_CHARGE        = ELEMENTARY_CHARGE_SI / UNIT_CHARGE;  // C
const double SPEED_OF_LIGHT           = SPEED_OF_LIGHT_SI / UNIT_SPEED;  // m/s
const double ELECTRON_CHARGE          = ELECTRON_CHARGE_SI / UNIT_CHARGE;
const double ELECTRON_MASS            = ELECTRON_MASS_SI / UNIT_MASS;
const double MUE0                     = MUE0_SI / (UNIT_MASS * UNIT_LENGTH / (UNIT_CHARGE*UNIT_CHARGE));    //! magnetic constant must be double 3.92907e-39
const double EPS0                     = 1. / ( MUE0*SPEED_OF_LIGHT*SPEED_OF_LIGHT );  //! electric constant must be double 2.54513e+38

//particleConfig.unitless
const double ELECTRON_TEMPERATURE_SI  = ELECTRON_TEMPERATURE_keV * UNITCONV_keV_to_Joule; // J
const double ELECTRON_TEMPERATURE     = ELECTRON_TEMPERATURE_SI / UNIT_ENERGY;
const double ION_TEMPERATURE_SI       = ION_TEMPERATURE_keV * UNITCONV_keV_to_Joule;      // J
const double ION_TEMPERATURE          = ION_TEMPERATURE_SI / UNIT_ENERGY;
const double ION_CHARGE               = ION_CHARGE_SI / UNIT_CHARGE;
const double ION_MASS                 = ION_MASS_SI / UNIT_MASS;
const uint32_t NUMBER_OF_PARTICLES    = NUMBER_OF_PARTICLES_PER_CELL * NUMBER_OF_CELLS_X * NUMBER_OF_CELLS_Y * NUMBER_OF_CELLS_Z;

//GridConfig.unitless
const double DELTA_T                  = DELTA_T_SI / UNIT_TIME;
const double CELL_SIZE_X              = CELL_SIZE_X_SI / UNIT_LENGTH;
const double CELL_SIZE_Y              = CELL_SIZE_Y_SI / UNIT_LENGTH;
const double CELL_SIZE_Z              = CELL_SIZE_Z_SI / UNIT_LENGTH;
const uint32_t CELL_SIZE[3]           = { CELL_SIZE_X, CELL_SIZE_Y, CELL_SIZE_Z };
const uint32_t NUMBER_OF_CELLS[3]     = { NUMBER_OF_CELLS_X, NUMBER_OF_CELLS_Y, NUMBER_OF_CELLS_Z };
