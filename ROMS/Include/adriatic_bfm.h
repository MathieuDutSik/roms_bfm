/*
**  Options for 2km ADRIATIC RUN
*/
#define DEBUG
#define UV_ADV 		/* turn ON or OFF advection terms */
#define UV_COR 		/* turn ON or OFF Coriolis term */
#undef DIAGNOSTICS_UV 	/* define if writing out momentum diagnostics */
#define CURVGRID  	/* define if using  curvilinear coordinate grid*/
#define DJ_GRADPS 	/* Splines density  Jacobian (Shchepetkin, 2000) */
/* define one advection scheme here*/
#undef  TS_U3HADVECTION/* define if 3rd-order upstream horiz. advection */
#define TS_MPDATA	/* positive defi for TS */
#undef T_PASSIVE	/*using pasive tracer dye_01*/
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define RI_SPLINES
#define TS_DIF2 	/* turn ON or OFF Laplacian horizontal mixing */
#define NONLIN_EOS 	/* define if using nonlinear equation of state */
#define SALINITY 	/* define if using salinity */
#define MASKING 	/* define if there is land in the domain */
#undef TS_FIXED	/* this is diagnostic run!! */
#define SOLVE3D 	/* define if solving 3D primitive equations */

#define STATIONS 	/* define if writing out station data */
#define FLOATS 		/* define drifters */

#undef READ_WATER
#undef WRITE_WATER
#undef OUT_DOUBLE
#define PERFECT_RESTART
#undef ANA_INITIAL
#undef RST_SINGLE

#define BFM_COUPLING

#undef ANA_TAIR
#undef ANA_PAIR
#undef ANA_HUMIDITY

#undef DIAGNOSTICS_UV
#define STOKES_DRIFT_OLD_VERSION

#define WIND_STABILITY_CORRECTION

#undef SWAN_COUPLING
#ifdef SWAN_COUPLING
# define NEARSHORE_MELLOR
# undef NEARSHORE_MELLOR05
# define MCT_LIB
#endif
#undef WWM_COUPLING
#undef WWM_MPI
#define ROMS_WWM_PGMCL_COUPLING
#undef DUMMY_COUPLING
#ifdef WWM_COUPLING
# define NEARSHORE_ARDHUIN
# undef NEARSHORE_MELLOR05
# undef NEARSHORE_LONGUETHIGGINS
# ifdef NEARSHORE_ARDHUIN
#  define WAVE_ADVECTION
#  define WAVE_ADVECTION_TRACER
#  define WAVE_ADVECTION_MOMENTUM
#  define WAVE_ADVECTION_TURBULENCE
#  define APPLY_VORTEX
#  define APPLY_PRESSURE
#  define WAVE_BOUNDARY
#  define FIRST_ORDER_ARDHUIN
#  define WAVE_COR
#  define STOKES_CORR_FLOAT
#  define STOKES_DRIFT_USING_INTEGRAL
#  define USE_STRESS_FROM_WAVE
#  define USE_WAVE_PRESSURE_INTEGRAL
#  define GET_CD_UFRIC_ROUGHNESS
#  define ZOS_HSIG       /* have to have Hsig */
#  define CRAIG_BANNER
# endif
# ifdef NEARSHORE_MELLOR05
#  define WAVE_ADVECTION
#  define WAVE_ADVECTION_TRACER
#  define WAVE_ADVECTION_MOMENTUM
#  undef WAVE_BOUNDARY
#  undef FIRST_ORDER_ARDHUIN
#  define WAVE_COR
#  define STOKES_CORR_FLOAT
#  define STOKES_DRIFT_USING_INTEGRAL
#  define USE_STRESS_FROM_WAVE
#  define USE_WAVE_PRESSURE_INTEGRAL
#  define GET_CD_UFRIC_ROUGHNESS
#  define ZOS_HSIG       /* have to have Hsig */
#  define CRAIG_BANNER
# endif
#endif

#undef COARE_TAYLOR_YELLAND
#undef COARE_OOST




#undef ASSIMILATION_SST
#undef NUDGING_SST

#define  NOSEDBBL
#ifdef NOSEDBBL
/* BBL run with no water column sediment */
# undef SEDIMENT
# undef SUSPLOAD
# undef ANA_SEDIMENT
# undef ANA_WWAVE
# undef SWAN
# undef RIVER_SEDIMENT
#else
/* BBL run with full sediment run */
# define SEDIMENT
# define SUSPLOAD
# undef  ANA_SEDIMENT
# undef  ANA_WWAVE
# define SWAN
# define RIVER_SEDIMENT
#endif

/* define only one of the four following */
#define UV_LOGDRAG /* turn ON or OFF logarithmic bottom friction */
#undef MB_BBL
#undef SG_BBL
#undef SSW_BBL

#ifdef SG_BBL
#define SG_CALC_ZNOT
#undef  SG_LOGINT
#endif

#ifdef MB_BBL
#define MB_CALC_ZNOT
#undef  MB_Z0BIO
#undef  MB_Z0BL
#undef  MB_Z0RIP
#endif

#ifdef SSW_BBL
#define SSW_CALC_ZNOT
#undef  SSW_LOGINT
#endif

/* define one vertical mixing scheme here*/
#undef  LMD_MIXING
#undef  MY25_MIXING
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif
/* waves effects at the water column */
#ifndef ZOS_HSIG
# define CHARNOK
#endif
#undef TKE_WAVEDISS   /* have to have Wave_dissipation */

/* define atmo model WRF, ALADIN COAMPS or LAMI */


#define SPONGE
#define MIX_S_UV
#define TS_DIF2
#define UV_VIS2
#define TS_DIF2         /* turn ON or OFF Laplacian horizontal mixing */
#define MIX_GEO_TS 

#define ALADIN

#ifdef WRF
# define BULK_FLUXES
# define EMINUSP /* compute evaporation and effect on salinity */
# define COOL_SKIN
# define SOLAR_SOURCE /* define solar radiation source term */
# define ANA_SSFLUX /* analytical surface salinity flux */
# define ANA_BSFLUX /* analytical bottom salinity flux */
# define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
# define ANA_BTFLUX /* analytical bottom temperature flux */
# define ANA_SPFLUX /* analytical surface passive tracers fluxes */
# define LONGWAVE_OUT /* Compute net longwave radiation internally from Tair,Qair,Cloud,SST*/
# undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
# define ANA_RAIN /* analytical rain fall rate */
#elif defined ALADIN
#   define BULK_FLUXES
#   define EMINUSP /* compute evaporation and effect on salinity */
#   define COOL_SKIN
#   define SOLAR_SOURCE /* define solar radiation source term */
#   undef ANA_SSFLUX /* analytical surface salinity flux */
#   define ANA_BSFLUX /* analytical bottom salinity flux */
#   define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
#   define ANA_BTFLUX /* analytical bottom temperature flux */
#   define ANA_SPFLUX /* analytical surface passive tracers fluxes */
#   define LONGWAVE /* Compute net longwave radiation internally from Tair,Qair,Cloud,SST*/
#   undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
#   define ANA_RAIN /* analytical rain fall rate */
#elif defined LAMI_MODIF
# define BULK_FLUXES
# undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
# undef ALBEDO          /* use albedo equation for shortwave radiation */
# define ANA_SSFLUX /* analytical surface salinity flux */
# define ANA_BSFLUX /* analytical bottom salinity flux */
# define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
# define ANA_BTFLUX /* analytical bottom temperature flux */
# define ANA_SPFLUX /* analytical surface passive tracers fluxes */
# undef LONGWAVE /* Compute net longwave radiation internally from Tair,Qair,Cloud,SST*/
# define SOLAR_SOURCE /* define solar radiation source term */
# define ANA_RAIN /* analytical rain fall rate */
# undef  DIURNAL_SRFLUX /* impose shortwave radiation local diurnal cycle SWR is 1/per day*/
#elif defined LAMI
# define BULK_FLUXES
# undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
# undef ALBEDO          /* use albedo equation for shortwave radiation */
# define ANA_SSFLUX /* analytical surface salinity flux */
# define ANA_BSFLUX /* analytical bottom salinity flux */
# define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
# define ANA_BTFLUX /* analytical bottom temperature flux */
# define ANA_SPFLUX /* analytical surface passive tracers fluxes */
# define LONGWAVE /* Compute net longwave radiation internally from Tair,Qair,Cloud,SST*/
# undef SOLAR_SOURCE /* define solar radiation source term */
# define ANA_RAIN /* analytical rain fall rate */
/* define DIURNAL_SRFLUX for LAMI and ALADIN driven run*/
# define  DIURNAL_SRFLUX /* impose shortwave radiation local diurnal cycle SWR is 1/per day*/
#elif defined COAMPS
# define BULK_FLUXES
# undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
# undef ALBEDO          /* use albedo equation for shortwave radiation */
# define ANA_CLOUD
# define ANA_SSFLUX /* analytical surface salinity flux */
# define ANA_BSFLUX /* analytical bottom salinity flux */
# define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
# define ANA_BTFLUX /* analytical bottom temperature flux */
# define ANA_SPFLUX /* analytical surface passive tracers fluxes */
# define LONGWAVE /* Compute net longwave radiation internally from Tair,Qair,Cloud,SST*/
# undef SOLAR_SOURCE /* define solar radiation source term */
# define ANA_RAIN /* analytical rain fall rate */
#else
# define ANA_SMFLUX
# define ANA_SSFLUX /* analytical surface salinity flux */
# define ANA_BSFLUX /* analytical bottom salinity flux */
# define ANA_BPFLUX /* analytical bottom passive tracers fluxes */
# define ANA_BTFLUX /* analytical bottom temperature flux */
# define ANA_SPFLUX /* analytical surface passive tracers fluxes */
# define ANA_STFLUX 
#endif


/* define the boundary conditions */

#define SSH_TIDES
#define UV_TIDES

#define ADD_FSOBC
#define ADD_M2OBC

#define RADIATION_2D
