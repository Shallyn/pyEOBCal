/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "dyEvolution.h"
#include "myOptparser.h"
#include "myIO.h"

#include <time.h>

#define DEFAULT_ampO 0
#define DEFAULT_phiRef 0.0
#define DEFAULT_fRef 0.0
#define DEFAULT_m1 10.0
#define DEFAULT_m2 1.4
#define DEFAULT_f_min 40.0
#define DEFAULT_f_max 0.0
#define DEFAULT_distance 100.0
#define DEFAULT_s1x 0.0
#define DEFAULT_s1y 0.0
#define DEFAULT_s1z 0.0
#define DEFAULT_s2x 0.0
#define DEFAULT_s2y 0.0
#define DEFAULT_s2z 0.0
#define DEFAULT_eccentricity 0.0
#define DEFAULT_srate 16384
#define DEFAULT_eps 1e-5
#define DEFAULT_max_iterstep 100

typedef struct tagNQCGSParams {
    REAL8 deltaT;             /**< sampling interval */
    REAL8 m1;                 /**< mass of companion 1 */
    REAL8 m2;                 /**< mass of companion 2 */
    REAL8 f_min;              /**< start frequency */
    REAL8 e0;                 /**< eccentricity at start frequency */
    REAL8 s1x;                /**< (x,y,z) components of spin of m1 body */
    REAL8 s1y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s1z;                /**< dimensionless spin, Kerr bound: |s1| <= 1 */
    REAL8 s2x;                /**< (x,y,z) component ofs spin of m2 body */
    REAL8 s2y;                /**< z-axis along line of sight, L in x-z plane */
    REAL8 s2z;                /**< dimensionless spin, Kerr bound: |s2| <= 1 */
    FILE *file;
    REAL8 eps;
    INT maxiterstep;
    gboolean get_waveform;
} NQCPARAMS;

INT usage(const CHAR *program)
{
    INT a,c;
    print_err("usage: %s [options]\n", program);
    print_err("\t-h, --help\tprint this message and exit\n");
    print_err("\t-v, --verbose\tverbose output\n");
    print_err("\t-W --get-waveform\tWill generate waveform.\n");
    print_err("\t-e ECC, --eccentricity=ECC\n\t\torbital eccentricity at f_min [%g]\n", DEFAULT_eccentricity);
    print_err("\t-R SRATE, --sample-rate=SRATE\n\t\tsample rate in Hertz [%g]\n", DEFAULT_srate);
    print_err("\t-M M1, --m1=M1\n\t\tmass or primary in solar masses [%g]\n", DEFAULT_m1);
    print_err("\t-m M2, --m2=M2\n\t\tmass or secondary in solar masses [%g]\n", DEFAULT_m2);
    print_err( "\t-X S1X, --spin1x=S1X            \n\t\tx-component of dimensionless spin of primary [%g]\n", DEFAULT_s1x);
    print_err("\t-Y S1Y, --spin1y=S1Y            \n\t\ty-component of dimensionless spin of primary [%g]\n", DEFAULT_s1y);
    print_err("\t-Z S1Z, --spin1z=S1Z            \n\t\tz-component of dimensionless spin of primary [%g]\n", DEFAULT_s1z);
    print_err("\t-x S2X, --spin2x=S2X            \n\t\tx-component of dimensionless spin of secondary [%g]\n", DEFAULT_s2x);
    print_err("\t-y S2Y, --spin2y=S2Y            \n\t\ty-component of dimensionless spin of secondary [%g]\n", DEFAULT_s2y);
    print_err("\t-z S2Z, --spin2z=S2Z            \n\t\tz-component of dimensionless spin of secondary [%g]\n", DEFAULT_s2z);
    print_err("\t-f FMIN, --f-min=FMIN           \n\t\tfrequency to start waveform in Hertz [%g]\n", DEFAULT_f_min);
    print_err("\t-k KK, --KK=KK                 \n\t\tAdjustable Coefficient K.\n");
    print_err("\t-S dSS --dSS=dSS               \n\t\tAdjustable Coefficient dSS.\n");
    print_err("\t-O dSO, --dSO=dSO                 \n\t\tAdjustable Coefficient dSO.\n");
    print_err("\t-T dtPeak --dtPeak=dtPeak               \n\t\tAdjustable Coefficient dtPeak.\n");
    print_err("\t-F FILESXS --file-SXS=FILESXS               \n\t\tInput SXS waveform file.\n");
    print_err("\t-E EPS --eps=EPS                       \n\t\tError allowed for iteration [%g]\n", DEFAULT_eps);
    print_err("\t-I MAXITERSTEP --max-iterstep=MAXITERSTEP                  \n\t\tMaximum iteration step [%g]\n", DEFAULT_max_iterstep);
    return CEV_SUCCESS;
}

NQCPARAMS parseargs(INT argc, CHAR **argv, AdjParams *adjParams)
{
    NQCPARAMS p;
    p.m1 = DEFAULT_m1;
    p.m2 = DEFAULT_m2;
    p.f_min = DEFAULT_f_min;
    p.s1x = DEFAULT_s1x;
    p.s1y = DEFAULT_s1y;
    p.s1z = DEFAULT_s1z;
    p.s2x = DEFAULT_s2x;
    p.s2y = DEFAULT_s2y;
    p.s2z = DEFAULT_s2z;
    p.deltaT = 1./DEFAULT_srate;
    p.e0 = DEFAULT_eccentricity;
    p.maxiterstep = DEFAULT_max_iterstep;
    p.eps = DEFAULT_eps;
    p.file = NULL;
    p.get_waveform = FALSE;
    gboolean default_adj = TRUE;
    extern CHAR *EXT_optarg;
    extern INT EXT_optind;

    
    OPTION long_options[] = {
        {"help", opt_no_argument, 0, 'h'},
        {"verbose", opt_no_argument, 0, 'v'},
        {"get-waveform", opt_no_argument, 0, 'W'},
        {"eccentricity", opt_required_argument, 0, 'e'},
        {"sample-rate", opt_required_argument, 0, 'R'},
        {"m1", opt_required_argument, 0, 'M'},
        {"m2", opt_required_argument, 0, 'm'},
        {"spin1x", opt_required_argument, 0, 'X'},
        {"spin1y", opt_required_argument, 0, 'Y'},
        {"spin1z", opt_required_argument, 0, 'Z'},
        {"spin2x", opt_required_argument, 0, 'x'},
        {"spin2y", opt_required_argument, 0, 'y'},
        {"spin2z", opt_required_argument, 0, 'z'},
        {"f-min", opt_required_argument, 0, 'f'},
        {"KK", opt_required_argument , 0 ,'K'},
        {"dSS", opt_required_argument , 0 ,'S'},
        {"dSO", opt_required_argument , 0,'O'},
        {"dtPeak", opt_required_argument ,0 ,'T'},
        {"file-SXS", opt_required_argument, 0, 'F'},
        {"eps", opt_required_argument, 0, 'E'},
        {"max-iterstep", opt_required_argument, 0, 'I'},
        {0, 0, 0, 0}
    };
    CHAR args[] =
    "h:v:W:e:R:M:m:X:Y:Z:x:y:z:f:K:S:O:T:F:E:I";
    while (1)
    {
        INT option_index = 0;
        INT c;
        c = getopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
            case 0:
                if (long_options[option_index].flag)
                    break;
                else
                {
                    print_err("error parsing option %s with argument %s\n", long_options[option_index].name, EXT_optarg);
                    exit(1);
                }
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'W':
                p.get_waveform = TRUE;
                break;
            case 'e':
                p.e0 = atof(EXT_optarg);
                break;
            case 'R':
                p.deltaT = 1./atof(EXT_optarg);
                break;
            case 'M':
                p.m1 = atof(EXT_optarg);
                break;
            case 'm':
                p.m2 = atof(EXT_optarg);
                break;
            case 'X':
                p.s1x = atof(EXT_optarg);
                break;
            case 'Y':      /* spin1y */
                p.s1y = atof(EXT_optarg);
                break;
            case 'Z':      /* spin1z */
                p.s1z = atof(EXT_optarg);
                break;
            case 'x':      /* spin2x */
                p.s2x = atof(EXT_optarg);
                break;
            case 'y':      /* spin2y */
                p.s2y = atof(EXT_optarg);
                break;
            case 'z':      /* spin2z */
                p.s2z = atof(EXT_optarg);
                break;
            case 'f':
                p.f_min = atof(EXT_optarg);
                break;
            case 'K':
                adjParams->KK = atof(EXT_optarg);
                default_adj = FALSE;
                break;
            case 'S':
                adjParams->dSS = atof(EXT_optarg);
                break;
            case 'O':
                adjParams->dSO = atof(EXT_optarg);
                default_adj = FALSE;
                break;
            case 'T':
                adjParams->dtPeak = atof(EXT_optarg);
                default_adj = FALSE;
                break;
            case 'F':
                p.file = fopen(EXT_optarg, "r");
                if(!p.file)
                {
                    print_warning("Could not open SXS file %s\n", EXT_optarg);
                    exit(1);
                }
                break;
            case 'E':
                p.eps = atof(EXT_optarg);
                break;
            case 'I':
                p.maxiterstep = atoi(EXT_optarg);
                break;
            default:
                print_err("unknown error while parsing options\n");
                exit(1);
        }
    }
    if (EXT_optind < argc)
    {
        print_err("extraneous command line arguments:\n");
        while (EXT_optind < argc)
            print_err("%s\n", argv[EXT_optind++]);
        exit(1);
    }
    if (default_adj)
    {
        applyDefaultAdjustableParameters(adjParams, p.m1, p.m2, p.s1z, p.s2z);
    }
    if(!p.file)
    {
        print_warning("You have not specified SXS file, will use fit one.\n");
    }
    return p;
}

void printNQC(EOBNonQCCoeffs *nqcParams);
void printhLMAll(COMPLEX16TimeSeries *hLM);

INT main(INT argc, CHAR **argv)
{
    NQCPARAMS p;
    INT status;
    AdjParams adjParams;
    CtrlParams ctrlParams;
    p = parseargs(argc, argv, &adjParams);
    EOBNonQCCoeffs nqcParams;
#if DEBUG
print_debug("CMD: --m1 %f --m2 %f --f-min %f --e0 %f --spin1z %f --spin2z %f --KK %f --dSS %f --dSO %f --dtPeak %f\n", 
    p.m1, p.m2, p.f_min, p.e0, p.s1z, p.s2z, adjParams.KK, adjParams.dSS, adjParams.dSO, adjParams.dtPeak);
#endif
    REAL8Vector *time = NULL;
    REAL8Vector *hSXSreal = NULL;
    REAL8Vector *hSXSimag = NULL;
    if(p.file)
    {
        status = read_waveform(&time, &hSXSreal, &hSXSimag, p.file);
        if(status != CEV_SUCCESS)
        {
            print_warning("Failed to load SXS file, quit.\n");
            if(time)
                DestroyREAL8Vector(time);
            if(hSXSreal)
                DestroyREAL8Vector(hSXSreal);
            if(hSXSimag)
                DestroyREAL8Vector(hSXSimag);
            return -1;
        }
    }
    COMPLEX16TimeSeries *hLMAll = NULL;
    status = IterateNQCCorrectionCoeffs(p.m1, p.m2, p.f_min, p.e0, p.deltaT, 
        p.s1z, p.s2z, time, hSXSreal, hSXSimag, &nqcParams, &hLMAll, &adjParams, &ctrlParams, p.eps, p.maxiterstep, p.get_waveform);
    if(status != CEV_SUCCESS)
    {
        print_warning("Failed! return code = %d", status);
        if(time)
            DestroyREAL8Vector(time);
        if(hSXSreal)
            DestroyREAL8Vector(hSXSreal);
        if(hSXSimag)
            DestroyREAL8Vector(hSXSimag);
        if(hLMAll)
            DestroyCOMPLEX16TimeSeries(hLMAll);
        return -1;
    }

    if(hLMAll)
    {
        printhLMAll(hLMAll);
    }
    else
    {
        printNQC(&nqcParams);
    }

    DestroyREAL8Vector(time);
    DestroyREAL8Vector(hSXSreal);
    DestroyREAL8Vector(hSXSimag);
    if(hLMAll)
        DestroyCOMPLEX16TimeSeries(hLMAll);
    return 0;
}

void printNQC(EOBNonQCCoeffs *nqcParams)
{
    print_out("#a1 #a2 #a3 #a3S #a4 #a5 #b1 #b2 #b3 #b4\n");
    print_out("%le %le %le %le %le %le %le %le %le\n", 
                nqcParams->a1, nqcParams->a2, nqcParams->a3,
                nqcParams->a3S, nqcParams->a4, nqcParams->a5,
                nqcParams->b1, nqcParams->b2, nqcParams->b3,
                nqcParams->b4);
    return;
}

void printhLMAll(COMPLEX16TimeSeries *hLM)
{
    UINT length, i;
    REAL8 t0, dt;
    t0 = hLM->epoch;
    dt = hLM->deltaT;
    length = hLM->data->length;
    print_out("#time #hreal #himag\n", t0 + i*dt, creal(hLM->data->data[i]),cimag(hLM->data->data[i]));

    for (i=0;i<length;i++)
    {
        print_out("%.16e %.16e %.16e\n", t0 + i*dt, creal(hLM->data->data[i]),cimag(hLM->data->data[i]));
    }

    return;
}
