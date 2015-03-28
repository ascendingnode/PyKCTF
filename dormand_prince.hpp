#ifndef __DORMAND_PRINCE__
#define __DORMAND_PRINCE__

#include <cstdio>
#include <cmath>
#include <vector>

// Base class; classes to be integrated should inherit this
class Dormand_Prince_Base { public:
    double x;
    std::vector<double> y,dy;
    Dormand_Prince_Base() {}
    void dt(const double &xx,const std::vector<double> &yy,std::vector<double> &dyy) const {}
    bool is_sane() const { return true; }
    void renorm() {}
};

// The integrator class, template argument is the tye of class to be integrated
template<class A> class Dormand_Prince { public:

    bool first;
    double C2,C3,C4,C5, H;
    double A21,A31,A32,A41,A42,A43,A51,A52,A53,A54;
    double A61,A62,A63,A64,A65,A71,A73,A74,A75,A76;
    double E1,E3,E4,E5,E6,E7;
    double ATOLI,RTOLI,UROUND,SAFE,FAC1,FAC2,BETA;
    unsigned nvar,NFCN,NSTEP,NACCPT,NREJCT,NMAX;
    std::vector<double> K1,K2,K3,K4,K5,K6,Y1,YSTI;

    Dormand_Prince() {
        ATOLI = 1e-10;
        RTOLI = 1e-10;
        init();
    }

    Dormand_Prince(double tol) {
        ATOLI = RTOLI = tol;
        init();
    }

    // Initialize common parameters
    void init() {

        first = true;
        H = 0;

        // Coefficients
        C2 = 0.2;
        C3 = 0.3;
        C4 = 0.8;
        C5 = 8.0/9.0;
        A21 = 0.2;
        A31 = 3.0/40.0;
        A32 = 9.0/40.0;
        A41 = 44.0/45.0;
        A42 = -56.0/15.0;
        A43 = 32.0/9.0;
        A51 = 19372.0/6561.0;
        A52 = -25360.0/2187.0;
        A53 = 64448.0/6561.0;
        A54 = -212.0/729.0;
        A61 = 9017.0/3168.0;
        A62 = -355.0/33.0;
        A63 = 46732.0/5247.0;
        A64 = 49.0/176.0;
        A65 = -5103.0/18656.0;
        A71 = 35.0/384.0;
        A73 = 500.0/1113.0;
        A74 = 125.0/192.0;
        A75 = -2187.0/6784.0;
        A76 = 11.0/84.0;
        E1 = 71.0/57600.0;
        E3 = -71.0/16695.0;
        E4 = 71.0/1920.0;
        E5 = -17253.0/339200.0;
        E6 = 22.0/525.0;
        E7 = -1.0/40.0;

        // SETTING THE PARAMETERS 
        nvar = 0;
        NFCN = 0;
        NSTEP = 0;
        NACCPT = 0;
        NREJCT = 0;
        // -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- 
        NMAX=100000;
        // -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
        UROUND = 1e-35;
        // -------  SAFETY FACTOR -------------
        SAFE = 0.9;
        // -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
        FAC1=0.2;
        FAC2=10.0;
        // --------- BETA FOR STEP CONTROL STABILIZATION -----------
        BETA = 0.04;
    }

    // Allocate the storage arrays
    void allocate_arrays(unsigned n) {
        K1.resize(n); K2.resize(n); K3.resize(n);
        K4.resize(n); K5.resize(n); K6.resize(n);
        Y1.resize(n); YSTI.resize(n);
    }

    //  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
    void integrate(A &sys,double XEND) {
    
        if(nvar!=sys.y.size()) {
            nvar = sys.y.size();
            allocate_arrays(nvar);
        }

        double FACOLD=1.e-4;
        double EXPO1=0.20-BETA*0.750;
        double FACC1=1.0/FAC1;
        double FACC2=1.0/FAC2;
        double POSNEG = (XEND-sys.x)/fabs(XEND-sys.x);
        
        // --- INITIAL PREPARATIONS   
        bool LAST = false;
        sys.dt(sys.x,sys.y,K1);  

        double HMAX = fabs(XEND-sys.x);
        if(first) {
            H = initial_H(sys,HMAX,POSNEG);
            first = false;
        }
        NFCN += 2;
        bool REJECT = false;
        unsigned IRTRN = 0;

        // --- BASIC INTEGRATION STEP  
        while(1) {

            if(!sys.is_sane()) return;
        
            //IF (NSTEP.GT.NMAX) GOTO 78
            H = std::max(fabs(H),fabs(sys.x)*UROUND*10.)*POSNEG;
            if ((sys.x+1.01*H-XEND)*POSNEG > 0.0) {
                H = XEND-sys.x;
                LAST = true;
            }
            NSTEP = NSTEP+1;
            
            // --- THE FIRST 6 STAGES
            if(IRTRN>=2) sys.dt(sys.x,sys.y,K1);
            for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*A21*K1[i];
            sys.dt(sys.x+C2*H,Y1,K2);
            for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*(A31*K1[i] + A32*K2[i]);
            sys.dt(sys.x+C3*H,Y1,K3);
            for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*(A41*K1[i] + A42*K2[i] + A43*K3[i]);
            sys.dt(sys.x+C4*H,Y1,K4);
            for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*(A51*K1[i] + A52*K2[i] + A53*K3[i] + A54*K4[i]);
            sys.dt(sys.x+C5*H,Y1,K5);
            for(int i=0;i<nvar;i++) YSTI[i] = sys.y[i] + H*(A61*K1[i] + A62*K2[i] + A63*K3[i] + A64*K4[i] + A65*K5[i]);
            double XPH = sys.x+H;
            sys.dt(XPH,YSTI,K6);
            for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*(A71*K1[i] + A73*K3[i] + A74*K4[i] + A75*K5[i] + A76*K6[i]);
            sys.dt(XPH,Y1,K2);
            for(int i=0;i<nvar;i++) K4[i] = (E1*K1[i] + E3*K3[i] + E4*K4[i] + E5*K5[i] + E6*K6[i] + E7*K2[i])*H;
            NFCN += 6;
            
            // --- ERROR ESTIMATION  
            double ERR=0.0;
            //printf("%g ",H);
            for(unsigned i=0;i<nvar;i++) {
                double SK = ATOLI + RTOLI*std::max(fabs(sys.y[i]),fabs(Y1[i]));
                double K4SK = K4[i]/SK;
                //printf("%g ",K4SK*K4SK);
                ERR += K4SK*K4SK;
            }
            //printf("\n");
            ERR = sqrt(ERR/(double)(nvar));

            // --- COMPUTATION OF HNEW
            double FAC11 = pow(ERR,EXPO1);
            
            // --- LUND-STABILIZATION
            double FAC = FAC11/pow(FACOLD,BETA);
    
            // --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
            FAC = std::max(FACC2,std::min(FACC1,FAC/SAFE));
            double HNEW = H/FAC;  
      
            // --- STEP IS ACCEPTED
            if (ERR <= 1.0) {
        
                FACOLD = std::max(ERR,1.0e-4);
                NACCPT = NACCPT+1;
        
                K1 = K2;
                sys.y = Y1;
                sys.x = XPH;
                sys.renorm();

                // ------- NORMAL EXIT
                if (LAST) {
                    H = HNEW;
                    return;
                }
                if(fabs(HNEW)>HMAX) HNEW = POSNEG*HMAX;
                if(REJECT) HNEW = POSNEG*std::min(fabs(HNEW),fabs(H));
                REJECT = false;

            } 
            // --- STEP IS REJECTED
            else {

                HNEW = H/std::min(FACC1,FAC11/SAFE);
                REJECT = true;  
                if(NACCPT>=1) NREJCT = NREJCT+1;
                LAST = false;
            }
            H = HNEW;
        }
    }

    // ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
    double initial_H(A &sys, double HMAX, double POSNEG) {
        
        // ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
        // ----   H = 0.01 * NORM (Y0) / NORM (F0)
        // ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
        // ---- COMPARED TO THE SOLUTION
        std::vector<double> SK(nvar);
        double DNF=0, DNY=0;
        for(int i=0;i<nvar;i++) {
            SK[i] = ATOLI + RTOLI*sys.y[i];
            DNF += (   K1[i]/SK[i])*(   K1[i]/SK[i]);
            DNY += (sys.y[i]/SK[i])*(sys.y[i]/SK[i]);
        }
        double H = (DNF<=1e-10 or DNY<=1.e-10) ?
            1.0e-6 :
            sqrt(DNY/DNF)*0.01;
        H = std::min(H,HMAX)*POSNEG; 
        
        // ---- PERFORM AN EXPLICIT EULER STEP
        //Y1 = Y + H*K1;
        for(int i=0;i<nvar;i++) Y1[i] = sys.y[i] + H*K1[i];
        sys.dt(sys.x+H,Y1,K2);
        
        // ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
        double DER2 = 0;
        for(int i=0;i<nvar;i++) {
            double diff = (K2[i]-K1[i])/SK[i];
            DER2 += diff*diff;
        }
        DER2 = sqrt(DER2)/H;

        
        // ---- STEP SIZE IS COMPUTED SUCH THAT
        // ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
        double DER12 = std::max(fabs(DER2),sqrt(DNF));
        double H1 = (DER12<=1e-15) ?
            std::max(1.0e-6,fabs(H)*1.0e-3) :
            pow(0.01/DER12, 1.0/5.0);
        return std::min(100.*fabs(H),std::min(H1,HMAX))*POSNEG;
    }
};

#endif
