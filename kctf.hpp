#ifndef __KCTF_HPP__
#define __KCTF_HPP__

#include "dormand_prince.hpp"

namespace kctf_const {
    const double Gc = 6.67384e-11;
    const double m_sun = 1.98892e30;
    const double AU = 149597870700.0;
    const double d2r = M_PI/180.;
    const double day = 24.*3600.;
    const double year = 365.25*day;
    const double c = 299792458.;
}

class Tidal_Body { public:
    double m,P,theta0,phi0,R,rg,k2,dt;
    Tidal_Body() {}
};

class KCTF : public Dormand_Prince_Base { public:

    double J20_const,J21_const;
    double m0,m1,m2, mu1,mu2,beta1,beta2, mindist,mine;
    double R0,R1, C0,C1, k20,k21, K0c,K1c, gammac,degrc;
    double J20c,J21c, alpha10c,alpha11c, alpha20c,alpha21c;
    std::vector<double> e2,G0;
    Dormand_Prince<KCTF> dp;

    KCTF() {}

    KCTF(const Tidal_Body &b0,const Tidal_Body &b1,double m20,double a1,double e1m,double I,double w,double a2,double e2m) {
        using namespace kctf_const;
        J20_const = J21_const = 0;
        m0 = b0.m; m1 = b1.m; m2 = m20;
        mu1 = Gc*(m0+m1);
        mu2 = Gc*(m0+m1+m2);
        beta1 = (m0*m1)/(m0+m1);
        beta2 = ((m0+m1)*m2)/(m0+m1+m2);
        R0 = b0.R; R1 = b1.R;
        C0 = b0.rg*m0*R0*R0;
        C1 = b1.rg*m1*R1*R1;
        k20 = b0.k2; k21 = b1.k2;
        K0c = b0.dt*3.*k20*Gc*m1*m1;
        K1c = b1.dt*3.*k21*Gc*m0*m0;
        J20c = (k20*R0*R0*R0)/(3.*Gc*m0);
        J21c = (k21*R1*R1*R1)/(3.*Gc*m1);
        alpha10c = 1.5*Gc*m0*m1*R0*R0;
        alpha11c = 1.5*Gc*m0*m1*R1*R1;
        alpha20c = 1.5*Gc*m2*m0*R0*R0;
        alpha21c = 1.5*Gc*m2*m1*R1*R1;
        gammac = 0.75*Gc*m2*beta1;
        degrc = (3.*mu1)/(c*c);

        double ele1[] = {a1*(1-e1m),e1m,I,0,w};
        double ele2[] = {a2*(1-e2m),e2m,0,0,0};
        std::vector<double> G1(3),e1(3), G2(3); 
        e2.resize(3); G0.resize(3);
        ele2Ge(ele1,m0   ,m1,G1,e1);
        ele2Ge(ele2,m0+m1,m2,G2,e2);
        double Cw0 = (C0*2.*M_PI)/b0.P;
        double Cw1 = (C1*2.*M_PI)/b1.P;
        std::vector<double> k1 = div(G1,norm(G1));
        std::vector<double> L0 = mult(Cw0,make_spin(k1,I,b0.theta0,b0.phi0));
        std::vector<double> L1 = mult(Cw1,make_spin(k1,I,b1.theta0,b1.phi0));
        y.resize(15); dy.resize(15);
        for(unsigned i=0;i<3;i++) {
            y[  i] = G1[i]; y[ 3+i] = G2[i]; y[6+i] = e1[i];
            y[9+i] = L0[i]; y[12+i] = L1[i];
            G0[i] = G1[i] + G2[i] + L0[i] + L1[i];
        }

        x = 0.;
        dt(x,y,dy);

        double roche0 = 1.26*R0*pow(m1/m0,1./3.);
        double roche1 = 1.26*R1*pow(m0/m1,1./3.);
        mindist = std::max(roche0,roche1);
        mine = 0;
    }

    static double dot(const std::vector<double> &u,const std::vector<double> &v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    }

    static double norm(const std::vector<double> &u) {
        return sqrt(dot(u,u));
    }

    static std::vector<double> cross(const std::vector<double> &u,const std::vector<double> &v) {
        std::vector<double> ret(3);
        ret[0] = u[1]*v[2] - u[2]*v[1];
        ret[1] = u[2]*v[0] - u[0]*v[2];
        ret[2] = u[0]*v[1] - u[1]*v[0];
        return ret;
    }

    static std::vector<double> add(const std::vector<double> &u,const std::vector<double> &v) {
        std::vector<double> ret(3);
        ret[0] = u[0]+v[0]; ret[1] = u[1]+v[1]; ret[2] = u[2]+v[2];
        return ret;
    }

    static std::vector<double> sub(const std::vector<double> &u,const std::vector<double> &v) {
        std::vector<double> ret(3);
        ret[0] = u[0]-v[0]; ret[1] = u[1]-v[1]; ret[2] = u[2]-v[2];
        return ret;
    }

    static std::vector<double> mult(double s,const std::vector<double> &u) {
        std::vector<double> ret(3);
        ret[0] = s*u[0]; ret[1] = s*u[1]; ret[2] = s*u[2];
        return ret;                                  
    }

    static std::vector<double> mult(const std::vector<double> &u,double s) {
        return mult(s,u);
    }

    static std::vector<double> div(const std::vector<double> &u,double s) {
        std::vector<double> ret(3);                        
        ret[0] = u[0]/s; ret[1] = u[1]/s; ret[2] = u[2]/s;
        return ret;                                      
    }

    void ele2Ge(double ele[5],double ma,double mb,std::vector<double> &Ga,std::vector<double> &ea) {
        using namespace kctf_const;
        double beta = (ma*mb)/(ma+mb);
        double mu = Gc*(ma+mb);
        double a = ele[0]/(1.-ele[1]);
        double e=ele[1], i=ele[2], O=ele[3], w=ele[4];
        double muh = sqrt(mu/(a*(1-e*e)));
        std::vector<double> r(3),v(3);
        double rm = ele[0];
        double si=sin(i), ci=cos(i);
        double sO=sin(O), cO=cos(O);
        double sw=sin(w), cw=cos(w);
        r[0] = rm*(cO*cw - sO*sw*ci);
        r[1] = rm*(sO*cw + cO*sw*ci);
        r[2] = rm*            sw*si;
        v[0] = -muh*(cO*(sw + e*sw) + sO*(cw + e*cw)*ci);
        v[1] = -muh*(sO*(sw + e*sw) - cO*(cw + e*cw)*ci);
        v[2] =  muh*                     (cw + e*cw)*si;
        Ga = mult(beta,cross(r,v));
        ea = add( div(cross(v,Ga),beta*mu), div(r,-norm(r)) );
    }

    std::vector<double> make_spin(const std::vector<double>&k,double I,double theta,double phi) {
        double ct = cos(theta), ce = cos(I)*ct + sin(I)*sin(theta)*cos(phi);
        double kix=k[0], kiy=k[1], kiz = k[2], six;
        double ce2=ce*ce, kix2=kix*kix, kiy2=kiy*kiy, kiz2=kiz*kiz;
        double rad1 = -ce2*kiz2 - (ce2 - 1.)*kix2 - (ce2 - 1.)*kiy2 + 2.*ce*ct*kiz - ct*ct;
        if(rad1>0)
            six = -(ce*kix*kiz - ct*kix + sqrt(rad1)*kiy)/(kix2 + kiy2);
        else {
            rad1 = -ce2*kiz2 - (ce2 - 1)*kix2 - (ce2 - 1)*kiy2 + 2*ce*ct*kiz - ct*ct;
            six = -(ce*kix*kiz - ct*kix - sqrt(rad1)*kiy)/(kix2 + kiy2); }
        double siy = sqrt(1.-six*six-ce2);
        std::vector<double> si(3); si[0]=six; si[1]=siy; si[2]=ce;
        if(fabs(dot(si,k)-ct)>1e-3) si[1]=-siy;
        return si;
    }

    void dt(double xx,const std::vector<double> &yy,std::vector<double> &dyy) const {

        static unsigned i;
        static double G1m,G2m,Cw0,Cw1,cosI,cost0,cost1,cose0,cose1,de1k2,de1s0,de1s1,
            e12,e14,e16,one21,rone21,a1,n1,Ra0,Ra1,Ra05,Ra15,K0,K1,
            ae1,ae13,a2,ae2,ae23,gamma,one21_5,rone21_3,rone21_9,rone21_13,f1,f2,f4,f5,
            w0,w1,J20,J21,alpha10,alpha11,alpha20,alpha21,w0n1,w1n1,
            dGa,dea,dL0a,dL1a,dL0b,dL1b,det0a,det1a,det0b,det1b,degr,dL0t,dL1t;
        static std::vector<double> G1(3),G2(3),e1(3),L0(3),L1(3), k1(3),k2(3),s0(3),s1(3),
            ck1k2(3),ce1k1(3),ce1k2(3),ck1s0(3),ck1s1(3),ck2s0(3),ck2s1(3);

        for(i=0;i<3;i++) {
            G1[i] = yy[  i]; G2[i] = yy[ 3+i]; e1[i] = yy[6+i];
            L0[i] = yy[9+i]; L1[i] = yy[12+i];
        }

        G1m = norm(G1); G2m = norm(G2);
        Cw0 = norm(L0); Cw1 = norm(L1);
        k1 = div(G1,G1m); k2 = div(G2,G2m);
        s0 = div(L0,Cw0); s1 = div(L1,Cw1);
        cosI  = dot(k1,k2); de1s1 = dot(e1,s1);
        cost0 = dot(k1,s0); cost1 = dot(k1,s1);
        cose0 = dot(k2,s0); cose1 = dot(k2,s1);
        de1k2 = dot(e1,k2); de1s0 = dot(e1,s0);
        ce1k1 = cross(e1,k1); ce1k2 = cross(e1,k2);
        ck1s0 = cross(k1,s0); ck1s1 = cross(k1,s1);
        ck2s0 = cross(k2,s0); ck2s1 = cross(k2,s1);
        ck1k2 = cross(k1,k2);

        e12 = dot(e1,e1); e14 = e12*e12; e16 = e12*e14;
        one21 = 1.-e12; rone21 = sqrt(one21);
        a1 = (G1m*G1m)/(beta1*beta1*mu1*one21);
        n1 = sqrt(mu1/(a1*a1*a1));
        Ra0 = R0/a1; Ra05 = Ra0*Ra0*Ra0*Ra0*Ra0;
        Ra1 = R1/a1; Ra15 = Ra1*Ra1*Ra1*Ra1*Ra1;
        K0 = K0c*Ra05/a1; K1 = K1c*Ra15/a1;
   
        ae1 = a1*rone21; ae13 = ae1*ae1*ae1;
        a2 = (G2m*G2m)/(beta2*beta2*mu2*(1.-dot(e2,e2)));
        ae2 = a2*sqrt(1.-dot(e2,e2)); ae23 = ae2*ae2*ae2;
        gamma = (gammac*a1*a1)/ae23;
        one21_5 = one21*one21*one21* one21*one21;
        rone21_3 = rone21*rone21*rone21;
        rone21_9 = rone21_3*rone21_3*rone21_3;
        rone21_13 = rone21_9 * one21*one21;
        f1 = (1. + 3.*e12 + 0.375*e14)/rone21_9;
        f2 = (1. + 7.5*e12 + 5.625*e14 + 0.3125*e16)/(one21_5*one21);
        f4 = (1. + 1.5*e12 + 0.125*e14)/one21_5;
        f5 = (1. + 3.75*e12 + 1.875*e14 + 0.078125*e16)/rone21_13;

        w0 = Cw0/C0; w1 = Cw1/C1;
        if(J20_const>0) J20 = J20_const;
        else J20 = J20c*w0*w0; 
        if(J21_const>0) J21 = J21_const;
        else J21 = J21c*w1*w1;
        alpha10 = (alpha10c*J20)/ae13;  alpha11 = (alpha11c*J21)/ae13;
        alpha20 = (alpha20c*J20)/ae23;  alpha21 = (alpha21c*J21)/ae23;
        w0n1=w0/n1; w1n1=w1/n1;

        for(i=0;i<3;i++) {
            dGa =  (gamma*one21*cosI*ck1k2[i]) - (5.*gamma*de1k2*ce1k2[i]);
            dea = (gamma*one21/G1m)*( (cosI*ce1k2[i]) - (2.*ce1k1[i]) - (5.*de1k2*ck1k2[i]) );
            dL0a = -alpha10*cost0*ck1s0[i];
            dL1a = -alpha11*cost1*ck1s1[i];
            dL0b = -alpha20*cose0*ck2s0[i];
            dL1b = -alpha21*cose1*ck2s1[i];
            det0a = (-7.5*k20*n1)*(m1/m0)*Ra05*f4*ce1k1[i];
            det1a = (-7.5*k21*n1)*(m0/m1)*Ra15*f4*ce1k1[i];
            det0b = (-K0/(beta1*a1*a1))*( (f4*w0n1*0.5*de1s0*k1[i]) - ((5.5*f4*cost0*w0n1 - 9.*f5)*e1[i]) );
            det1b = (-K1/(beta1*a1*a1))*( (f4*w1n1*0.5*de1s1*k1[i]) - ((5.5*f4*cost1*w1n1 - 9.*f5)*e1[i]) );
            dL0t = K0*n1*(((f4*rone21*0.5*w0n1*(s0[i]-(cost0*k1[i])))-(f1*w0n1*s0[i]))+((f2*k1[i])+((de1s0*(6.+e12)*w0n1)/(4.*rone21_9)*e1[i])));
            dL1t = K1*n1*(((f4*rone21*0.5*w1n1*(s1[i]-(cost1*k1[i])))-(f1*w1n1*s1[i]))+((f2*k1[i])+((de1s1*(6.+e12)*w1n1)/(4.*rone21_9)*e1[i])));
            degr = ((-degrc*n1)/(a1*one21))*ce1k1[i];
            dyy[   i] =  dGa - dL0a-dL1a - dL0t-dL1t;
            dyy[ 3+i] = -dGa - dL0b-dL0b;
            dyy[ 6+i] =  dea + det0a+det1a + det0b+det1b + degr;
            dyy[ 9+i] = dL0a+dL0b + dL0t;
            dyy[12+i] = dL1a+dL1b + dL1t;
        }
    }

    void print_line() const {
        using namespace kctf_const;
        static std::vector<double> G1(3),G2(3),e1(3),L0(3),L1(3), k1(3),k2(3),s0(3),s1(3);
        for(int i=0;i<3;i++) {
            G1[i] = y[  i]; G2[i] = y[ 3+i]; e1[i] = y[6+i];
            L0[i] = y[9+i]; L1[i] = y[12+i];
        }
        double G1m = norm(G1), G2m = norm(G2);
        double a1 = dot(G1,G1)/(beta1*beta1*mu1*(1.-dot(e1,e1)));
        k1 = div(G1,G1m); k2 = div(G2,G2m);
        double I = acos(dot(k1,k2));
        double Cw0 =norm(L0);
        double P0 = (2.*M_PI*C0)/(Cw0*day);
        s0 = div(L0,Cw0);
        double theta0 = acos(dot(k1,s0));
        double phi0 = acos((dot(k2,s0)-cos(I)*cos(theta0))/(sin(I)*sin(theta0)));
        double Cw1 = norm(L1);
        double P1 = (2.*M_PI*C1)/(Cw1*day);
        s1 = div(L1,Cw1);
        double theta1 = acos(dot(k1,s1));
        double phi1 = acos((dot(k2,s1)-cos(I)*cos(theta1))/(sin(I)*sin(theta1)));
        double leak = norm(sub(G0,add(add(G1,G2),add(L0,L1))))/norm(G0);
        printf("%.5e %.5e %8.5f %8.6f %8.5f %7.4f %7.4f %8.4f %8.4f %8.4f %8.4f %.5e\n",
                x/year,leak,
                a1/AU,norm(e1),I/d2r,
                P0,P1,theta0/d2r,theta1/d2r,phi0/d2r,phi1/d2r
                ,J21c*Cw1*Cw1/(C1*C1)
              );
    }

    bool is_sane(bool verbose=false) {
        double e1m = sqrt(y[6]*y[6] + y[7]*y[7] + y[8]*y[8]);
        if(mine>0 and mine>e1m) {
            if(verbose) printf("# Inner orbit circularized!\n");
            return false;
        }
        double G12 = (y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
        double a1 = G12/(beta1*beta1*mu1*(1.-(e1m*e1m)));
        double rp = a1*(1.-e1m);
        if(rp<mindist) {
            if(verbose) printf("# Inner orbit hit roche limit!\n");
            return false;
        }
        return true;
    }

    void evolve(double tgoal) {
        dp.integrate(*this,tgoal);
    }

};

#endif
