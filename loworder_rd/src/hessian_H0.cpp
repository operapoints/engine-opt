#include "LORd.hpp"

void LORd_model::make_H0(const LORd_model_params_t params){
    auto& H0 = this->H0;
    // Local copies allow register hits interleaved with L1 writes
    // Restrict condition is trivially true

    constexpr double pi = 3.14159265358979323846;

    // ---------------- Fields copied directly from params ----------------
    const double EC = params.EC;       // Compressor shaft Young's modulus
    // params.ES is intentionally unused (center shaft Young's modulus)
    const double ET = params.ET;       // Turbine shaft Young's modulus
    const double LC = params.LC;       // Compressor overhang length
    const double LS = params.LS;       // Center shaft length
    const double LT = params.LT;       // Turbine overhang length
    const double mC = params.mC;       // Compressor mass
    const double mT = params.mT;       // Turbine mass
    const double JdC = params.JdC;     // Compressor diametral moment
    const double JdT = params.JdT;     // Turbine diametral moment
    const double JpC = params.JpC;     // Compressor polar moment
    const double JpT = params.JpT;     // Turbine polar moment
    const double KbyC = params.KbyC;   // Compressor end bearing y axis stiffness
    const double KbyT = params.KbyT;   // Turbine end bearing y axis stiffness
    const double KbzC = params.KbzC;   // Compressor end bearing z axis stiffness
    const double KbzT = params.KbzT;   // Turbine end bearing z axis stiffness
    const double Omega = params.Omega; // Shaft speed, positive means CW looking at the compressor

    // ---------------- Local geometry copies ----------------
    const double rC = params.rC;       // Compressor shaft radius
    const double rS = params.rS;       // Center shaft radius
    const double rT = params.rT;       // Turbine shaft radius
    const double rhoS = params.rhoS;   // Center shaft density

    // ---------------- Derived quantities (solid circular cylinder shafts) ----------------
    // Second moment of area (bending MoI) of a solid circular section:
    //   I = pi * r^4 / 4
    const double IC = 0.25 * pi * rC*rC*rC*rC;   // Compressor shaft MoI
    const double IT = 0.25 * pi * rT*rT*rT*rT;   // Turbine shaft MoI

    // Center shaft mass: density * cross-section area * length
    //   mS = rhoS * pi * rS^2 * LS
    const double mS = rhoS * pi * rS*rS * LS;    // Center shaft mass

    // Center shaft polar moment of inertia (mass MoI about spin axis), solid cylinder:
    //   JpS = (1/2) * mS * rS^2
    const double JpS = 0.5 * mS * rS*rS;         // Center shaft polar moment

    // Center shaft diametral moment of inertia (about a transverse axis through its
    // centre of mass), solid cylinder of length LS:
    //   JdS = mS * (3*rS^2 + LS^2) / 12
    const double JdS = mS * (3.0*rS*rS + LS*LS) / 12.0; // Center shaft diametral moment


    const double uc = LC/LS;
    const double ut = LT/LS;
    const double rc = 1. + uc;
    const double rc3 = 3. + uc;
    const double rt = 1. + ut;
    const double rt3 = 3. + ut;
    const double kc = (3.*EC*IC)/(LC*LC*LC);
    const double kt = (3.*ET*IT)/(LT*LT*LT);
    const double gc = (0.375*JpC)/(LC*LC);
    const double gt = (0.375*JpT)/(LT*LT);
    const double dc = (0.75*JdC)/(LC*LC);
    const double dt = (0.75*JdT)/(LT*LT);
    const double gs = (0.5*JpS)/(LS*LS);
    const double ds = JdS/(LS*LS);
    const double ms4 = 0.25*mS;
    const double kYY = kc*(rc*rc) + kt*(ut*ut);
    const double kYZ = kc*rc*uc + kt*rt*ut;
    const double kZZ = kt*(rt*rt) + kc*(uc*uc);
    const double gpYY = 0.33333333333333333*(gc*(rc3*rc3) + gt*(ut*ut));
    const double gpYZ = 0.33333333333333333*(gc*rc3*uc + gt*rt3*ut);
    const double gpZZ = 0.33333333333333333*(gt*(rt3*rt3) + gc*(uc*uc));
    const double dmYY = 0.33333333333333333*(dc*(rc3*rc3) + dt*(ut*ut));
    const double dmYZ = 0.33333333333333333*(dc*rc3*uc + dt*rt3*ut);
    const double dmZZ = 0.33333333333333333*(dt*(rt3*rt3) + dc*(uc*uc));

    const double c1 = 0.33333333333333333*Omega*(3.*gs + gc*(rc3*rc3) + gt*(ut*ut));
    const double c2 = -0.33333333333333333*Omega*(3.*gs + gc*rc3*uc + gt*rt3*ut);
    const double c3 = 0.33333333333333333*Omega*(3.*gs + gt*(rt3*rt3) + gc*(uc*uc));
    const double c4 = -0.33333333333333333*Omega*(3.*gs + gc*(rc3*rc3) + gt*(ut*ut));
    const double c5 = 0.33333333333333333*Omega*(3.*gs + gc*rc3*uc + gt*rt3*ut);
    const double c6 = -0.33333333333333333*Omega*(3.*gs + gt*(rt3*rt3) + gc*(uc*uc));
    const double c7 = ds + 0.25*mS + 0.33333333333333333*dc*(rc3*rc3) + 0.33333333333333333*dt*(ut*ut);
    const double c8 = -1.*ds + 0.083333333333333333*(3.*mS - 4.*(dc*rc3*uc + dt*rt3*ut));
    const double c9 = ds + 0.083333333333333333*(3.*mS + 4.*dt*(rt3*rt3) + 4.*dc*(uc*uc));

    H0(0, 0) = -1.*kc;
    H0(0, 1) = kc*rc;
    H0(0, 2) = -1.*kc*uc;
    H0(0, 12) = 3.*gc*Omega;
    H0(0, 13) = -1.*gc*Omega*rc3;
    H0(0, 14) = gc*Omega*uc;
    H0(1, 0) = kc*rc;
    H0(1, 1) = -1.*KbyC - 1.*kc*(rc*rc) - 1.*kt*(ut*ut);
    H0(1, 2) = kYZ;
    H0(1, 3) = -1.*kt*ut;
    H0(1, 12) = -1.*gc*Omega*rc3;
    H0(1, 13) = c1;
    H0(1, 14) = c2;
    H0(1, 15) = gt*Omega*ut;
    H0(2, 0) = -1.*kc*uc;
    H0(2, 1) = kYZ;
    H0(2, 2) = -1.*KbyT - 1.*kt*(rt*rt) - 1.*kc*(uc*uc);
    H0(2, 3) = kt*rt;
    H0(2, 12) = gc*Omega*uc;
    H0(2, 13) = c2;
    H0(2, 14) = c3;
    H0(2, 15) = -1.*gt*Omega*rt3;
    H0(3, 1) = -1.*kt*ut;
    H0(3, 2) = kt*rt;
    H0(3, 3) = -1.*kt;
    H0(3, 13) = gt*Omega*ut;
    H0(3, 14) = -1.*gt*Omega*rt3;
    H0(3, 15) = 3.*gt*Omega;
    H0(4, 4) = -1.*kc;
    H0(4, 5) = kc*rc;
    H0(4, 6) = -1.*kc*uc;
    H0(4, 8) = -3.*gc*Omega;
    H0(4, 9) = gc*Omega*rc3;
    H0(4, 10) = -1.*gc*Omega*uc;
    H0(5, 4) = kc*rc;
    H0(5, 5) = -1.*KbzC - 1.*kc*(rc*rc) - 1.*kt*(ut*ut);
    H0(5, 6) = kYZ;
    H0(5, 7) = -1.*kt*ut;
    H0(5, 8) = gc*Omega*rc3;
    H0(5, 9) = c4;
    H0(5, 10) = c5;
    H0(5, 11) = -1.*gt*Omega*ut;
    H0(6, 4) = -1.*kc*uc;
    H0(6, 5) = kYZ;
    H0(6, 6) = -1.*KbzT - 1.*kt*(rt*rt) - 1.*kc*(uc*uc);
    H0(6, 7) = kt*rt;
    H0(6, 8) = -1.*gc*Omega*uc;
    H0(6, 9) = c5;
    H0(6, 10) = c6;
    H0(6, 11) = gt*Omega*rt3;
    H0(7, 5) = -1.*kt*ut;
    H0(7, 6) = kt*rt;
    H0(7, 7) = -1.*kt;
    H0(7, 9) = -1.*gt*Omega*ut;
    H0(7, 10) = gt*Omega*rt3;
    H0(7, 11) = -3.*gt*Omega;
    H0(8, 4) = -3.*gc*Omega;
    H0(8, 5) = gc*Omega*rc3;
    H0(8, 6) = -1.*gc*Omega*uc;
    H0(8, 8) = 3.*dc + mC;
    H0(8, 9) = -1.*dc*rc3;
    H0(8, 10) = dc*uc;
    H0(9, 4) = gc*Omega*rc3;
    H0(9, 5) = c4;
    H0(9, 6) = c5;
    H0(9, 7) = -1.*gt*Omega*ut;
    H0(9, 8) = -1.*dc*rc3;
    H0(9, 9) = c7;
    H0(9, 10) = c8;
    H0(9, 11) = dt*ut;
    H0(10, 4) = -1.*gc*Omega*uc;
    H0(10, 5) = c5;
    H0(10, 6) = c6;
    H0(10, 7) = gt*Omega*rt3;
    H0(10, 8) = dc*uc;
    H0(10, 9) = c8;
    H0(10, 10) = c9;
    H0(10, 11) = -1.*dt*rt3;
    H0(11, 5) = -1.*gt*Omega*ut;
    H0(11, 6) = gt*Omega*rt3;
    H0(11, 7) = -3.*gt*Omega;
    H0(11, 9) = dt*ut;
    H0(11, 10) = -1.*dt*rt3;
    H0(11, 11) = 3.*dt + mT;
    H0(12, 0) = 3.*gc*Omega;
    H0(12, 1) = -1.*gc*Omega*rc3;
    H0(12, 2) = gc*Omega*uc;
    H0(12, 12) = 3.*dc + mC;
    H0(12, 13) = -1.*dc*rc3;
    H0(12, 14) = dc*uc;
    H0(13, 0) = -1.*gc*Omega*rc3;
    H0(13, 1) = c1;
    H0(13, 2) = c2;
    H0(13, 3) = gt*Omega*ut;
    H0(13, 12) = -1.*dc*rc3;
    H0(13, 13) = c7;
    H0(13, 14) = c8;
    H0(13, 15) = dt*ut;
    H0(14, 0) = gc*Omega*uc;
    H0(14, 1) = c2;
    H0(14, 2) = c3;
    H0(14, 3) = -1.*gt*Omega*rt3;
    H0(14, 12) = dc*uc;
    H0(14, 13) = c8;
    H0(14, 14) = c9;
    H0(14, 15) = -1.*dt*rt3;
    H0(15, 1) = gt*Omega*ut;
    H0(15, 2) = -1.*gt*Omega*rt3;
    H0(15, 3) = 3.*gt*Omega;
    H0(15, 13) = dt*ut;
    H0(15, 14) = -1.*dt*rt3;
    H0(15, 15) = 3.*dt + mT;
}