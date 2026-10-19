// Standalone reference test: no project header required.
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <cassert>
#include <cmath>
#include <iostream>

struct LORd_model_params_t {
    double EC;
    double ES;
    double ET;
    double rC;
    double rS;
    double rT;
    double LC;
    double LS;
    double LT;
    double rhoS;
    double mC;
    double mT;
    double JdC;
    double JdT;
    double JpC;
    double JpT;
    double KbyC;
    double KbyT;
    double KbzC;
    double KbzT;
    double Omega;
};

struct LORd_model {
    Eigen::Matrix<double,16,16> H0;
    void make_H0(const LORd_model_params_t params);
};

void LORd_model::make_H0(const LORd_model_params_t params) {
    auto& H0 = this->H0;
    const double EC = params.EC;
    const double ES = params.ES;
    const double ET = params.ET;
    const double rC = params.rC;
    const double rS = params.rS;
    const double rT = params.rT;
    const double LC = params.LC;
    const double LS = params.LS;
    const double LT = params.LT;
    const double rhoS = params.rhoS;
    const double mC = params.mC;
    const double mT = params.mT;
    const double JdC = params.JdC;
    const double JdT = params.JdT;
    const double JpC = params.JpC;
    const double JpT = params.JpT;
    const double KbyC = params.KbyC;
    const double KbyT = params.KbyT;
    const double KbzC = params.KbzC;
    const double KbzT = params.KbzT;
    const double Omega = params.Omega;

    // Geometry, in dependency order.
    const double AS = (3.141592653589793238462643383279502884 * ((rS) * (rS)));
    const double IC = ((1.0 / 4.0) * 3.141592653589793238462643383279502884 * ((rC) * (rC) * (rC) * (rC)));
    const double IS = ((1.0 / 4.0) * 3.141592653589793238462643383279502884 * ((rS) * (rS) * (rS) * (rS)));
    const double IT = ((1.0 / 4.0) * 3.141592653589793238462643383279502884 * ((rT) * (rT) * (rT) * (rT)));
    const double muS = (AS * rhoS);
    const double jdS = (IS * rhoS);
    const double jpS = (2.0 * jdS);

    // Slope stiffness blocks; no symbolic inversion.
    Eigen::Matrix4d Ktt = Eigen::Matrix4d::Zero();
    Eigen::Matrix4d Ktu = Eigen::Matrix4d::Zero();
    const double lorPre0 = (1.0 / LC);
    const double lorPre1 = (1.0 / ((LC) * (LC)));
    const double lorPre2 = (1.0 / LS);
    const double lorPre3 = (1.0 / ((LS) * (LS)));
    const double lorPre4 = (1.0 / LT);
    const double lorPre5 = (1.0 / ((LT) * (LT)));
    const double lorPre6 = (ES * IS * lorPre2);
    const double lorPre7 = (2.0 * lorPre6);
    const double lorPre8 = (2.0 * EC * IC * lorPre0);
    const double lorPre9 = (2.0 * ET * IT * lorPre4);
    const double lorPre10 = (6.0 * EC * IC * lorPre1);
    const double lorPre11 = (-6.0 * ET * IT * lorPre5);
    Ktt(0, 0) = (4.0 * EC * IC * lorPre0);
    Ktt(0, 1) = lorPre8;
    Ktt(1, 0) = lorPre8;
    Ktt(1, 1) = (4.0 * ((EC * IC * lorPre0) + lorPre6));
    Ktt(1, 2) = lorPre7;
    Ktt(2, 1) = lorPre7;
    Ktt(2, 2) = (4.0 * ((ET * IT * lorPre4) + lorPre6));
    Ktt(2, 3) = lorPre9;
    Ktt(3, 2) = lorPre9;
    Ktt(3, 3) = (4.0 * ET * IT * lorPre4);
    Ktu(0, 0) = lorPre10;
    Ktu(0, 1) = (-6.0 * EC * IC * lorPre1);
    Ktu(1, 0) = lorPre10;
    Ktu(1, 1) = (-6.0 * ((EC * IC * lorPre1) + (-1.0 * ES * IS * lorPre3)));
    Ktu(1, 2) = (-6.0 * ES * IS * lorPre3);
    Ktu(2, 1) = (6.0 * ES * IS * lorPre3);
    Ktu(2, 2) = (-6.0 * ((ES * IS * lorPre3) + (-1.0 * ET * IT * lorPre5)));
    Ktu(2, 3) = lorPre11;
    Ktu(3, 2) = (6.0 * ET * IT * lorPre5);
    Ktu(3, 3) = lorPre11;
    Eigen::LLT<Eigen::Matrix4d> slopeFactor(Ktt);
    assert(slopeFactor.info() == Eigen::Success);
    const Eigen::Matrix4d S = -slopeFactor.solve(Ktu);

    const double s00 = S(0, 0);
    const double s01 = S(0, 1);
    const double s02 = S(0, 2);
    const double s03 = S(0, 3);
    const double s10 = S(1, 0);
    const double s11 = S(1, 1);
    const double s12 = S(1, 2);
    const double s13 = S(1, 3);
    const double s20 = S(2, 0);
    const double s21 = S(2, 1);
    const double s22 = S(2, 2);
    const double s23 = S(2, 3);
    const double s30 = S(3, 0);
    const double s31 = S(3, 1);
    const double s32 = S(3, 2);
    const double s33 = S(3, 3);

    // CSE across all 30 independent matrix coefficients.
    const double lorCoeff0 = (1.0 / LC);
    const double lorCoeff1 = (1.0 / ((LC) * (LC)));
    const double lorCoeff2 = (1.0 / ((LC) * (LC) * (LC)));
    const double lorCoeff3 = (1.0 / LS);
    const double lorCoeff4 = (1.0 / ((LS) * (LS)));
    const double lorCoeff5 = ((LS) * (LS));
    const double lorCoeff6 = (1.0 / ((LS) * (LS) * (LS)));
    const double lorCoeff7 = ((LS) * (LS) * (LS));
    const double lorCoeff8 = (1.0 / LT);
    const double lorCoeff9 = (1.0 / ((LT) * (LT)));
    const double lorCoeff10 = (1.0 / ((LT) * (LT) * (LT)));
    const double lorCoeff11 = ((s00) * (s00));
    const double lorCoeff12 = ((s01) * (s01));
    const double lorCoeff13 = ((s02) * (s02));
    const double lorCoeff14 = ((s03) * (s03));
    const double lorCoeff15 = ((s10) * (s10));
    const double lorCoeff16 = ((s11) * (s11));
    const double lorCoeff17 = ((s12) * (s12));
    const double lorCoeff18 = ((s13) * (s13));
    const double lorCoeff19 = ((s20) * (s20));
    const double lorCoeff20 = ((s21) * (s21));
    const double lorCoeff21 = ((s22) * (s22));
    const double lorCoeff22 = ((s23) * (s23));
    const double lorCoeff23 = ((s30) * (s30));
    const double lorCoeff24 = ((s31) * (s31));
    const double lorCoeff25 = ((s32) * (s32));
    const double lorCoeff26 = ((s33) * (s33));
    const double lorCoeff27 = (18.0 * jpS * lorCoeff3);
    const double lorCoeff28 = (252.0 * jdS * lorCoeff3);
    const double lorCoeff29 = (78.0 * LS * muS);
    const double lorCoeff30 = (3.0 * EC * IC * lorCoeff2);
    const double lorCoeff31 = (3.0 * ES * IS * lorCoeff6);
    const double lorCoeff32 = (3.0 * ET * IT * lorCoeff10);
    const double lorCoeff33 = (3.0 * EC * IC * lorCoeff1 * s00);
    const double lorCoeff34 = (-3.0 * EC * IC * lorCoeff1 * s01);
    const double lorCoeff35 = (3.0 * EC * IC * lorCoeff1 * s02);
    const double lorCoeff36 = (3.0 * EC * IC * lorCoeff1 * s03);
    const double lorCoeff37 = (3.0 * EC * IC * lorCoeff1 * s10);
    const double lorCoeff38 = (-3.0 * EC * IC * lorCoeff1 * s11);
    const double lorCoeff39 = (3.0 * EC * IC * lorCoeff1 * s12);
    const double lorCoeff40 = (3.0 * EC * IC * lorCoeff1 * s13);
    const double lorCoeff41 = (-3.0 * ES * IS * lorCoeff4 * s10);
    const double lorCoeff42 = (3.0 * ES * IS * lorCoeff4 * s11);
    const double lorCoeff43 = (-3.0 * ES * IS * lorCoeff4 * s12);
    const double lorCoeff44 = (-3.0 * ES * IS * lorCoeff4 * s13);
    const double lorCoeff45 = (-3.0 * ES * IS * lorCoeff4 * s20);
    const double lorCoeff46 = (3.0 * ES * IS * lorCoeff4 * s21);
    const double lorCoeff47 = (-3.0 * ES * IS * lorCoeff4 * s22);
    const double lorCoeff48 = (-3.0 * ES * IS * lorCoeff4 * s23);
    const double mc00 = ((1.0 / 210.0) * ((210.0 * JdC * lorCoeff11) + (210.0 * JdT * lorCoeff23) + (28.0 * jdS * lorCoeff15 * LS) + (28.0 * jdS * lorCoeff19 * LS) + (210.0 * mC) + (2.0 * lorCoeff15 * lorCoeff7 * muS) + (2.0 * lorCoeff19 * lorCoeff7 * muS) + (-14.0 * jdS * LS * s10 * s20) + (-3.0 * lorCoeff7 * muS * s10 * s20)));
    const double mc01 = ((1.0 / 420.0) * ((420.0 * JdC * s00 * s01) + (42.0 * jdS * s10) + (22.0 * lorCoeff5 * muS * s10) + (56.0 * jdS * LS * s10 * s11) + (4.0 * lorCoeff7 * muS * s10 * s11) + (42.0 * jdS * s20) + (-13.0 * lorCoeff5 * muS * s20) + (-14.0 * jdS * LS * s11 * s20) + (-3.0 * lorCoeff7 * muS * s11 * s20) + (-14.0 * jdS * LS * s10 * s21) + (-3.0 * lorCoeff7 * muS * s10 * s21) + (56.0 * jdS * LS * s20 * s21) + (4.0 * lorCoeff7 * muS * s20 * s21) + (420.0 * JdT * s30 * s31)));
    const double mc02 = ((1.0 / 420.0) * ((420.0 * JdC * s00 * s02) + (-42.0 * jdS * s10) + (13.0 * lorCoeff5 * muS * s10) + (56.0 * jdS * LS * s10 * s12) + (4.0 * lorCoeff7 * muS * s10 * s12) + (-42.0 * jdS * s20) + (-22.0 * lorCoeff5 * muS * s20) + (-14.0 * jdS * LS * s12 * s20) + (-3.0 * lorCoeff7 * muS * s12 * s20) + (-14.0 * jdS * LS * s10 * s22) + (-3.0 * lorCoeff7 * muS * s10 * s22) + (56.0 * jdS * LS * s20 * s22) + (4.0 * lorCoeff7 * muS * s20 * s22) + (420.0 * JdT * s30 * s32)));
    const double mc03 = ((1.0 / 420.0) * ((420.0 * JdC * s00 * s03) + (56.0 * jdS * LS * s10 * s13) + (4.0 * lorCoeff7 * muS * s10 * s13) + (-14.0 * jdS * LS * s13 * s20) + (-3.0 * lorCoeff7 * muS * s13 * s20) + (-14.0 * jdS * LS * s10 * s23) + (-3.0 * lorCoeff7 * muS * s10 * s23) + (56.0 * jdS * LS * s20 * s23) + (4.0 * lorCoeff7 * muS * s20 * s23) + (420.0 * JdT * s30 * s33)));
    const double mc11 = ((1.0 / 210.0) * ((210.0 * JdC * lorCoeff12) + (210.0 * JdT * lorCoeff24) + lorCoeff28 + lorCoeff29 + (28.0 * jdS * lorCoeff16 * LS) + (28.0 * jdS * lorCoeff20 * LS) + (2.0 * lorCoeff16 * lorCoeff7 * muS) + (2.0 * lorCoeff20 * lorCoeff7 * muS) + (42.0 * jdS * s11) + (22.0 * lorCoeff5 * muS * s11) + (42.0 * jdS * s21) + (-13.0 * lorCoeff5 * muS * s21) + (-14.0 * jdS * LS * s11 * s21) + (-3.0 * lorCoeff7 * muS * s11 * s21)));
    const double mc12 = ((1.0 / 420.0) * ((-504.0 * jdS * lorCoeff3) + (54.0 * LS * muS) + (420.0 * JdC * s01 * s02) + (-42.0 * jdS * s11) + (13.0 * lorCoeff5 * muS * s11) + (42.0 * jdS * s12) + (22.0 * lorCoeff5 * muS * s12) + (56.0 * jdS * LS * s11 * s12) + (4.0 * lorCoeff7 * muS * s11 * s12) + (-42.0 * jdS * s21) + (-22.0 * lorCoeff5 * muS * s21) + (-14.0 * jdS * LS * s12 * s21) + (-3.0 * lorCoeff7 * muS * s12 * s21) + (42.0 * jdS * s22) + (-13.0 * lorCoeff5 * muS * s22) + (-14.0 * jdS * LS * s11 * s22) + (-3.0 * lorCoeff7 * muS * s11 * s22) + (56.0 * jdS * LS * s21 * s22) + (4.0 * lorCoeff7 * muS * s21 * s22) + (420.0 * JdT * s31 * s32)));
    const double mc13 = ((1.0 / 420.0) * ((420.0 * JdC * s01 * s03) + (42.0 * jdS * s13) + (22.0 * lorCoeff5 * muS * s13) + (56.0 * jdS * LS * s11 * s13) + (4.0 * lorCoeff7 * muS * s11 * s13) + (-14.0 * jdS * LS * s13 * s21) + (-3.0 * lorCoeff7 * muS * s13 * s21) + (42.0 * jdS * s23) + (-13.0 * lorCoeff5 * muS * s23) + (-14.0 * jdS * LS * s11 * s23) + (-3.0 * lorCoeff7 * muS * s11 * s23) + (56.0 * jdS * LS * s21 * s23) + (4.0 * lorCoeff7 * muS * s21 * s23) + (420.0 * JdT * s31 * s33)));
    const double mc22 = ((1.0 / 210.0) * ((210.0 * JdC * lorCoeff13) + (210.0 * JdT * lorCoeff25) + lorCoeff28 + lorCoeff29 + (28.0 * jdS * lorCoeff17 * LS) + (28.0 * jdS * lorCoeff21 * LS) + (2.0 * lorCoeff17 * lorCoeff7 * muS) + (2.0 * lorCoeff21 * lorCoeff7 * muS) + (-42.0 * jdS * s12) + (13.0 * lorCoeff5 * muS * s12) + (-42.0 * jdS * s22) + (-22.0 * lorCoeff5 * muS * s22) + (-14.0 * jdS * LS * s12 * s22) + (-3.0 * lorCoeff7 * muS * s12 * s22)));
    const double mc23 = ((1.0 / 420.0) * ((420.0 * JdC * s02 * s03) + (-42.0 * jdS * s13) + (13.0 * lorCoeff5 * muS * s13) + (56.0 * jdS * LS * s12 * s13) + (4.0 * lorCoeff7 * muS * s12 * s13) + (-14.0 * jdS * LS * s13 * s22) + (-3.0 * lorCoeff7 * muS * s13 * s22) + (-42.0 * jdS * s23) + (-22.0 * lorCoeff5 * muS * s23) + (-14.0 * jdS * LS * s12 * s23) + (-3.0 * lorCoeff7 * muS * s12 * s23) + (56.0 * jdS * LS * s22 * s23) + (4.0 * lorCoeff7 * muS * s22 * s23) + (420.0 * JdT * s32 * s33)));
    const double mc33 = ((1.0 / 210.0) * ((210.0 * JdC * lorCoeff14) + (210.0 * JdT * lorCoeff26) + (28.0 * jdS * lorCoeff18 * LS) + (28.0 * jdS * lorCoeff22 * LS) + (210.0 * mT) + (2.0 * lorCoeff18 * lorCoeff7 * muS) + (2.0 * lorCoeff22 * lorCoeff7 * muS) + (-14.0 * jdS * LS * s13 * s23) + (-3.0 * lorCoeff7 * muS * s13 * s23)));
    const double kc00 = (4.0 * ((EC * IC * lorCoeff0 * lorCoeff11) + (EC * IC * lorCoeff0 * lorCoeff15) + (ES * IS * lorCoeff15 * lorCoeff3) + (ES * IS * lorCoeff19 * lorCoeff3) + lorCoeff30 + lorCoeff33 + lorCoeff37 + (ET * IT * lorCoeff19 * lorCoeff8) + (ET * IT * lorCoeff23 * lorCoeff8) + (EC * IC * lorCoeff0 * s00 * s10) + (ES * IS * lorCoeff3 * s10 * s20) + (ET * IT * lorCoeff8 * s20 * s30)));
    const double kc01 = (-2.0 * ((6.0 * EC * IC * lorCoeff2) + lorCoeff33 + lorCoeff34 + lorCoeff37 + lorCoeff38 + lorCoeff41 + lorCoeff45 + (-2.0 * EC * IC * lorCoeff0 * s00 * s01) + (-1.0 * EC * IC * lorCoeff0 * s01 * s10) + (-1.0 * EC * IC * lorCoeff0 * s00 * s11) + (-2.0 * EC * IC * lorCoeff0 * s10 * s11) + (-2.0 * ES * IS * lorCoeff3 * s10 * s11) + (-1.0 * ES * IS * lorCoeff3 * s11 * s20) + (-1.0 * ES * IS * lorCoeff3 * s10 * s21) + (-2.0 * ES * IS * lorCoeff3 * s20 * s21) + (-2.0 * ET * IT * lorCoeff8 * s20 * s21) + (-1.0 * ET * IT * lorCoeff8 * s21 * s30) + (-1.0 * ET * IT * lorCoeff8 * s20 * s31) + (-2.0 * ET * IT * lorCoeff8 * s30 * s31)));
    const double kc02 = (2.0 * (lorCoeff35 + lorCoeff39 + lorCoeff41 + lorCoeff45 + (2.0 * EC * IC * lorCoeff0 * s00 * s02) + (EC * IC * lorCoeff0 * s02 * s10) + (EC * IC * lorCoeff0 * s00 * s12) + (2.0 * EC * IC * lorCoeff0 * s10 * s12) + (2.0 * ES * IS * lorCoeff3 * s10 * s12) + (3.0 * ET * IT * lorCoeff9 * s20) + (ES * IS * lorCoeff3 * s12 * s20) + (ES * IS * lorCoeff3 * s10 * s22) + (2.0 * ES * IS * lorCoeff3 * s20 * s22) + (2.0 * ET * IT * lorCoeff8 * s20 * s22) + (3.0 * ET * IT * lorCoeff9 * s30) + (ET * IT * lorCoeff8 * s22 * s30) + (ET * IT * lorCoeff8 * s20 * s32) + (2.0 * ET * IT * lorCoeff8 * s30 * s32)));
    const double kc03 = (2.0 * (lorCoeff36 + lorCoeff40 + (2.0 * EC * IC * lorCoeff0 * s00 * s03) + (EC * IC * lorCoeff0 * s03 * s10) + (EC * IC * lorCoeff0 * s00 * s13) + (2.0 * EC * IC * lorCoeff0 * s10 * s13) + (2.0 * ES * IS * lorCoeff3 * s10 * s13) + (-3.0 * ET * IT * lorCoeff9 * s20) + (ES * IS * lorCoeff3 * s13 * s20) + (ES * IS * lorCoeff3 * s10 * s23) + (2.0 * ES * IS * lorCoeff3 * s20 * s23) + (2.0 * ET * IT * lorCoeff8 * s20 * s23) + (-3.0 * ET * IT * lorCoeff9 * s30) + (ET * IT * lorCoeff8 * s23 * s30) + (ET * IT * lorCoeff8 * s20 * s33) + (2.0 * ET * IT * lorCoeff8 * s30 * s33)));
    const double kc11 = (4.0 * ((EC * IC * lorCoeff0 * lorCoeff12) + (EC * IC * lorCoeff0 * lorCoeff16) + (ES * IS * lorCoeff16 * lorCoeff3) + (ES * IS * lorCoeff20 * lorCoeff3) + lorCoeff30 + lorCoeff31 + lorCoeff34 + lorCoeff38 + lorCoeff42 + lorCoeff46 + (ET * IT * lorCoeff20 * lorCoeff8) + (ET * IT * lorCoeff24 * lorCoeff8) + (EC * IC * lorCoeff0 * s01 * s11) + (ES * IS * lorCoeff3 * s11 * s21) + (ET * IT * lorCoeff8 * s21 * s31)));
    const double kc12 = (-2.0 * (lorCoeff35 + lorCoeff39 + lorCoeff42 + lorCoeff43 + lorCoeff46 + lorCoeff47 + (6.0 * ES * IS * lorCoeff6) + (-2.0 * EC * IC * lorCoeff0 * s01 * s02) + (-1.0 * EC * IC * lorCoeff0 * s02 * s11) + (-1.0 * EC * IC * lorCoeff0 * s01 * s12) + (-2.0 * EC * IC * lorCoeff0 * s11 * s12) + (-2.0 * ES * IS * lorCoeff3 * s11 * s12) + (-3.0 * ET * IT * lorCoeff9 * s21) + (-1.0 * ES * IS * lorCoeff3 * s12 * s21) + (-1.0 * ES * IS * lorCoeff3 * s11 * s22) + (-2.0 * ES * IS * lorCoeff3 * s21 * s22) + (-2.0 * ET * IT * lorCoeff8 * s21 * s22) + (-3.0 * ET * IT * lorCoeff9 * s31) + (-1.0 * ET * IT * lorCoeff8 * s22 * s31) + (-1.0 * ET * IT * lorCoeff8 * s21 * s32) + (-2.0 * ET * IT * lorCoeff8 * s31 * s32)));
    const double kc13 = (-2.0 * (lorCoeff36 + lorCoeff40 + lorCoeff44 + lorCoeff48 + (-2.0 * EC * IC * lorCoeff0 * s01 * s03) + (-1.0 * EC * IC * lorCoeff0 * s03 * s11) + (-1.0 * EC * IC * lorCoeff0 * s01 * s13) + (-2.0 * EC * IC * lorCoeff0 * s11 * s13) + (-2.0 * ES * IS * lorCoeff3 * s11 * s13) + (3.0 * ET * IT * lorCoeff9 * s21) + (-1.0 * ES * IS * lorCoeff3 * s13 * s21) + (-1.0 * ES * IS * lorCoeff3 * s11 * s23) + (-2.0 * ES * IS * lorCoeff3 * s21 * s23) + (-2.0 * ET * IT * lorCoeff8 * s21 * s23) + (3.0 * ET * IT * lorCoeff9 * s31) + (-1.0 * ET * IT * lorCoeff8 * s23 * s31) + (-1.0 * ET * IT * lorCoeff8 * s21 * s33) + (-2.0 * ET * IT * lorCoeff8 * s31 * s33)));
    const double kc22 = (4.0 * ((EC * IC * lorCoeff0 * lorCoeff13) + (EC * IC * lorCoeff0 * lorCoeff17) + (ES * IS * lorCoeff17 * lorCoeff3) + (ES * IS * lorCoeff21 * lorCoeff3) + lorCoeff31 + lorCoeff32 + lorCoeff43 + lorCoeff47 + (ET * IT * lorCoeff21 * lorCoeff8) + (ET * IT * lorCoeff25 * lorCoeff8) + (EC * IC * lorCoeff0 * s02 * s12) + (3.0 * ET * IT * lorCoeff9 * s22) + (ES * IS * lorCoeff3 * s12 * s22) + (3.0 * ET * IT * lorCoeff9 * s32) + (ET * IT * lorCoeff8 * s22 * s32)));
    const double kc23 = (2.0 * ((-6.0 * ET * IT * lorCoeff10) + lorCoeff44 + lorCoeff48 + (2.0 * EC * IC * lorCoeff0 * s02 * s03) + (EC * IC * lorCoeff0 * s03 * s12) + (EC * IC * lorCoeff0 * s02 * s13) + (2.0 * EC * IC * lorCoeff0 * s12 * s13) + (2.0 * ES * IS * lorCoeff3 * s12 * s13) + (-3.0 * ET * IT * lorCoeff9 * s22) + (ES * IS * lorCoeff3 * s13 * s22) + (3.0 * ET * IT * lorCoeff9 * s23) + (ES * IS * lorCoeff3 * s12 * s23) + (2.0 * ES * IS * lorCoeff3 * s22 * s23) + (2.0 * ET * IT * lorCoeff8 * s22 * s23) + (-3.0 * ET * IT * lorCoeff9 * s32) + (ET * IT * lorCoeff8 * s23 * s32) + (3.0 * ET * IT * lorCoeff9 * s33) + (ET * IT * lorCoeff8 * s22 * s33) + (2.0 * ET * IT * lorCoeff8 * s32 * s33)));
    const double kc33 = (4.0 * ((EC * IC * lorCoeff0 * lorCoeff14) + (EC * IC * lorCoeff0 * lorCoeff18) + (ES * IS * lorCoeff18 * lorCoeff3) + (ES * IS * lorCoeff22 * lorCoeff3) + lorCoeff32 + (ET * IT * lorCoeff22 * lorCoeff8) + (ET * IT * lorCoeff26 * lorCoeff8) + (EC * IC * lorCoeff0 * s03 * s13) + (-3.0 * ET * IT * lorCoeff9 * s23) + (ES * IS * lorCoeff3 * s13 * s23) + (-3.0 * ET * IT * lorCoeff9 * s33) + (ET * IT * lorCoeff8 * s23 * s33)));
    const double wc00 = ((1.0 / 15.0) * ((15.0 * JpC * lorCoeff11) + (15.0 * JpT * lorCoeff23) + (2.0 * jpS * lorCoeff15 * LS) + (2.0 * jpS * lorCoeff19 * LS) + (-1.0 * jpS * LS * s10 * s20)));
    const double wc01 = ((1.0 / 30.0) * ((30.0 * JpC * s00 * s01) + (3.0 * jpS * s10) + (4.0 * jpS * LS * s10 * s11) + (3.0 * jpS * s20) + (-1.0 * jpS * LS * s11 * s20) + (-1.0 * jpS * LS * s10 * s21) + (4.0 * jpS * LS * s20 * s21) + (30.0 * JpT * s30 * s31)));
    const double wc02 = ((1.0 / 30.0) * ((30.0 * JpC * s00 * s02) + (-3.0 * jpS * s10) + (4.0 * jpS * LS * s10 * s12) + (-3.0 * jpS * s20) + (-1.0 * jpS * LS * s12 * s20) + (-1.0 * jpS * LS * s10 * s22) + (4.0 * jpS * LS * s20 * s22) + (30.0 * JpT * s30 * s32)));
    const double wc03 = ((1.0 / 30.0) * ((30.0 * JpC * s00 * s03) + (4.0 * jpS * LS * s10 * s13) + (-1.0 * jpS * LS * s13 * s20) + (-1.0 * jpS * LS * s10 * s23) + (4.0 * jpS * LS * s20 * s23) + (30.0 * JpT * s30 * s33)));
    const double wc11 = ((1.0 / 15.0) * ((15.0 * JpC * lorCoeff12) + (15.0 * JpT * lorCoeff24) + lorCoeff27 + (2.0 * jpS * lorCoeff16 * LS) + (2.0 * jpS * lorCoeff20 * LS) + (3.0 * jpS * s11) + (3.0 * jpS * s21) + (-1.0 * jpS * LS * s11 * s21)));
    const double wc12 = ((1.0 / 30.0) * ((-36.0 * jpS * lorCoeff3) + (30.0 * JpC * s01 * s02) + (-3.0 * jpS * s11) + (3.0 * jpS * s12) + (4.0 * jpS * LS * s11 * s12) + (-3.0 * jpS * s21) + (-1.0 * jpS * LS * s12 * s21) + (3.0 * jpS * s22) + (-1.0 * jpS * LS * s11 * s22) + (4.0 * jpS * LS * s21 * s22) + (30.0 * JpT * s31 * s32)));
    const double wc13 = ((1.0 / 30.0) * ((30.0 * JpC * s01 * s03) + (3.0 * jpS * s13) + (4.0 * jpS * LS * s11 * s13) + (-1.0 * jpS * LS * s13 * s21) + (3.0 * jpS * s23) + (-1.0 * jpS * LS * s11 * s23) + (4.0 * jpS * LS * s21 * s23) + (30.0 * JpT * s31 * s33)));
    const double wc22 = ((1.0 / 15.0) * ((15.0 * JpC * lorCoeff13) + (15.0 * JpT * lorCoeff25) + lorCoeff27 + (2.0 * jpS * lorCoeff17 * LS) + (2.0 * jpS * lorCoeff21 * LS) + (-3.0 * jpS * s12) + (-3.0 * jpS * s22) + (-1.0 * jpS * LS * s12 * s22)));
    const double wc23 = ((1.0 / 30.0) * ((30.0 * JpC * s02 * s03) + (-3.0 * jpS * s13) + (4.0 * jpS * LS * s12 * s13) + (-1.0 * jpS * LS * s13 * s22) + (-3.0 * jpS * s23) + (-1.0 * jpS * LS * s12 * s23) + (4.0 * jpS * LS * s22 * s23) + (30.0 * JpT * s32 * s33)));
    const double wc33 = ((1.0 / 15.0) * ((15.0 * JpC * lorCoeff14) + (15.0 * JpT * lorCoeff26) + (2.0 * jpS * lorCoeff18 * LS) + (2.0 * jpS * lorCoeff22 * LS) + (-1.0 * jpS * LS * s13 * s23)));

    // Zero all entries, then fill the upper triangle and mirror.
    H0.setZero();
    const double lorStore0 = (-1.0 * kc00);
    const double lorStore1 = (-1.0 * kc01);
    const double lorStore2 = (-1.0 * kc02);
    const double lorStore3 = (-1.0 * kc03);
    const double lorStore4 = (-1.0 * kc11);
    const double lorStore5 = (-1.0 * kc12);
    const double lorStore6 = (-1.0 * kc13);
    const double lorStore7 = (-1.0 * kc22);
    const double lorStore8 = (-1.0 * kc23);
    const double lorStore9 = (-1.0 * kc33);
    const double lorStore10 = ((-1.0 / 2.0) * Omega * wc01);
    const double lorStore11 = ((-1.0 / 2.0) * Omega * wc02);
    const double lorStore12 = ((-1.0 / 2.0) * Omega * wc03);
    const double lorStore13 = ((-1.0 / 2.0) * Omega * wc12);
    const double lorStore14 = ((-1.0 / 2.0) * Omega * wc13);
    const double lorStore15 = ((-1.0 / 2.0) * Omega * wc23);
    const double lorStore16 = ((1.0 / 2.0) * Omega * wc01);
    const double lorStore17 = ((1.0 / 2.0) * Omega * wc02);
    const double lorStore18 = ((1.0 / 2.0) * Omega * wc03);
    const double lorStore19 = ((1.0 / 2.0) * Omega * wc12);
    const double lorStore20 = ((1.0 / 2.0) * Omega * wc13);
    const double lorStore21 = ((1.0 / 2.0) * Omega * wc23);
    H0(0, 0) = lorStore0;
    H0(0, 1) = lorStore1;
    H0(1, 0) = H0(0, 1);
    H0(0, 2) = lorStore2;
    H0(2, 0) = H0(0, 2);
    H0(0, 3) = lorStore3;
    H0(3, 0) = H0(0, 3);
    H0(0, 12) = ((1.0 / 2.0) * Omega * wc00);
    H0(12, 0) = H0(0, 12);
    H0(0, 13) = lorStore16;
    H0(13, 0) = H0(0, 13);
    H0(0, 14) = lorStore17;
    H0(14, 0) = H0(0, 14);
    H0(0, 15) = lorStore18;
    H0(15, 0) = H0(0, 15);
    H0(1, 1) = ((-1.0 * KbyC) + lorStore4);
    H0(1, 2) = lorStore5;
    H0(2, 1) = H0(1, 2);
    H0(1, 3) = lorStore6;
    H0(3, 1) = H0(1, 3);
    H0(1, 12) = lorStore16;
    H0(12, 1) = H0(1, 12);
    H0(1, 13) = ((1.0 / 2.0) * Omega * wc11);
    H0(13, 1) = H0(1, 13);
    H0(1, 14) = lorStore19;
    H0(14, 1) = H0(1, 14);
    H0(1, 15) = lorStore20;
    H0(15, 1) = H0(1, 15);
    H0(2, 2) = ((-1.0 * KbyT) + lorStore7);
    H0(2, 3) = lorStore8;
    H0(3, 2) = H0(2, 3);
    H0(2, 12) = lorStore17;
    H0(12, 2) = H0(2, 12);
    H0(2, 13) = lorStore19;
    H0(13, 2) = H0(2, 13);
    H0(2, 14) = ((1.0 / 2.0) * Omega * wc22);
    H0(14, 2) = H0(2, 14);
    H0(2, 15) = lorStore21;
    H0(15, 2) = H0(2, 15);
    H0(3, 3) = lorStore9;
    H0(3, 12) = lorStore18;
    H0(12, 3) = H0(3, 12);
    H0(3, 13) = lorStore20;
    H0(13, 3) = H0(3, 13);
    H0(3, 14) = lorStore21;
    H0(14, 3) = H0(3, 14);
    H0(3, 15) = ((1.0 / 2.0) * Omega * wc33);
    H0(15, 3) = H0(3, 15);
    H0(4, 4) = lorStore0;
    H0(4, 5) = lorStore1;
    H0(5, 4) = H0(4, 5);
    H0(4, 6) = lorStore2;
    H0(6, 4) = H0(4, 6);
    H0(4, 7) = lorStore3;
    H0(7, 4) = H0(4, 7);
    H0(4, 8) = ((-1.0 / 2.0) * Omega * wc00);
    H0(8, 4) = H0(4, 8);
    H0(4, 9) = lorStore10;
    H0(9, 4) = H0(4, 9);
    H0(4, 10) = lorStore11;
    H0(10, 4) = H0(4, 10);
    H0(4, 11) = lorStore12;
    H0(11, 4) = H0(4, 11);
    H0(5, 5) = ((-1.0 * KbzC) + lorStore4);
    H0(5, 6) = lorStore5;
    H0(6, 5) = H0(5, 6);
    H0(5, 7) = lorStore6;
    H0(7, 5) = H0(5, 7);
    H0(5, 8) = lorStore10;
    H0(8, 5) = H0(5, 8);
    H0(5, 9) = ((-1.0 / 2.0) * Omega * wc11);
    H0(9, 5) = H0(5, 9);
    H0(5, 10) = lorStore13;
    H0(10, 5) = H0(5, 10);
    H0(5, 11) = lorStore14;
    H0(11, 5) = H0(5, 11);
    H0(6, 6) = ((-1.0 * KbzT) + lorStore7);
    H0(6, 7) = lorStore8;
    H0(7, 6) = H0(6, 7);
    H0(6, 8) = lorStore11;
    H0(8, 6) = H0(6, 8);
    H0(6, 9) = lorStore13;
    H0(9, 6) = H0(6, 9);
    H0(6, 10) = ((-1.0 / 2.0) * Omega * wc22);
    H0(10, 6) = H0(6, 10);
    H0(6, 11) = lorStore15;
    H0(11, 6) = H0(6, 11);
    H0(7, 7) = lorStore9;
    H0(7, 8) = lorStore12;
    H0(8, 7) = H0(7, 8);
    H0(7, 9) = lorStore14;
    H0(9, 7) = H0(7, 9);
    H0(7, 10) = lorStore15;
    H0(10, 7) = H0(7, 10);
    H0(7, 11) = ((-1.0 / 2.0) * Omega * wc33);
    H0(11, 7) = H0(7, 11);
    H0(8, 8) = mc00;
    H0(8, 9) = mc01;
    H0(9, 8) = H0(8, 9);
    H0(8, 10) = mc02;
    H0(10, 8) = H0(8, 10);
    H0(8, 11) = mc03;
    H0(11, 8) = H0(8, 11);
    H0(9, 9) = mc11;
    H0(9, 10) = mc12;
    H0(10, 9) = H0(9, 10);
    H0(9, 11) = mc13;
    H0(11, 9) = H0(9, 11);
    H0(10, 10) = mc22;
    H0(10, 11) = mc23;
    H0(11, 10) = H0(10, 11);
    H0(11, 11) = mc33;
    H0(12, 12) = mc00;
    H0(12, 13) = mc01;
    H0(13, 12) = H0(12, 13);
    H0(12, 14) = mc02;
    H0(14, 12) = H0(12, 14);
    H0(12, 15) = mc03;
    H0(15, 12) = H0(12, 15);
    H0(13, 13) = mc11;
    H0(13, 14) = mc12;
    H0(14, 13) = H0(13, 14);
    H0(13, 15) = mc13;
    H0(15, 13) = H0(13, 15);
    H0(14, 14) = mc22;
    H0(14, 15) = mc23;
    H0(15, 14) = H0(14, 15);
    H0(15, 15) = mc33;
}

int main() {
    LORd_model model;
    {
        LORd_model_params_t params{};
        params.EC = 2.e11;
        params.ES = 1.9e11;
        params.ET = 1.8e11;
        params.rC = 0.008;
        params.rS = 0.006;
        params.rT = 0.007;
        params.LC = 0.04;
        params.LS = 0.15;
        params.LT = 0.05;
        params.rhoS = 7800.;
        params.mC = 0.2;
        params.mT = 0.15;
        params.JdC = 0.0001;
        params.JdT = 0.00008;
        params.JpC = 0.0002;
        params.JpT = 0.00016;
        params.KbyC = 1.e7;
        params.KbyT = 9.e6;
        params.KbzC = 1.2e7;
        params.KbzT = 1.1e7;
        params.Omega = 1000.;
        model.make_H0(params);
        const double reference[256] = {
            -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, 0., 0., 0., 0., 0., 0., 0., 0., 68.58251166382975, -69.82710869156826, 4.250016776272202, -3.0054197485337064, 3.830768523996932e6, -1.5426023172553241e7, 3.316403775027614e6, -1.7211491264713265e6, 0., 0., 0., 0., 0., 0., 0., 0., -69.82710869156826, 71.22505809244268, -6.127204408993079, 4.729255008118641, -1.9862811966009932e6, 3.316403775027614e6, -1.273146535642568e7, 2.40134277799907e6, 0., 0., 0., 0., 0., 0., 0., 0., 4.250016776272202, -6.127204408993079, 42.31281573948085, -40.43562810675997, 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, 0., 0., 0., 0., 0., 0., 0., 0., -3.0054197485337064, 4.729255008118641, -40.43562810675997, 38.711792847175055, 0., 0., 0., 0., -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, -68.58251166382975, 69.82710869156826, -4.250016776272202, 3.0054197485337064, 0., 0., 0., 0., 0., 0., 0., 0., 3.830768523996932e6, -1.742602317255324e7, 3.316403775027614e6, -1.7211491264713265e6, 69.82710869156826, -71.22505809244268, 6.127204408993079, -4.729255008118641, 0., 0., 0., 0., 0., 0., 0., 0., -1.9862811966009932e6, 3.316403775027614e6, -1.473146535642568e7, 2.40134277799907e6, -4.250016776272202, 6.127204408993079, -42.31281573948085, 40.43562810675997, 0., 0., 0., 0., 0., 0., 0., 0., 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, 3.0054197485337064, -4.729255008118641, 40.43562810675997, -38.711792847175055, 0., 0., 0., 0., 0., 0., 0., 0., -68.58251166382975, 69.82710869156826, -4.250016776272202, 3.0054197485337064, 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, 0., 0., 0., 0., 0., 0., 0., 0., 69.82710869156826, -71.22505809244268, 6.127204408993079, -4.729255008118641, -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, 0., 0., 0., 0., 0., 0., 0., 0., -4.250016776272202, 6.127204408993079, -42.31281573948085, 40.43562810675997, -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, 0., 0., 0., 0., 0., 0., 0., 0., 3.0054197485337064, -4.729255008118641, 40.43562810675997, -38.711792847175055, 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625, 0., 0., 0., 0., 68.58251166382975, -69.82710869156826, 4.250016776272202, -3.0054197485337064, 0., 0., 0., 0., 0., 0., 0., 0., 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, -69.82710869156826, 71.22505809244268, -6.127204408993079, 4.729255008118641, 0., 0., 0., 0., 0., 0., 0., 0., -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, 4.250016776272202, -6.127204408993079, 42.31281573948085, -40.43562810675997, 0., 0., 0., 0., 0., 0., 0., 0., -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, -3.0054197485337064, 4.729255008118641, -40.43562810675997, 38.711792847175055, 0., 0., 0., 0., 0., 0., 0., 0., 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625
        };
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                const double expected = reference[16*r+c];
                const double actual = model.H0(r,c);
                const double error = std::abs(actual - expected);
                const double tolerance = 1e-8 + 1e-10*std::abs(expected);
                if (!std::isfinite(actual) || error > tolerance) {
                    std::cerr << "Hessian mismatch: case 1, entry (" << r << "," << c << ")\n";
                    return 1;
                }
            }
        }
    }
    {
        LORd_model_params_t params{};
        params.EC = 2.e11;
        params.ES = 1.9e11;
        params.ET = 1.8e11;
        params.rC = 0.008;
        params.rS = 0.006;
        params.rT = 0.007;
        params.LC = 0.04;
        params.LS = 0.15;
        params.LT = 0.05;
        params.rhoS = 7800.;
        params.mC = 0.2;
        params.mT = 0.15;
        params.JdC = 0.0001;
        params.JdT = 0.00008;
        params.JpC = 0.0002;
        params.JpT = 0.00016;
        params.KbyC = 1.e7;
        params.KbyT = 9.e6;
        params.KbzC = 1.2e7;
        params.KbzT = 1.1e7;
        params.Omega = 0.;
        model.make_H0(params);
        const double reference[256] = {
            -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3.830768523996932e6, -1.5426023172553241e7, 3.316403775027614e6, -1.7211491264713265e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1.9862811966009932e6, 3.316403775027614e6, -1.273146535642568e7, 2.40134277799907e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3.830768523996932e6, -1.742602317255324e7, 3.316403775027614e6, -1.7211491264713265e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1.9862811966009932e6, 3.316403775027614e6, -1.473146535642568e7, 2.40134277799907e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625
        };
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                const double expected = reference[16*r+c];
                const double actual = model.H0(r,c);
                const double error = std::abs(actual - expected);
                const double tolerance = 1e-8 + 1e-10*std::abs(expected);
                if (!std::isfinite(actual) || error > tolerance) {
                    std::cerr << "Hessian mismatch: case 2, entry (" << r << "," << c << ")\n";
                    return 1;
                }
            }
        }
    }
    {
        LORd_model_params_t params{};
        params.EC = 2.e11;
        params.ES = 1.9e11;
        params.ET = 1.8e11;
        params.rC = 0.008;
        params.rS = 0.006;
        params.rT = 0.007;
        params.LC = 0.04;
        params.LS = 0.15;
        params.LT = 0.05;
        params.rhoS = 7800.;
        params.mC = 0.2;
        params.mT = 0.15;
        params.JdC = 0.0001;
        params.JdT = 0.00008;
        params.JpC = 0.0002;
        params.JpT = 0.00016;
        params.KbyC = 1.e7;
        params.KbyT = 9.e6;
        params.KbzC = 1.2e7;
        params.KbzT = 1.1e7;
        params.Omega = -1700.;
        model.make_H0(params);
        const double reference[256] = {
            -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, 0., 0., 0., 0., 0., 0., 0., 0., -116.59026982851056, 118.70608477566601, -7.225028519662741, 5.109213572507301, 3.830768523996932e6, -1.5426023172553241e7, 3.316403775027614e6, -1.7211491264713265e6, 0., 0., 0., 0., 0., 0., 0., 0., 118.70608477566601, -121.08259875715257, 10.416247495288234, -8.03973351380169, -1.9862811966009932e6, 3.316403775027614e6, -1.273146535642568e7, 2.40134277799907e6, 0., 0., 0., 0., 0., 0., 0., 0., -7.225028519662741, 10.416247495288234, -71.93178675711744, 68.74056778149193, 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, 0., 0., 0., 0., 0., 0., 0., 0., 5.109213572507301, -8.03973351380169, 68.74056778149193, -65.81004784019758, 0., 0., 0., 0., -2.7784985207055937e6, 3.830768523996932e6, -1.9862811966009932e6, 934011.193309633, 116.59026982851056, -118.70608477566601, 7.225028519662741, -5.109213572507301, 0., 0., 0., 0., 0., 0., 0., 0., 3.830768523996932e6, -1.742602317255324e7, 3.316403775027614e6, -1.7211491264713265e6, -118.70608477566601, 121.08259875715257, -10.416247495288234, 8.03973351380169, 0., 0., 0., 0., 0., 0., 0., 0., -1.9862811966009932e6, 3.316403775027614e6, -1.473146535642568e7, 2.40134277799907e6, 7.225028519662741, -10.416247495288234, 71.93178675711744, -68.74056778149193, 0., 0., 0., 0., 0., 0., 0., 0., 934011.193309633, -1.7211491264713265e6, 2.40134277799907e6, -1.6142048448373785e6, -5.109213572507301, 8.03973351380169, -68.74056778149193, 65.81004784019758, 0., 0., 0., 0., 0., 0., 0., 0., 116.59026982851056, -118.70608477566601, 7.225028519662741, -5.109213572507301, 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, 0., 0., 0., 0., 0., 0., 0., 0., -118.70608477566601, 121.08259875715257, -10.416247495288234, 8.03973351380169, -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, 0., 0., 0., 0., 0., 0., 0., 0., 7.225028519662741, -10.416247495288234, 71.93178675711744, -68.74056778149193, -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, 0., 0., 0., 0., 0., 0., 0., 0., -5.109213572507301, 8.03973351380169, -68.74056778149193, 65.81004784019758, 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625, 0., 0., 0., 0., -116.59026982851056, 118.70608477566601, -7.225028519662741, 5.109213572507301, 0., 0., 0., 0., 0., 0., 0., 0., 0.2855522382507991, -0.11225772018709444, -0.02094304090007367, 0.006314051366847807, 118.70608477566601, -121.08259875715257, 10.416247495288234, -8.03973351380169, 0., 0., 0., 0., 0., 0., 0., 0., -0.11225772018709444, 0.18888192025534403, 0.04846750106129806, -0.015841633186786844, -7.225028519662741, 10.416247495288234, -71.93178675711744, 68.74056778149193, 0., 0., 0., 0., 0., 0., 0., 0., -0.02094304090007367, 0.04846750106129806, 0.12951888991249833, -0.06482845236805868, 5.109213572507301, -8.03973351380169, 68.74056778149193, -65.81004784019758, 0., 0., 0., 0., 0., 0., 0., 0., 0.006314051366847807, -0.015841633186786844, -0.06482845236805868, 0.19654942257829625
        };
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                const double expected = reference[16*r+c];
                const double actual = model.H0(r,c);
                const double error = std::abs(actual - expected);
                const double tolerance = 1e-8 + 1e-10*std::abs(expected);
                if (!std::isfinite(actual) || error > tolerance) {
                    std::cerr << "Hessian mismatch: case 3, entry (" << r << "," << c << ")\n";
                    return 1;
                }
            }
        }
    }
    {
        LORd_model_params_t params{};
        params.EC = 2.e11;
        params.ES = 1.2e11;
        params.ET = 1.8e11;
        params.rC = 0.008;
        params.rS = 0.009;
        params.rT = 0.007;
        params.LC = 0.06;
        params.LS = 0.11;
        params.LT = 0.05;
        params.rhoS = 7800.;
        params.mC = 0.2;
        params.mT = 0.15;
        params.JdC = 0.0001;
        params.JdT = 0.00008;
        params.JpC = 0.0002;
        params.JpT = 0.00016;
        params.KbyC = 1.e7;
        params.KbyT = 9.e6;
        params.KbzC = 1.2e7;
        params.KbzT = 1.1e7;
        params.Omega = 1000.;
        model.make_H0(params);
        const double reference[256] = {
            -3.376317971950258e6, 5.721655870528685e6, -3.453499709110663e6, 1.108161810532236e6, 0., 0., 0., 0., 0., 0., 0., 0., 39.488502129617174, -44.170341794695034, 9.402671687387288, -4.720832022309414, 5.721655870528685e6, -2.045858555824785e7, 8.292187956066915e6, -3.555258268347749e6, 0., 0., 0., 0., 0., 0., 0., 0., -44.170341794695034, 50.43223892242247, -17.21730114729713, 10.95540401956969, -3.453499709110663e6, 8.292187956066915e6, -2.0339602739327207e7, 6.500914492370959e6, 0., 0., 0., 0., 0., 0., 0., 0., 9.402671687387288, -17.21730114729713, 64.12902232784872, -56.31439286793888, 1.108161810532236e6, -3.555258268347749e6, 6.500914492370959e6, -4.053818034555448e6, 0., 0., 0., 0., 0., 0., 0., 0., -4.720832022309414, 10.95540401956969, -56.31439286793888, 50.0798208706786, 0., 0., 0., 0., -3.376317971950258e6, 5.721655870528685e6, -3.453499709110663e6, 1.108161810532236e6, -39.488502129617174, 44.170341794695034, -9.402671687387288, 4.720832022309414, 0., 0., 0., 0., 0., 0., 0., 0., 5.721655870528685e6, -2.245858555824785e7, 8.292187956066915e6, -3.555258268347749e6, 44.170341794695034, -50.43223892242247, 17.21730114729713, -10.95540401956969, 0., 0., 0., 0., 0., 0., 0., 0., -3.453499709110663e6, 8.292187956066915e6, -2.2339602739327207e7, 6.500914492370959e6, -9.402671687387288, 17.21730114729713, -64.12902232784872, 56.31439286793888, 0., 0., 0., 0., 0., 0., 0., 0., 1.108161810532236e6, -3.555258268347749e6, 6.500914492370959e6, -4.053818034555448e6, 4.720832022309414, -10.95540401956969, 56.31439286793888, -50.0798208706786, 0., 0., 0., 0., 0., 0., 0., 0., -39.488502129617174, 44.170341794695034, -9.402671687387288, 4.720832022309414, 0.24344487065609874, -0.06341522079051264, -0.004806568754091427, -0.0014216366849242037, 0., 0., 0., 0., 0., 0., 0., 0., 44.170341794695034, -50.43223892242247, 17.21730114729713, -10.95540401956969, -0.06341522079051264, 0.16438523321465023, 0.04077436333114514, -0.0031093022570448624, 0., 0., 0., 0., 0., 0., 0., 0., -9.402671687387288, 17.21730114729713, -64.12902232784872, 56.31439286793888, -0.004806568754091427, 0.04077436333114514, 0.16740179721256154, -0.0732259293586014, 0., 0., 0., 0., 0., 0., 0., 0., 4.720832022309414, -10.95540401956969, 56.31439286793888, -50.0798208706786, -0.0014216366849242037, -0.0031093022570448624, -0.0732259293586014, 0.20351109418393165, 0., 0., 0., 0., 39.488502129617174, -44.170341794695034, 9.402671687387288, -4.720832022309414, 0., 0., 0., 0., 0., 0., 0., 0., 0.24344487065609874, -0.06341522079051264, -0.004806568754091427, -0.0014216366849242037, -44.170341794695034, 50.43223892242247, -17.21730114729713, 10.95540401956969, 0., 0., 0., 0., 0., 0., 0., 0., -0.06341522079051264, 0.16438523321465023, 0.04077436333114514, -0.0031093022570448624, 9.402671687387288, -17.21730114729713, 64.12902232784872, -56.31439286793888, 0., 0., 0., 0., 0., 0., 0., 0., -0.004806568754091427, 0.04077436333114514, 0.16740179721256154, -0.0732259293586014, -4.720832022309414, 10.95540401956969, -56.31439286793888, 50.0798208706786, 0., 0., 0., 0., 0., 0., 0., 0., -0.0014216366849242037, -0.0031093022570448624, -0.0732259293586014, 0.20351109418393165
        };
        for (int r = 0; r < 16; ++r) {
            for (int c = 0; c < 16; ++c) {
                const double expected = reference[16*r+c];
                const double actual = model.H0(r,c);
                const double error = std::abs(actual - expected);
                const double tolerance = 1e-8 + 1e-10*std::abs(expected);
                if (!std::isfinite(actual) || error > tolerance) {
                    std::cerr << "Hessian mismatch: case 4, entry (" << r << "," << c << ")\n";
                    return 1;
                }
            }
        }
    }
    std::cout << "All generated Hessian tests passed.\n";
    return 0;
}
