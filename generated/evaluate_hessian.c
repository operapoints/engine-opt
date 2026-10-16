#include "math.h"

#include "WolframRTL.h"

static WolframCompileLibrary_Functions funStructCompile;

static mint I0_0;

static mreal R0_21;

static mreal R0_27;

static mreal R0_31;

static mreal R0_40;

static mreal R0_45;

static mreal R0_61;

static mreal R0_71;

static mreal R0_74;

static mreal R0_77;

static mreal R0_80;

static mreal R0_87;

static mreal R0_90;

static mreal R0_96;

static mreal R0_98;

static mreal R0_116;

static mbool initialize = 1;

#include "evaluate_hessian.h"

DLLEXPORT int Initialize_evaluate_hessian(WolframLibraryData libData)
{
if( initialize)
{
funStructCompile = libData->compileLibraryFunctions;
R0_45 = (mreal) -0.125;
R0_77 = (mreal) 4.;
R0_116 = (mreal) 0.;
R0_87 = (mreal) -1.;
R0_40 = (mreal) 0.375;
R0_27 = (mreal) -0.375;
R0_61 = (mreal) 0.125;
R0_74 = (mreal) -0.75;
R0_31 = (mreal) -3.;
R0_96 = (mreal) 1.125;
I0_0 = (mint) 3;
R0_21 = (mreal) 3.;
R0_80 = (mreal) 0.75;
R0_71 = (mreal) -1.125;
R0_98 = (mreal) 2.25;
R0_90 = (mreal) 0.25;
initialize = 0;
}
return 0;
}

DLLEXPORT void Uninitialize_evaluate_hessian(WolframLibraryData libData)
{
if( !initialize)
{
initialize = 1;
}
}

DLLEXPORT int evaluate_hessian(WolframLibraryData libData, mreal A1, mreal A2, mreal A3, mreal A4, mreal A5, mreal A6, mreal A7, mreal A8, mreal A9, mreal A10, mreal A11, mreal A12, mreal A13, mreal A14, mreal A15, mreal A16, mreal A17, MTensor *Res)
{
mreal R0_0;
mreal R0_1;
mreal R0_2;
mreal R0_3;
mreal R0_4;
mreal R0_5;
mreal R0_6;
mreal R0_7;
mreal R0_8;
mreal R0_9;
mreal R0_10;
mreal R0_11;
mreal R0_12;
mreal R0_13;
mreal R0_14;
mreal R0_15;
mreal R0_16;
mreal R0_17;
mreal R0_18;
mreal R0_19;
mreal R0_20;
mreal R0_22;
mreal R0_23;
mreal R0_24;
mreal R0_25;
mreal R0_26;
mreal R0_28;
mreal R0_29;
mreal R0_30;
mreal R0_32;
mreal R0_33;
mreal R0_34;
mreal R0_35;
mreal R0_36;
mreal R0_37;
mreal R0_38;
mreal R0_39;
mreal R0_41;
mreal R0_42;
mreal R0_43;
mreal R0_44;
mreal R0_46;
mreal R0_47;
mreal R0_48;
mreal R0_49;
mreal R0_50;
mreal R0_51;
mreal R0_52;
mreal R0_53;
mreal R0_54;
mreal R0_55;
mreal R0_56;
mreal R0_57;
mreal R0_58;
mreal R0_59;
mreal R0_60;
mreal R0_62;
mreal R0_63;
mreal R0_64;
mreal R0_65;
mreal R0_66;
mreal R0_67;
mreal R0_68;
mreal R0_69;
mreal R0_70;
mreal R0_72;
mreal R0_73;
mreal R0_75;
mreal R0_76;
mreal R0_78;
mreal R0_79;
mreal R0_81;
mreal R0_82;
mreal R0_83;
mreal R0_84;
mreal R0_85;
mreal R0_86;
mreal R0_88;
mreal R0_89;
mreal R0_91;
mreal R0_92;
mreal R0_93;
mreal R0_94;
mreal R0_95;
mreal R0_97;
mreal R0_99;
mreal R0_100;
mreal R0_101;
mreal R0_102;
mreal R0_103;
mreal R0_104;
mreal R0_105;
mreal R0_106;
mreal R0_107;
mreal R0_108;
mreal R0_109;
mreal R0_110;
mreal R0_111;
mreal R0_112;
mreal R0_113;
mreal R0_114;
mreal R0_115;
mreal R0_117;
mreal R0_118;
mreal R0_119;
mreal R0_120;
mreal R0_121;
mreal R0_122;
MTensor* T0_0;
MTensorInitializationData Tinit;
mreal *P0;
int err = 0;
Tinit = funStructCompile->GetInitializedMTensors(libData, 1);
T0_0 = MTensorInitializationData_getTensor(Tinit, 0);
R0_0 = A1;
R0_1 = A2;
R0_2 = A3;
R0_3 = A4;
R0_4 = A5;
R0_5 = A6;
R0_6 = A7;
R0_7 = A8;
R0_8 = A9;
R0_9 = A10;
R0_10 = A11;
R0_11 = A12;
R0_12 = A13;
R0_13 = A14;
R0_14 = A15;
R0_15 = A16;
R0_16 = A17;
if( I0_0 == 0)
{
if( R0_6 == 0)
{
err = 1;
goto error_label;
}
else
{
R0_17 = 1;
}
}
else
{
mint S0 = I0_0;
mreal S1 = R0_6;
mbool S2 = 0;
if( S0 < 0)
{
S2 = 1;
S0 = -S0;
}
R0_17 = 1;
while( S0)
{
if( S0 & 1)
{
R0_17 = S1 * R0_17;
}
S1 = S1 * S1;
S0 = S0 >> 1;
}
if( S2)
{
R0_17 = 1 / R0_17;
}
}
R0_18 = 1 / R0_17;
R0_17 = 1 / R0_8;
R0_19 = R0_6 * R0_6;
R0_20 = 1 / R0_19;
R0_19 = R0_6 + R0_8;
R0_22 = R0_21 * R0_2 * R0_3 * R0_18 * R0_17 * R0_19;
R0_23 = R0_8 * R0_8;
R0_24 = 1 / R0_23;
R0_23 = R0_7 * R0_7;
R0_25 = 1 / R0_23;
R0_23 = R0_21 * R0_8;
R0_26 = R0_6 + R0_23;
R0_28 = R0_27 * R0_12 * R0_20 * R0_17 * R0_26;
R0_29 = 1 / R0_6;
R0_30 = 1 / R0_7;
R0_32 = R0_31 * R0_2 * R0_3 * R0_20 * R0_17;
R0_33 = R0_2 * R0_3 * R0_20 * R0_19;
R0_34 = R0_8 + R0_7;
R0_35 = R0_4 * R0_5 * R0_25 * R0_34;
R0_36 = R0_33 + R0_35;
R0_37 = R0_21 * R0_24 * R0_36;
if( I0_0 == 0)
{
if( R0_7 == 0)
{
err = 1;
goto error_label;
}
else
{
R0_38 = 1;
}
}
else
{
mint S0 = I0_0;
mreal S1 = R0_7;
mbool S2 = 0;
if( S0 < 0)
{
S2 = 1;
S0 = -S0;
}
R0_38 = 1;
while( S0)
{
if( S0 & 1)
{
R0_38 = S1 * R0_38;
}
S1 = S1 * S1;
S0 = S0 >> 1;
}
if( S2)
{
R0_38 = 1 / R0_38;
}
}
R0_39 = 1 / R0_38;
R0_38 = R0_40 * R0_12 * R0_29 * R0_17;
R0_41 = R0_12 * R0_26 * R0_7;
R0_42 = R0_23 + R0_7;
R0_43 = R0_13 * R0_6 * R0_42;
R0_44 = R0_41 + R0_43;
R0_46 = R0_45 * R0_29 * R0_24 * R0_30 * R0_44;
R0_47 = R0_31 * R0_4 * R0_5 * R0_17 * R0_25;
R0_48 = R0_21 * R0_4 * R0_5 * R0_17 * R0_39 * R0_34;
R0_49 = R0_40 * R0_13 * R0_17 * R0_30;
R0_50 = R0_27 * R0_13 * R0_17 * R0_25 * R0_42;
R0_51 = R0_31 * R0_2 * R0_3 * R0_18;
R0_52 = R0_19 * R0_19;
R0_53 = R0_40 * R0_12 * R0_20 * R0_17 * R0_26;
R0_54 = R0_6 * R0_6;
R0_55 = R0_13 * R0_54;
R0_56 = R0_26 * R0_26;
R0_57 = R0_12 * R0_56;
R0_58 = R0_55 + R0_57;
R0_59 = R0_34 * R0_34;
R0_60 = R0_27 * R0_12 * R0_29 * R0_17;
R0_62 = R0_61 * R0_29 * R0_24 * R0_30 * R0_44;
R0_63 = R0_7 * R0_7;
R0_64 = R0_12 * R0_63;
R0_65 = R0_42 * R0_42;
R0_66 = R0_13 * R0_65;
R0_67 = R0_64 + R0_66;
R0_68 = R0_31 * R0_4 * R0_5 * R0_39;
R0_69 = R0_27 * R0_13 * R0_17 * R0_30;
R0_70 = R0_40 * R0_13 * R0_17 * R0_25 * R0_42;
R0_72 = R0_71 * R0_12 * R0_20;
R0_73 = R0_45 * R0_20 * R0_24 * R0_58;
R0_75 = R0_74 * R0_14 * R0_20 * R0_17 * R0_26;
R0_76 = R0_8 * R0_8;
R0_78 = R0_77 * R0_16;
R0_79 = R0_45 * R0_24 * R0_25 * R0_67;
R0_81 = R0_80 * R0_14 * R0_29 * R0_17;
R0_82 = R0_21 * R0_15 * R0_6 * R0_8;
R0_83 = R0_14 + R0_78 + R0_15;
R0_84 = R0_83 * R0_6 * R0_7;
R0_85 = R0_21 * R0_14 * R0_8 * R0_7;
R0_86 = R0_82 + R0_84 + R0_85;
R0_88 = R0_87 * R0_29 * R0_24 * R0_30 * R0_86;
R0_89 = R0_88 + R0_9;
R0_91 = R0_90 * R0_89;
R0_92 = R0_76 * R0_9;
R0_93 = R0_71 * R0_13 * R0_25;
R0_94 = R0_80 * R0_15 * R0_17 * R0_30;
R0_95 = R0_74 * R0_15 * R0_17 * R0_25 * R0_42;
R0_97 = R0_96 * R0_12 * R0_20;
R0_99 = R0_98 * R0_14 * R0_20;
R0_100 = R0_99 + R0_11;
R0_101 = R0_61 * R0_20 * R0_24 * R0_58;
R0_102 = R0_14 * R0_56;
R0_103 = R0_78 + R0_15 + R0_92;
R0_104 = R0_54 * R0_103;
R0_105 = R0_102 + R0_104;
R0_106 = R0_90 * R0_20 * R0_24 * R0_105;
R0_107 = R0_61 * R0_24 * R0_25 * R0_67;
R0_108 = R0_15 * R0_65;
R0_109 = R0_14 + R0_78 + R0_92;
R0_110 = R0_63 * R0_109;
R0_111 = R0_108 + R0_110;
R0_112 = R0_90 * R0_24 * R0_25 * R0_111;
R0_113 = R0_96 * R0_13 * R0_25;
R0_114 = R0_98 * R0_15 * R0_25;
R0_115 = R0_114 + R0_10;
R0_117 = R0_87 * R0_0;
R0_118 = R0_87 * R0_2 * R0_3 * R0_18 * R0_52;
R0_119 = R0_87 * R0_4 * R0_5 * R0_30;
R0_118 = R0_118 + R0_119;
R0_119 = R0_21 * R0_24 * R0_118;
R0_117 = R0_117 + R0_119;
R0_119 = R0_87 * R0_1;
if( I0_0 == 0)
{
if( R0_7 == 0)
{
err = 1;
goto error_label;
}
else
{
R0_118 = 1;
}
}
else
{
mint S0 = I0_0;
mreal S1 = R0_7;
mbool S2 = 0;
if( S0 < 0)
{
S2 = 1;
S0 = -S0;
}
R0_118 = 1;
while( S0)
{
if( S0 & 1)
{
R0_118 = S1 * R0_118;
}
S1 = S1 * S1;
S0 = S0 >> 1;
}
if( S2)
{
R0_118 = 1 / R0_118;
}
}
R0_120 = R0_2 * R0_3 * R0_118;
R0_118 = R0_4 * R0_5 * R0_6 * R0_59;
R0_120 = R0_120 + R0_118;
R0_118 = R0_31 * R0_29 * R0_24 * R0_39 * R0_120;
R0_119 = R0_119 + R0_118;
if( I0_0 == 0)
{
if( R0_6 == 0)
{
err = 1;
goto error_label;
}
else
{
R0_118 = 1;
}
}
else
{
mint S0 = I0_0;
mreal S1 = R0_6;
mbool S2 = 0;
if( S0 < 0)
{
S2 = 1;
S0 = -S0;
}
R0_118 = 1;
while( S0)
{
if( S0 & 1)
{
R0_118 = S1 * R0_118;
}
S1 = S1 * S1;
S0 = S0 >> 1;
}
if( S2)
{
R0_118 = 1 / R0_118;
}
}
R0_120 = R0_4 * R0_5 * R0_118;
R0_118 = R0_2 * R0_3 * R0_52 * R0_7;
R0_120 = R0_120 + R0_118;
R0_118 = R0_31 * R0_18 * R0_24 * R0_30 * R0_120;
R0_120 = R0_21 * R0_2 * R0_3 * R0_29;
R0_121 = R0_0 + R0_1;
R0_121 = R0_121 * R0_76;
R0_122 = R0_21 * R0_4 * R0_5 * R0_39 * R0_59;
R0_120 = R0_120 + R0_121 + R0_122;
R0_121 = R0_87 * R0_24 * R0_120;
{
mint S0[2] = {16, 16};
err = funStructCompile->MTensor_allocate(T0_0, 3, 2, S0);
if( err)
{
goto error_label;
}
P0 = MTensor_getRealDataMacro(*T0_0);
P0[0] = R0_51;
P0[1] = R0_22;
P0[2] = R0_32;
P0[3] = R0_116;
P0[4] = R0_116;
P0[5] = R0_116;
P0[6] = R0_116;
P0[7] = R0_116;
P0[8] = R0_116;
P0[9] = R0_116;
P0[10] = R0_116;
P0[11] = R0_116;
P0[12] = R0_97;
P0[13] = R0_28;
P0[14] = R0_38;
P0[15] = R0_116;
P0[16] = R0_22;
P0[17] = R0_117;
P0[18] = R0_37;
P0[19] = R0_47;
P0[20] = R0_116;
P0[21] = R0_116;
P0[22] = R0_116;
P0[23] = R0_116;
P0[24] = R0_116;
P0[25] = R0_116;
P0[26] = R0_116;
P0[27] = R0_116;
P0[28] = R0_28;
P0[29] = R0_101;
P0[30] = R0_46;
P0[31] = R0_49;
P0[32] = R0_32;
P0[33] = R0_37;
P0[34] = R0_119;
P0[35] = R0_48;
P0[36] = R0_116;
P0[37] = R0_116;
P0[38] = R0_116;
P0[39] = R0_116;
P0[40] = R0_116;
P0[41] = R0_116;
P0[42] = R0_116;
P0[43] = R0_116;
P0[44] = R0_38;
P0[45] = R0_46;
P0[46] = R0_107;
P0[47] = R0_50;
P0[48] = R0_116;
P0[49] = R0_47;
P0[50] = R0_48;
P0[51] = R0_68;
P0[52] = R0_116;
P0[53] = R0_116;
P0[54] = R0_116;
P0[55] = R0_116;
P0[56] = R0_116;
P0[57] = R0_116;
P0[58] = R0_116;
P0[59] = R0_116;
P0[60] = R0_116;
P0[61] = R0_49;
P0[62] = R0_50;
P0[63] = R0_113;
P0[64] = R0_116;
P0[65] = R0_116;
P0[66] = R0_116;
P0[67] = R0_116;
P0[68] = R0_51;
P0[69] = R0_22;
P0[70] = R0_32;
P0[71] = R0_116;
P0[72] = R0_72;
P0[73] = R0_53;
P0[74] = R0_60;
P0[75] = R0_116;
P0[76] = R0_116;
P0[77] = R0_116;
P0[78] = R0_116;
P0[79] = R0_116;
P0[80] = R0_116;
P0[81] = R0_116;
P0[82] = R0_116;
P0[83] = R0_116;
P0[84] = R0_22;
P0[85] = R0_118;
P0[86] = R0_37;
P0[87] = R0_47;
P0[88] = R0_53;
P0[89] = R0_73;
P0[90] = R0_62;
P0[91] = R0_69;
P0[92] = R0_116;
P0[93] = R0_116;
P0[94] = R0_116;
P0[95] = R0_116;
P0[96] = R0_116;
P0[97] = R0_116;
P0[98] = R0_116;
P0[99] = R0_116;
P0[100] = R0_32;
P0[101] = R0_37;
P0[102] = R0_121;
P0[103] = R0_48;
P0[104] = R0_60;
P0[105] = R0_62;
P0[106] = R0_79;
P0[107] = R0_70;
P0[108] = R0_116;
P0[109] = R0_116;
P0[110] = R0_116;
P0[111] = R0_116;
P0[112] = R0_116;
P0[113] = R0_116;
P0[114] = R0_116;
P0[115] = R0_116;
P0[116] = R0_116;
P0[117] = R0_47;
P0[118] = R0_48;
P0[119] = R0_68;
P0[120] = R0_116;
P0[121] = R0_69;
P0[122] = R0_70;
P0[123] = R0_93;
P0[124] = R0_116;
P0[125] = R0_116;
P0[126] = R0_116;
P0[127] = R0_116;
P0[128] = R0_116;
P0[129] = R0_116;
P0[130] = R0_116;
P0[131] = R0_116;
P0[132] = R0_72;
P0[133] = R0_53;
P0[134] = R0_60;
P0[135] = R0_116;
P0[136] = R0_100;
P0[137] = R0_75;
P0[138] = R0_81;
P0[139] = R0_116;
P0[140] = R0_116;
P0[141] = R0_116;
P0[142] = R0_116;
P0[143] = R0_116;
P0[144] = R0_116;
P0[145] = R0_116;
P0[146] = R0_116;
P0[147] = R0_116;
P0[148] = R0_53;
P0[149] = R0_73;
P0[150] = R0_62;
P0[151] = R0_69;
P0[152] = R0_75;
P0[153] = R0_106;
P0[154] = R0_91;
P0[155] = R0_94;
P0[156] = R0_116;
P0[157] = R0_116;
P0[158] = R0_116;
P0[159] = R0_116;
P0[160] = R0_116;
P0[161] = R0_116;
P0[162] = R0_116;
P0[163] = R0_116;
P0[164] = R0_60;
P0[165] = R0_62;
P0[166] = R0_79;
P0[167] = R0_70;
P0[168] = R0_81;
P0[169] = R0_91;
P0[170] = R0_112;
P0[171] = R0_95;
P0[172] = R0_116;
P0[173] = R0_116;
P0[174] = R0_116;
P0[175] = R0_116;
P0[176] = R0_116;
P0[177] = R0_116;
P0[178] = R0_116;
P0[179] = R0_116;
P0[180] = R0_116;
P0[181] = R0_69;
P0[182] = R0_70;
P0[183] = R0_93;
P0[184] = R0_116;
P0[185] = R0_94;
P0[186] = R0_95;
P0[187] = R0_115;
P0[188] = R0_116;
P0[189] = R0_116;
P0[190] = R0_116;
P0[191] = R0_116;
P0[192] = R0_97;
P0[193] = R0_28;
P0[194] = R0_38;
P0[195] = R0_116;
P0[196] = R0_116;
P0[197] = R0_116;
P0[198] = R0_116;
P0[199] = R0_116;
P0[200] = R0_116;
P0[201] = R0_116;
P0[202] = R0_116;
P0[203] = R0_116;
P0[204] = R0_100;
P0[205] = R0_75;
P0[206] = R0_81;
P0[207] = R0_116;
P0[208] = R0_28;
P0[209] = R0_101;
P0[210] = R0_46;
P0[211] = R0_49;
P0[212] = R0_116;
P0[213] = R0_116;
P0[214] = R0_116;
P0[215] = R0_116;
P0[216] = R0_116;
P0[217] = R0_116;
P0[218] = R0_116;
P0[219] = R0_116;
P0[220] = R0_75;
P0[221] = R0_106;
P0[222] = R0_91;
P0[223] = R0_94;
P0[224] = R0_38;
P0[225] = R0_46;
P0[226] = R0_107;
P0[227] = R0_50;
P0[228] = R0_116;
P0[229] = R0_116;
P0[230] = R0_116;
P0[231] = R0_116;
P0[232] = R0_116;
P0[233] = R0_116;
P0[234] = R0_116;
P0[235] = R0_116;
P0[236] = R0_81;
P0[237] = R0_91;
P0[238] = R0_112;
P0[239] = R0_95;
P0[240] = R0_116;
P0[241] = R0_49;
P0[242] = R0_50;
P0[243] = R0_113;
P0[244] = R0_116;
P0[245] = R0_116;
P0[246] = R0_116;
P0[247] = R0_116;
P0[248] = R0_116;
P0[249] = R0_116;
P0[250] = R0_116;
P0[251] = R0_116;
P0[252] = R0_116;
P0[253] = R0_94;
P0[254] = R0_95;
P0[255] = R0_115;
}
funStructCompile->MTensor_copy(Res, *T0_0);
error_label:
funStructCompile->ReleaseInitializedMTensors(Tinit);
funStructCompile->WolframLibraryData_cleanUp(libData, 1);
return err;
}

