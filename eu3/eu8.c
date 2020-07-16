
double Betta_W_table4_1(const char *Category)
{
    if(strcmp(Category,"C235") == 0)
        return 0.8;
    else
        {print_name("Error in Betta_W_table4_1");FATAL();}
    return 1.e30;
}

enum CODE_NORM{Code_SNiP,Code_Eurocode} CODE;
void SNiP(const char *Category, double thk, double &PY, double &PU)
{
    if(CODE == Code_SNiP)
    {
        if(strcmp(Category,"C235") == 0)//(Category == "C235")
        {
            if(2e-3 <= thk && thk <=20e-3)
            {
                PY = 230e6;
                PU = 350e6;
            }
            else if(20e-3 <= thk && thk <=40e-3)
            {
                PY = 220e6;
                PU = 350e6;
            }
            else if(40e-3 <= thk && thk <=100e-3)
            {
                PY = 210e6;
                PU = 350e6;
            }
            else if(100e-3 <= thk)
            {
                PY = 190e6;
                PU = 350e6;
            }
            else
            {
                print_name("You can`t use C235");
                printf("Category = %s\n",Category);
                printf("thk = %f\n",thk);
                FATAL();
            }
        }
        if(strcmp(Category,"C245") == 0)//(Category == "C245")
        {
            if(2e-3 <= thk && thk <=20e-3)
            {
                PY = 240e6;
                PU = 360e6;
            }
            else
            {
                print_name("You can`t use C245");
                FATAL();
            }
        }
    }
    else if(CODE == Code_Eurocode)
    {
        if(strcmp(Category,"C235") == 0)//(Category == "C235")
        {
            if(thk <=40e-3)
            {
                PY = 235e6;
                PU = 360e6;
            }
            else
            {
                PY = 215e6;
                PU = 340e6;
            }
        }
        else if(strcmp(Category,"C275") == 0)//(Category == "C275")
        {
            if(thk <=40e-3)
            {
                PY = 275e6;
                PU = 430e6;
            }
            else
            {
                PY = 255e6;
                PU = 410e6;
            }
        }
        else
        {
            print_name("You can`t use this with Eurocode");
            FATAL();
        }
    }
    else
    {
        print_name("You can`t use this: What is this code");
        FATAL();
    }
};

///
/// UNITS
///

//    MPA   -   %5.1f MPa
//    MM    -   %5.1f mm
//    MM2   -   %5.1f mm2

///
/// CHOOSE STANDART
///

enum STANDART {STANDART_SNiP,STANDART_Eurocode};

///
/// FUNCTIONS FOR HOLE
///

enum HOLE_TYPE {NORMAL, BIG_HOLE, STUPID};

HOLE_TYPE EN1090_2_TABLE_11 (double DiameterBolts, double DiameterOfHole)
{
    HOLE_TYPE HT;
    //printf ("Diameter Of Bolt = %3.1 mm\n", DiameterBolts  *0.001);
    //printf ("Diameter Of Hole = %3.1 mm\n", DiameterOfHole *0.001);
    if(DiameterBolts<0.015)
    {
             if(DiameterOfHole-DiameterBolts<=0.001001)
                HT = NORMAL;
        else if(DiameterOfHole-DiameterBolts<=0.003001)
                HT = BIG_HOLE;
        else
                HT = STUPID;
    }
    else if(DiameterBolts<0.023)
    {
             if(DiameterOfHole-DiameterBolts<=0.002001)
                HT = NORMAL;
        else if(DiameterOfHole-DiameterBolts<=0.004001)
                HT = BIG_HOLE;
        else
                HT = STUPID;
    }
    else if(DiameterBolts<0.026)
    {
             if(DiameterOfHole-DiameterBolts<=0.002001)
                HT = NORMAL;
        else if(DiameterOfHole-DiameterBolts<=0.006001)
                HT = BIG_HOLE;
        else
                HT = STUPID;
    }
    else
    {
             if(DiameterOfHole-DiameterBolts<=0.003001)
                HT = NORMAL;
        else if(DiameterOfHole-DiameterBolts<=0.008001)
                HT = BIG_HOLE;
        else
                HT = STUPID;
    }
    if(HT == NORMAL)
    {
        printf("Hole is Normal round holes by table 11 EN1090-2/n");
    }
    else if(HT == BIG_HOLE)
    {
         printf("Hole is Oversize round holes by table 11 EN1090-2/n");
    }
    else
    {
        print_name("HOLE IS STUPID");
        FATAL();
    }
    return HT;
}

///
/// TENSION RESISTANCE FOR INDIVIDUAL BOLT
///

void EN1993_1_8_TABLE_3_4_FtRd(double &FtRd, double Pub, double As, double gamma_M2, bool OUT)
{
    double k2 = 0.9;
         FtRd = k2 * Pub * As / gamma_M2;
    if(OUT)
    {
        printf("Tension resistance FtRd:\n");
        printf("FtRd = k2 * Pub * As / gamma_M2\n");
        printf("FtRd = %.1f * %.1f * %.1f / %.2f\n",k2,Pub*1e-6,As*1e-6,gamma_M2);
        printf("FtRd = %.1fkN\n",FtRd*1e-3);
    }
}

///
/// SHEAR RESISTANCE FOR INDIVIDUAL BOLT
///

void EN1993_1_8_TABLE_3_4_FvRd(double &FvRd, double Pub, double A, double gamma_M2, BOLT_CLASS BS, bool OUT)
{
    double alphaV = 0.0;
    switch (BS)
    {
        case g3_6:
        case g_Ct3:
        case g_09G2S:
        case g_10G2S1:
            alphaV = 0.6;
            break;
        case g5_6:
            alphaV = 0.6;
            break;
        case g8_8:
            alphaV = 0.6;
            break;
        case g10_9:
            alphaV = 0.5;
            break;
        default:
            print_name("ERROR: EN1993_1_8_TABLE_3_4_FvRd ");
            FATAL();
    }
    FvRd = alphaV * Pub * A / gamma_M2;
    if(OUT)
    {
        printf("Shear resistance FvRd:\n");
        printf("alphaV = %1.1f\n",alphaV);
        printf("FvRd = alphaV * Pub * A / gamma_M2\n");
        printf("FvRd = %.1f * %.1f * %.1f / %.2f\n",alphaV,Pub*1e-6,A*1e-6,gamma_M2);
        printf("FvRd = %.1fkN\n",FvRd*1e-3);
    }
}

///
/// PUNCHING SHEAR RESISTANCE FOR INDIVIDUAL BOLT
///
void EN1993_1_8_TABLE_3_4_BpRd(double &BpRd, double DiameterBolt, double tp, double Fu, double gamma_M2, bool OUT)
{
    BpRd = 0.6 * CONST_M_PI * BOLT_Dm(DiameterBolt) * tp * Fu / gamma_M2;
    if(OUT)
    {
        printf("Punching shear resistance BpRd:\n");
        printf("BpRd = 0.6 * PI * Dm * tp * Fu / gamma_M2\n");
        printf("BpRd = 0.6 * PI * %.1f * %.1f * %.1f / %.1f\n",BOLT_Dm(DiameterBolt)*1e3,tp*1e3,Fu*1e-6,gamma_M2);
        printf("BpRd = %.1fkN\n",BpRd*1e-3);
    }
}

///
/// MATERIAL OF STEEL
///

enum STEEL {C235,C245,C275} ;

void STEEL_STRESS( STEEL Steel, double thk, double &Py, double &Pu, STANDART STD, bool OUT)
{
    if(STD == STANDART_SNiP)
    {
        switch(Steel)
        {
            case C235:
                if(2e-3 <= thk && thk <=20e-3)
                { Py = 230e6; Pu = 350e6; }
                else if(20e-3 <= thk && thk <=40e-3)
                { Py = 220e6; Pu = 350e6; }
                else if(40e-3 <= thk && thk <=100e-3)
                { Py = 210e6; Pu = 350e6; }
                else if(100e-3 <= thk)
                { Py = 190e6; Pu = 350e6; }
                else
                {
                    print_name("You can`t use C235");
                    FATAL();
                }
                break;
            default:
                print_name("ERROR: STEEL_STRESS");
                FATAL();
        }
    }
    else if(STD == STANDART_SNiP)
    {
        switch(Steel)
        {
            case C235:
                if(thk <=40e-3)
                { Py = 235e6; Pu = 360e6; }
                else
                { Py = 215e6; Pu = 340e6; }
            default:
                print_name("ERROR: STEEL_STRESS");
                FATAL();
        }
    }
    else
        FATAL();
    if(OUT)
    {
        printf("Py = %5.1 MPa\n",Py*1e-6);
        printf("Pu = %5.1 MPa\n",Pu*1e-6);
    }
}

///
/// BOLT CONNECTION: COEFFICIENT FROM PICTURE 6.1
///
double EN1993_1_8_Picture_6_1_APLHA(double lambda1, double lambda2, bool OUT=true)
{
    bool DEBUG =false;// true;//
    if(DEBUG) printf("lambda1=%.2f,lambda2=%.2f\n",lambda1,lambda2);
    if(lambda1<0 || lambda2<0)
    {
        line();
        printf("Stupid check\n");
        printf("lambda1=%+.2f,lambda2=%+.2f\n",lambda1,lambda2);
        line();
    }
    double ALPHA=1e30;
//    double value1 = 1e30;
//    double
//       F1 = + 0.99477448
//            - 2.45848503*lambda2
//            + 3.15497168*pow(lambda2,2.)
//            - 2.23017434*pow(lambda2,3.)
//            + 0.52850212*pow(lambda2,4.);
//    double
//       F2 = + 1.04213142
//            - 0.85759182*lambda2
//            + 1.15828063*pow(lambda2,2.)
//            - 0.79910192*pow(lambda2,3.)
//            + 0.21398139*pow(lambda2,4.);
//    double
//       F3 = + 8.130283
//            + 4.488295 * lambda1
//            - 3.441231 * lambda2
//            - 16.699661* pow(lambda1,2.)
//            + 4.657641 * pow(lambda2,2.)
//            - 6.802532 * lambda1*lambda2
//            + 8.747474 * pow(lambda1,3.)
//            - 1.197675 * pow(lambda2,3.)
//            - 1.227359 * lambda1*pow(lambda2,2.)
//            + 8.318217 * pow(lambda1,2.)*lambda2
//            + 0.000000 * pow(lambda1,4.)
//            + 0.000000 * pow(lambda2,4.);
//    double
//       F4 = + 1.245666
//            + 39.333003* lambda1
//            - 3.580332 * lambda2
//            - 55.940605* pow(lambda1,2.)
//            + 40.544586* pow(lambda2,2.)
//            - 55.343570* lambda1*lambda2
//            + 21.049463* pow(lambda1,3.)
//            - 33.001768* pow(lambda2,3.)
//            + 2.792410 * lambda1*pow(lambda2,2.)
//            + 44.062493* pow(lambda1,2.)*lambda2
//            + 0.000000 * pow(lambda1,4.)
//            + 0.000000 * pow(lambda2,4.);
//    double
//       F5 = - 86.505200
//            +478.588870* lambda1
//            + 79.430092* lambda2
//            -935.102794* pow(lambda1,2.)
//            -329.854733* pow(lambda2,2.)
//            - 68.228567* lambda1*lambda2
//            +809.056164* pow(lambda1,3.)
//            +531.672952* pow(lambda2,3.)
//            +252.193252* lambda1*pow(lambda2,2.)
//            - 44.242644* pow(lambda1,2.)*lambda2
//            -254.659837* pow(lambda1,4.)
//            -605.622885* pow(lambda2,4.);
//    double
//       F6 = - 226.979097
//            +1095.760732* lambda1
//            - 12.1186777* lambda2
//            -1818.467314* pow(lambda1,2.)
//            + 717.104423* pow(lambda2,2.)
//            - 264.307024* lambda1*lambda2
//            +1369.007748* pow(lambda1,3.)
//            -2120.516058* pow(lambda2,3.)
//            -  69.105002* lambda1*pow(lambda2,2.)
//            + 195.697905* pow(lambda1,2.)*lambda2
//            - 381.685783* pow(lambda1,4.)
//            +2562.146768* pow(lambda2,4.);
//
//    if(lambda1 <= F1) value1 = 2* CONST_M_PI;
////    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
//    if(lambda1 >= F2) value1 =  4.45 ;
////    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
//    if(F1 < lambda1 && lambda1 < F2)
//    {
//        if(lambda2 >= 0.45)
//            value1 = min(F3,2*CONST_M_PI);
////        printf("ALPHA3   = %.4e\n",ALPHA  );//DEBUG
//        if((0.2768*lambda1+0.14) <= lambda2 && lambda2 < 0.45)
//            value1 = min(F4,2*CONST_M_PI);
////        printf("ALPHA4   = %.4e\n",ALPHA  );//DEBUG
//        if((1.2971*lambda1-0.7782) <= lambda2 && lambda2 < (0.2768*lambda1+0.14))
//            value1 = min(F5,2*CONST_M_PI);
////        printf("ALPHA5   = %.4e\n",ALPHA  );//DEBUG
//        if(lambda2 < (1.2971*lambda1-0.7782))
//            value1 = min(F6,2*CONST_M_PI);
////        printf("ALPHA6   = %.4e\n",ALPHA  );//DEBUG
//    }
//    printf("value1 = %f\n",value1);


    /////////////////////////////////
    /// TABLE
    /////////////////////////////////
    double TABLE[9][5] =
    {
        {8       ,0.9625419,-3.6454601,5.0308278,3.0586037},
        {7       ,0.9636739,-2.6572564,3.9919931,2.9602062},
        {2*CONST_M_PI  ,0.9726341,-1.9932994,3.2341635,2.7782066},
        {6       ,0.9496498,-1.4803090,2.6267325,3.1078904},
        {5.5     ,0.9785034,-1.1933579,2.3203117,2.9238812},
        {5       ,0.9637229,-0.7743167,1.9115241,2.9607118},
        {4.75    ,0.9946669,-0.7516043,1.9434759,1.9056765},
        {4.5     ,1.2055126,-3.8888717,7.1149599,0.9088292},
        {4.45    ,1.2265766,-3.2360189,6.2854489,1.0927973}
    };
//    printf("DEBUG 05\n");//DEBUG
    double LAMBDA1[9];
    for(type_LLU i=0;i<9;i++)
    {
        double K = TABLE[i][1];
        double A = TABLE[i][2];
        double B = TABLE[i][3];
        double C = TABLE[i][4];
        LAMBDA1[i] = K + (A*lambda2)/pow(1+pow(B*lambda2,C),1./C);
        if(DEBUG)printf("[%+.5e, %+.5e, %+.5e, %+.5e] = %+.5e\n",K,A,B,C,LAMBDA1[i]);
    }
    {
//        double x[9], y[9];
//        for(int i=0;i<9;i++)
//        {
//            x[i] = TABLE[0][5*i];
//            y[i] = LAMBDA1[i];
//        }
//        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
//        gsl_spline *spline   ;
//        spline = ::gsl_spline_alloc (::gsl_interp_cspline/*linear*/, 9);
//        gsl_spline_init (spline, y, x, 9);
//        ALPHA = gsl_spline_eval (spline, lambda1, acc);
//        gsl_spline_free (spline);
//        gsl_interp_accel_free (acc);
        //Graph(Array<double> x, Array<double> y, double pX, double &pY)
        Array<double> x,y;
        for(int i=0;i<9;i++)
        {
            x.Add(TABLE[0][5*i]);
            y.Add(LAMBDA1[i]);
        }
        Graph(y,x,lambda1,ALPHA);
        x.Delete();
        y.Delete();
    }

//    printf("F1      = %.4e\n",F1);     //DEBUG
//    printf("F2      = %.4e\n",F2);     //DEBUG
//    printf("F3      = %.4e\n",F3);     //DEBUG
//    printf("lambda1 = %.4e\n",lambda1);//DEBUG
//    printf("lambda2 = %.4e\n",lambda2);//DEBUG
//    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
    if(OUT || DEBUG)
    {
        printf("Alpha(lambda1=%.2f,lambda2=%.2f) = %.2f by Figure 6.11 EN1993-1-8\n",lambda1,lambda2,ALPHA);
    }
    return ALPHA;
}

void   EN1993_1_8_Formula3_2(double Fu, double D, double Thk, double gammaM2,
                             double &Fb_Rd, bool OUT = true)
{
    ///
    /// Design bearing resistance
    ///
    Fb_Rd = 1.5*Fu*D*Thk/gammaM2;
    if(OUT)
    {
        printf("Design bearing resistance by Formula(3.2) EN1993-1-8\n");
        printf("Fb_Rd = 1.5*Fu*d*t/gammaM2\n");
        printf("Fb_Rd = 1.5*%.1f*%.1f*%.1f/%.2f\n",
                Fu*1e-6, D*1e3,Thk*1e3, gammaM2);
        printf("Fb_Rd = %.1f kN\n",Fb_Rd*1e-3);
    }
}

void   EN1993_1_8_Formula3_9(double Fu, double Ant, double gammaM2,
                             double Fy, double Anv, double gammaM0,
                             double &Veff_t_Rd, bool OUT = true)
{
    ///
    /// Block Tearing Verification
    ///
    Veff_t_Rd = Fu*Ant/gammaM2+(1/sqrt(3))*Fy*Anv/gammaM0;
    if(OUT)
    {
        printf("Block Tearing Verification\n");
        printf("Ant = %.1f mm2\n",Ant*1e6);
        printf("Anv = %.1f mm2\n",Anv*1e6);
        printf("Veff = Fu*Ant/gamma_M2 + (1/sqrt(3))*Fy*Anv/gamma_M0\n");
        printf("Veff = %.1f *%.1f /%.2f + 1/%.2f *%.1f *%.1f /%.2f \n",
                Fu*1e-6, Ant*1e+6, gammaM2,
                sqrt(3.),
                Fy*1e-6, Anv*1e+6, gammaM0);
        printf("Veff_Rd = %.1f kN\n",Veff_t_Rd*1e-3);
    }
}

double EN1993_1_8_Formula4_1(double Sigma_perp, double Tau_perp, double Tau_parr,
                             double Fu,const char *Category, double gamma_M2, bool OUT = true)
{
    double factor = 0;
    ///
    /// Weld Verification
    ///
    double SS       = sqrt(pow(Sigma_perp,2.)+3*(pow(Tau_perp,2.)+pow(Tau_parr,2.)));
    double SS_left  = Fu/Betta_W_table4_1(Category)/gamma_M2;
    double SS_left2 = 0.9*Fu/gamma_M2;
    if(Fu < 1e6 || gamma_M2 < 0.99) {print_name("STUPID in EN1993_1_8_Formula4_1");FATAL();}
    if(OUT)
    {
        printf("Weld Verification\n");
        printf("Sigma_perp = %.2f MPa\n",Sigma_perp*1e-6);
        printf("Tau_perp   = %.2f MPa\n",Tau_perp  *1e-6);
        printf("Tau_parr   = %.2f MPa\n",Tau_parr  *1e-6);
        printf("sqrt(Sigma_perp^2+3*(Tau_perp^2+Tau_parr^2)) <= Fu/Betta_W_table4_1/gamma_M2\n");
        printf("sqrt(%.2f^2+3*(%.2f^2+%.2f^2)) <= %.2f/%.2f/%.2f\n",Sigma_perp*1e-6,Tau_perp*1e-6,Tau_parr*1e-6,Fu*1e-6,Betta_W_table4_1(Category),gamma_M2);
        printf("%.2f MPa<= %.2f MPa\n",SS*1e-6,SS_left*1e-6);
    }
    Printf_CALC(SS/SS_left,factor,OUT);
    if(OUT)
    {
        printf("Sigma_perp <= 0.9*Fu/gamma_M2\n");
        printf("%.2f <= 0.9*%.2f/%.2f\n",Sigma_perp*1e-6,Fu*1e-6,gamma_M2);
        printf("%.2f MPa <= %.2f MPa\n",Sigma_perp*1e-6,SS_left2*1e-6);
    }
    Printf_CALC(Sigma_perp/SS_left2,factor,OUT);
    return factor;
}



///
/// TABLE 6.2 EN1993-1-8
///

double EN1993_1_8_Table6_2_M_pl_1_RD(double l_eff_1, double tf, double Fy, double gamma_M0, bool OUT = true)
{
    double M_pl_1_RD = 0.25*l_eff_1*pow(tf,2.)*Fy/gamma_M0;
    if(OUT)
    {
        printf("M_pl_1_Rd = 0.25*l_eff_1*tf^2*Fy/gamma_M0\n");
        printf("M_pl_1_Rd = 0.25*%.1f*%.1f^2*%.1f/%.2f\n",l_eff_1*1e3,tf*1e3,Fy*1e-6,gamma_M0);
        printf("M_pl_1_Rd = %.1f kN*m\n",M_pl_1_RD*1e-3);
    }
    return M_pl_1_RD;
};

double EN1993_1_8_Table6_2_M_pl_2_RD(double l_eff_2, double tf, double Fy, double gamma_M0, bool OUT = true)
{
    double M_pl_2_RD = 0.25*l_eff_2*pow(tf,2.)*Fy/gamma_M0;
    if(OUT)
    {
        printf("M_pl_2_RD = 0.25*l_eff_2*tf^2*Fy/gamma_M0\n");
        printf("M_pl_2_RD = 0.25*%.1f*%.1f^2*%.1f/%.2f\n",l_eff_2*1e3,tf*1e3,Fy*1e-6,gamma_M0);
        printf("M_pl_2_RD = %.1f kN*m\n",M_pl_2_RD*1e-3);
    }
    return M_pl_2_RD;
};

double EN1993_1_8_Table6_2_M_bp_RD(double l_eff_1, double tbp, double Fy_bp, double gamma_M0, bool OUT = true)
{
    double M_bp_RD = 0.25*l_eff_1*pow(tbp,2.)*Fy_bp/gamma_M0;
    if(OUT)
    {
        printf("M_bp_RD = 0.25*l_eff_1*tbp^2*Fy_bp/gamma_M0\n");
        printf("M_bp_RD = 0.25*%.1f*%.1f^2*%.1f/%.2f\n",l_eff_1*1e3,tbp*1e3,Fy_bp*1e-6,gamma_M0);
        printf("M_bp_RD = %.1f kN*m\n",M_bp_RD*1e-3);
    }
    return M_bp_RD;
};

double EN1993_1_8_Table6_2_F_t1_Rd(double M_pl_1_RD, double m, bool OUT = true)
{
    double F_t1_Rd = 4.*M_pl_1_RD/m;
    if(OUT)
    {
        printf("F_t1_Rd = 4*M_pl_1_RD/m\n");
        printf("F_t1_Rd = 4*%.1f/%.1f\n",M_pl_1_RD*1e-3,m*1e3);
        printf("F_t1_Rd = %.1f kN\n",F_t1_Rd*1e-3);
    }
    return F_t1_Rd;
};

double EN1993_1_8_Table6_2_F_t2_Rd(double M_pl_2_RD, double n, double F_t_Rd, double m, type_LLU number_bolts, bool OUT = true)
{
    double F_t2_Rd = (2*M_pl_2_RD+n*number_bolts*F_t_Rd)/(m+n);
    if(OUT)
    {
        printf("F_t2_Rd = (2*M_pl_2_RD+n*number_bolts*F_t_Rd)/(m+n)\n");
        printf("F_t2_Rd = (2*%.1f+%.1f*%u*%.1f)/(%.1f+%.1f)\n",M_pl_2_RD*1e-3,n*1e3,number_bolts,F_t_Rd*1e-3,m*1e3,n*1e3);
        printf("F_t2_Rd = %.1f kN\n",F_t2_Rd*1e-3);
    }
    return F_t2_Rd;
};

double EN1993_1_8_Table6_2_F_t3_Rd(double F_t_Rd, type_LLU number_bolts, bool OUT = true)
{
    double F_t3_Rd = number_bolts*F_t_Rd;
    if(OUT)
    {
        printf("F_t3_Rd = number_bolts*F_t_Rd\n");
        printf("F_t3_Rd = %u*%.1f\n",number_bolts,F_t_Rd*1e-3);
        printf("F_t3_Rd = %.1f kN\n",F_t3_Rd*1e-3);
    }
    return F_t3_Rd;
};

double EN1993_1_8_Table6_2_F_t_Rd(double F_t1_Rd,double F_t2_Rd,double F_t3_Rd, bool OUT = true)
{
    double F_t_Rd = min(F_t1_Rd,min(F_t2_Rd,F_t3_Rd));
    if(OUT)
    {
        printf("F_t_Rd  = min(F_t1_Rd,F_t2_Rd,F_t3_Rd)\n");
        printf("F_t_Rd  = min(%.1f,%.1f,%.1f)\n",F_t1_Rd*1e-3,F_t2_Rd*1e-3,F_t3_Rd*1e-3);
        printf("F_t_Rd  = %.1f kN\n",F_t_Rd*1e-3);
    }
    return F_t_Rd;
};



double EN1993_1_8_Formula6_22_F_t_wb_Rd(double Beff_t_wb,double twb,double Fy_wb, double gamma_M0, bool OUT = true)
{
    double F_t_wb_Rd = Beff_t_wb * twb * Fy_wb/gamma_M0;
    if(OUT)
    {
        printf("F_t_wb_Rd = Beff_t_wb * twb * Fy_wb/gamma_M0, see formula(6.22) EN1993-1-8\n");
        printf("F_t_wb_Rd = %.1f * %.1f * %.1f/%.2f\n",Beff_t_wb*1e3,twb*1e3,Fy_wb*1e-6,gamma_M0);
        printf("F_t_wb_Rd = %.1f kN\n",F_t_wb_Rd*1e-3);
    }
    return F_t_wb_Rd;
}

///
/// Table 6.6
///

void EN1993_1_8_Table_6_6_Leff_group(double &Leff1,double &Leff2, double Paramaters[4], type_LLU Row_table, bool OUT = true)
{
    double m      = Paramaters[0];
    double p      = Paramaters[1];
    double _alpha = Paramaters[2];
    double e      = Paramaters[3];
    for(type_LLU i=0;i<4;i++)
    {
        if(Paramaters[i] < 0)
        {
            printf("Stupid check in EN1993_1_8_Table_6_6_Leff_group\n");
            for(type_LLU j=0;j<4;j++)
                printf("Paramaters[%u] = %.1f mm",j,Paramaters[j]*1e3);
            FATAL();
        }
    }
    double Leff_cp = 0.0;
    double Leff_nc = 0.0;
    switch(Row_table)
    {
        case(1):
        {
            printf("Bolt-row location: Bolt-row outside tension flange of beam\n");
            FATAL();
            break;
        }
        case(2):
        {
            Leff_cp = CONST_M_PI*m+p;
            Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e);
            if(OUT)
            {
                printf("Bolt-row location: First bolt-row below tension flange of beam\n");
                printf(" Leff_cp = PI*m+p\n");
                printf(" Leff_cp = PI*%.1f+%.1f\n",m*1e3,p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e)\n");
                printf(" Leff_nc = 0.5*%.1f+%.2f*%.1f-(2*%.1f+0.625*%.1f)\n",p*1e3,_alpha,m*1e3,m*1e3,e*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        case(3):
        {
            Leff_cp = 2*p;
            Leff_nc = p;
            if(OUT)
            {
                printf("Bolt-row location: Other inner bolt-row\n");
                printf(" Leff_cp = 2*p\n");
                printf(" Leff_cp = 2*%.1f\n",p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = p\n");
                printf(" Leff_nc = %.1f\n",p*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        case(4):
        {
            Leff_cp = CONST_M_PI*m+p;
            Leff_nc = 2*m+0.625*e+0.5*p;
            if(OUT)
            {
                printf("Bolt-row location: Other end bolt-row\n");
                printf(" Leff_cp = PI*m+p\n");
                printf(" Leff_cp = PI*%.1f+%.1f\n",m*1e3,p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = 2*m+0.625*e+0.5*p\n");
                printf(" Leff_nc = 2*%.1f+0.625*%.1f+0.5*%.1f\n",m*1e3,e*1e3,p*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        default:
            print_name("Error in Leff_table6_6_group");
            FATAL();
    }
    Leff1 = min(Leff_nc,Leff_cp);
    if(OUT)printf("  Leff1 = min(Leff_nc,Leff_cp) = min(%.1f,%.1f)\n",Leff_nc*1e3,Leff_cp*1e3);
    Leff2 = Leff_nc;
    if(OUT)printf("  Leff1 = %.1f mm\n",Leff1*1e3);
    if(OUT)printf("  Leff2 = %.1f mm\n",Leff2*1e3);
}

void EN1993_1_8_Table_6_6_Leff_individ(double &Leff1,double &Leff2, double Paramaters[7], type_LLU Row_table, bool OUT = true)
{
    bool DEBUG = false;
    double mx     = Paramaters[0];
    double w      = Paramaters[1];
    double e      = Paramaters[2];
    double ex     = Paramaters[3];
    double bp     = Paramaters[4];
    double _alpha = Paramaters[5];
    double m      = Paramaters[6];
    for(type_LLU i=0;i<7;i++)
    {
        if(Paramaters[i] < 0)
        {
            printf("Stupid check in EN1993_1_8_Table_6_6_Leff_individ\n");
            for(type_LLU j=0;j<7;j++)
                printf("Paramaters[%u] = %.1f mm\n",j,Paramaters[j]*1e3);
            FATAL();
        }
    }
    if(DEBUG)
    {
        for(type_LLU i=0;i<7;i++)
            printf("Paramaters[%u] = %.1f mm\n",i,Paramaters[i]*1e3);
    }
    double Leff_cp = 0.0;
    double Leff_nc = 0.0;
    switch(Row_table)
    {
        case(1):
        {
            double Leff_cp1 = 2*CONST_M_PI*mx;
            double Leff_cp2 = CONST_M_PI*mx+w;
            double Leff_cp3 = CONST_M_PI*mx+2*e;
            Leff_cp = min(Leff_cp1,min(Leff_cp2,Leff_cp3));
            double Leff_nc1 = 4*mx+1.25*ex;
            double Leff_nc2 = e+2*mx+0.625*ex;
            double Leff_nc3 = 0.5*bp;
            double Leff_nc4 = 0.5*w+2*mx+0.625*ex;
            Leff_nc = min(Leff_nc1,min(Leff_nc2,min(Leff_nc3,Leff_nc4)));
            if(OUT)
            {
                printf("Bolt-row location: Bolt-row outside tension flange of beam\n");
                printf(" Leff_cp1 = 2*PI*mx = 2*PI* %.1f = %.1f mm\n",mx*1e3,Leff_cp1*1e3);
                printf(" Leff_cp2 = PI*mx+w = PI* %.1f + %.1f = %.1f mm\n",mx*1e3,w*1e3,Leff_cp2*1e3);
                printf(" Leff_cp3 = PI*mx+2*e = PI* %.1f + 2* %.1f = %.1f mm\n",mx*1e3,e*1e3,Leff_cp3*1e3);
                printf(" Leff_cp = min(Leff_cp1,Leff_cp2,Leff_cp3) = min(%.1f,%.1f,%.1f) = %.1f mm",
                     Leff_cp1*1e3,Leff_cp2*1e3,Leff_cp3*1e3,Leff_cp*1e3);
                printf(" Leff_nc1 = 4*mx+1.25*ex = 4*%.1f + 1.25*%.1f = %.1f mm\n",mx*1e3,ex*1e3,Leff_nc1*1e3);
                printf(" Leff_nc2 = e+2*mx+0.625*ex = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n",e*1e3,mx*1e3,ex*1e3,Leff_nc2*1e3);
                printf(" Leff_nc3 = 0.5*bp = 0.5 * %.1f = %.1f mm\n",bp*1e3,Leff_nc3*1e3);
                printf(" Leff_nc4 = 0.5*w+2*mx+0.625*ex = 0.5* %.1f +2* %.1f +0.625* %.1f = %.1f mm\n",w*1e3,mx*1e3,ex*1e3,Leff_nc4*1e3);
                printf(" Leff_nc = min(Leff_nc1,Leff_nc2,Leff_nc3,Leff_nc4)\n");
                printf(" Leff_nc = min(%.1f,%.1f,%.1f,%.1f) = %.1f mm\n",
                        Leff_nc1*1e3,Leff_nc2*1e3,Leff_nc3*1e3,Leff_nc4*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(2):
        {
            Leff_cp = 2*CONST_M_PI*m;
            Leff_nc = _alpha * m;
            if(OUT)
            {
                printf("Bolt-row location: First bolt-row below tension flange of beam\n");
                printf(" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp*1e3);
                printf(" Leff_nc = alpha * m = %.2f * %.1f = %.1f mm\n",_alpha,m*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(3):
        {
            Leff_cp = 2*CONST_M_PI*m;
            Leff_nc = 4*m+1.25*e;
            if(OUT)
            {
                printf("Bolt-row location: Other inner bolt-row\n");
                printf(" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp*1e3);
                printf(" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n",m*1e3,e*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(4):
        {
            Leff_cp = 2*CONST_M_PI*m;
            Leff_nc = 4*m+1.25*e;
            if(OUT)
            {
                printf("Bolt-row location: Other end bolt-row\n");
                printf(" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp*1e3);
                printf(" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n",m*1e3,e*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(5):
        {
            double Leff_cp1 = 2*CONST_M_PI*mx;
            double Leff_cp2 = 2*m+0.625*e+ex;
            double Leff_cp3 = 2*CONST_M_PI*m;
            Leff_cp = min(Leff_cp1,min(Leff_cp2,Leff_cp3));
            double Leff_nc1 = _alpha*mx-(2*mx+0.625*e)+ex;
            double Leff_nc2 = 4*m+1.25*e;
            double Leff_nc3 = e+2*mx+0.625*ex;
            double Leff_nc4 = 4*mx+1.25*ex;
            Leff_nc = min(Leff_nc1,min(Leff_nc2,min(Leff_nc3,Leff_nc4)));
            if(OUT)
            {
                printf("Bolt-row location: bolt-row with gusset between bolts\n");
                printf(" Leff_cp1 = 2*PI*mx = 2*PI* %.1f = %.1f mm\n",mx*1e3,Leff_cp1*1e3);
                printf(" Leff_cp2 = 2*m+0.625*e+ex = 2*%.1f+0.625*%.1f+%.1f = %.1f mm\n",m*1e3,e*1e3,ex*1e3,Leff_cp2*1e3);
                printf(" Leff_cp3 = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp3*1e3);
                printf(" Leff_nc1 = alpha*mx-(2*mx+0.625*e)+ex = %.1f*%.1f-(2*%.1f+0.625*%.1f)+%.1f = %.1f mm\n",_alpha*1e3,mx*1e3,mx*1e3,e*1e3,ex*1e3,Leff_nc1*1e3);
                printf(" Leff_nc2 = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n",m*1e3,e*1e3,Leff_nc2*1e3);
                printf(" Leff_nc3 = e+2*mx+0.625*ex = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n",e*1e3,mx*1e3,ex*1e3,Leff_nc3*1e3);
                printf(" Leff_nc4 = 4*mx+1.25*ex = 4*%.1f + 1.25*%.1f = %.1f mm\n",mx*1e3,ex*1e3,Leff_nc4*1e3);
            }
            break;
        }
        default:
            print_name("Error in Leff_table6_6_individ");
            FATAL();
    }
    Leff1 = min(Leff_nc,Leff_cp);
    if(OUT)printf("Leff1 = min(Leff_nc,Leff_cp) = min(%.1f,%.1f)\n",Leff_nc*1e3,Leff_cp*1e3);
    Leff2 = Leff_nc;
    if(OUT)printf("Leff1 = %.1f mm\n",Leff1*1e3);
    if(OUT)printf("Leff2 = %.1f mm\n",Leff2*1e3);
}

bool EN1993_1_8_Table_3_3_e1(double Do, double e1, double THK, bool OUT = true)
{
    bool ex = false;
    if(1.2*Do <= e1 && e1 <= 4.*THK+0.040)
        ex = true;
    if(OUT)
    {
        printf("Check end distanse e1 (on load direction) by table 3.3 EN1993-1-8:\n");
        printf("1.2 * Do <= e1\n");
        printf("1.2 * %.1f <= %.1f\n",Do*1e3,e1*1e3);
        printf("%.1f <= %.1f\n",1.2 * Do*1e3,e1*1e3);
        printf("e1 <= 4*THK+0.040\n");
        printf("%.1f <= 4*%.1f + 0.040\n",e1*1e3,THK*1e3);
        printf("%.1f <= %.1f ",e1*1e3,(4.*THK+0.040)*1e3);
        if(ex)
              printf("- OK\n");
        else  {printf("- None\n");FATAL();}
    }
    return ex;
}

bool EN1993_1_8_Table_3_3_e2(double Do, double e2, double THK, bool OUT = true)
{
    bool ex = false;
    if(1.2*Do <= e2 )//&& e2 <= 4*THK+0.040)
        ex = true;
    if(OUT)
    {
        printf("Check end distanse e2 (perpendicular on load direction) by table 3.3 EN1993-1-8:\n");
        printf("1.2 * Do <= e2\n");
        printf("1.2 * %.1f <= %.1f\n",Do*1e3,e2*1e3);
        printf("%.1f <= %.1f \n",1.2 * Do*1e3,e2*1e3);
        printf("e2 <= 4*THK+0.040\n");
        printf("%.1f <= 4*%.1f+0.040\n",e2*1e3,THK*1e3);
        printf("%.1f <= %.1f\n",e2*1e3,(4*THK+0.040)*1e3);
        if(ex)
              printf("- OK\n");
        else
        {
            //printf("- None\n");
            printf("1.2*Do <= e2 - ");
            if(1.2*Do <= e2) printf("OK\n"); else printf("None\n");
            printf("e2 <= 4*THK+0.040 - ");
            if(e2 <= 4*THK+0.040) printf("OK\n"); else printf("None\n");
            FATAL();
        }
    }
    return ex;
}

bool EN1993_1_8_Table_3_3_p1(double Do, double p1, double THK, bool OUT = true)
{
    bool ex = false;
    if(2.2*Do <= p1 && p1 <= min(14*THK,0.200))
        ex = true;
    if(OUT)
    {
        printf("Check spacing p1 by table 3.3 EN1993-1-8:\n");
        printf("2.2 * Do <= p1\n");
        printf("2.2 * %.1f <= %.1f\n",Do*1e3,p1*1e3);
        printf("%.1f <= %.1f ",2.2 * Do*1e3,p1*1e3);
        if(ex)
              printf("- OK\n");
        else  {printf("- None\n");FATAL();}
    }
    return ex;
}

bool EN1993_1_8_punkt_4_5_2(double a, bool OUT = true)
{
    bool ex = false;
    if(a >= 0.003)
        ex = true;
    if(OUT)
    {
        printf("Check minimum weld thickness by p.4.5.2(2) EN1993-1-8: ");
        if(ex) printf("OK\n");
        else
        {
            printf("None\n");
            FATAL();
        }
    }
    return ex;
}


bool EN1993_1_8_punkt_4_5_1(double a, double Leff, bool OUT = true)
{
    bool ex = false;
    if(Leff >= min(6*a,0.030))
        ex = true;
    if(OUT)
    {
        printf("Check lenght of weld  by p.4.5.1 EN1993-1-8: ");
        if(ex) printf("OK\n");
        else
        {
            printf("None\n");
            FATAL();
        }
    }
    return ex;
}


void EN1993_1_8_Formula3_11(double e2, double Do, double t, double Fu, double gamma_M2, double &Nu_Rd, bool OUT = true)
{
    Nu_Rd = 2.0*(e2 - 0.5*Do)*t*Fu/gamma_M2;
    if(OUT)
    {
        printf("Calculate by Formula 3.11 EN1993-1-8\n");
        printf("Tension force by a single row of one bolt in one leg\n");
        printf("Nu_Rd = 2.0*(e2-0.5*Do)*t*Fu/gamma_M2\n");
        printf("Nu_Rd = 2.0*(%.1f-0.5*%.1f)*%.1f*%.1f/%.2f\n",e2*1e3,Do*1e3,t*1e3,Fu*1e-6,gamma_M2);
        printf("Nu_Rd = %.1f kN\n",Nu_Rd*1e-3);
    }
}


void EN1993_1_8_Formula3_12(double p1, double Do, double Anet, double Fu, double gamma_M2, double &Nu_Rd, bool OUT = true)
{
    Array<double> x,y;
    x.AddList(2,2.5*Do,5.0*Do);
    y.AddList(2,0.4000,0.7000);
    double betta2;
    Graph(x,y,p1,betta2);
    Nu_Rd = betta2*Anet*Fu/gamma_M2;
    if(OUT)
    {
        printf("Calculate by Formula 3.12 EN1993-1-8\n");
        printf("Tension force by a single row of two bolts in one leg\n");
        printf("Nu_Rd = betta2*Anet*Fu/gamma_M2\n");
        printf("Nu_Rd = %.3f*%.1f*%.1f/%.2f\n",betta2,Anet*1e6,Fu*1e-6,gamma_M2);
        printf("Nu_Rd = %.1f kN\n",Nu_Rd*1e-3);
    }
    x.Delete();
    y.Delete();
}


void EN1993_1_8_Formula3_13(double p1, double Do, double Anet, double Fu, double gamma_M2, double &Nu_Rd, bool OUT = true)
{
    Array<double> x,y;
    x.AddList(2,2.5*Do,5.0*Do);
    y.AddList(2,0.5000,0.7000);
    double betta3;
    Graph(x,y,p1,betta3);
    Nu_Rd = betta3*Anet*Fu/gamma_M2;
    if(OUT)
    {
        printf("Calculate by Formula 3.13 EN1993-1-8\n");
        printf("Tension force by a single row of 3 or more bolts in one leg\n");
        printf("Nu_Rd = betta2*Anet*Fu/gamma_M2\n");
        printf("Nu_Rd = %.3f*%.1f*%.1f/%.2f\n",betta3,Anet*1e6,Fu*1e-6,gamma_M2);
        printf("Nu_Rd = %.1f kN\n",Nu_Rd*1e-3);
    }
    x.Delete();
    y.Delete();
}



///
/// Table 6.5
///

void EN1993_1_8_Table_6_5_Leff_group(double &Leff1,double &Leff2, double Paramaters[5], type_LLU Row_table, bool OUT = true)
{
    double m      = Paramaters[0];
    double p      = Paramaters[1];
    double _alpha = Paramaters[2];
    double e      = Paramaters[3];
    double e1     = Paramaters[4];
    for(type_LLU i=0;i<5;i++)
    {
        if(Paramaters[i] < 0)
        {
            printf("Stupid check in EN1993_1_8_Table_6_5_Leff_group\n");
            for(type_LLU j=0;j<5;j++)
                printf("Paramaters[%u] = %.1f mm",j,Paramaters[j]*1e3);
            FATAL();
        }
    }
    double Leff_cp = 0.0;
    double Leff_nc = 0.0;
    switch(Row_table)
    {
        case(1):
        {
            Leff_cp = CONST_M_PI*m+p;
            Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e);
            if(OUT)
            {
                printf("Bolt-row location: Bolt-row adjacent to a stiffener\n");
                printf(" Leff_cp = PI*m+p\n");
                printf(" Leff_cp = PI*%.1f+%.1f\n",m*1e3,p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e)\n");
                printf(" Leff_nc = 0.5*%.1f+%.2f*%.1f-(2*%.1f+0.625*%.1f)\n",p*1e3,_alpha,m*1e3,m*1e3,e*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        case(2):
        {
            Leff_cp = 2*p;
            Leff_nc = p;
            if(OUT)
            {
                printf("Bolt-row location: Other inner bolt-row\n");
                printf(" Leff_cp = 2*p\n");
                printf(" Leff_cp = 2*%.1f\n",p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = p\n");
                printf(" Leff_nc = %.1f\n",p*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        case(3):
        {
            Leff_cp = min(CONST_M_PI*m+p,2*e1+p);
            Leff_nc = min(2*m+0.625*e+0.5*p,e1+0.5*p);
            if(OUT)
            {
                printf("Bolt-row location: Other end bolt-row\n");
                printf(" Leff_cp = min(PI*m+p,2*e1+p)\n");
                printf(" Leff_cp = min(PI*%.1f+%.1f,2*%.1f+%.1f)\n",m*1e3,p*1e3,e1*1e3,p*1e3);
                printf(" Leff_cp = %.1f mm\n",Leff_cp*1e3);
                printf(" Leff_nc = min(2*m+0.625*e+0.5*p,e1+0.5*p)\n");
                printf(" Leff_nc = min(2*%.1f+0.625*%.1f+0.5*%.1f,%.1f+0.5*%.1f)\n",m*1e3,e*1e3,p*1e3,e1*1e3,p*1e3);
                printf(" Leff_nc = %.1f mm\n",Leff_nc*1e3);
            }
            break;
        }
        case(4):
        {
            print_name("NOT RELEVANT");
            FATAL();
            break;
        }
        default:
            print_name("Error in Leff_table6_5_group");
            FATAL();
    }
    Leff1 = min(Leff_nc,Leff_cp);
    if(OUT)printf("Leff1 = min(Leff_nc,Leff_cp) = min(%.1f,%.1f)\n",Leff_nc*1e3,Leff_cp*1e3);
    Leff2 = Leff_nc;
    if(OUT)printf("Leff1 = %.1f mm\n",Leff1*1e3);
    if(OUT)printf("Leff2 = %.1f mm\n",Leff2*1e3);
}

void EN1993_1_8_Table_6_5_Leff_individ(double &Leff1,double &Leff2, double Paramaters[4], type_LLU Row_table, bool OUT = true)
{
    bool DEBUG = false;
    double e      = Paramaters[0];
    double e1     = Paramaters[1];
    double _alpha = Paramaters[2];
    double m      = Paramaters[3];
    for(type_LLU i=0;i<6;i++)
    {
        if(Paramaters[i] < 0)
        {
            printf("Stupid check in EN1993_1_8_Table_6_5_Leff_individ\n");
            for(type_LLU j=0;j<6;j++)
                printf("Paramaters[%u] = %.1f mm\n",j,Paramaters[j]*1e3);
            FATAL();
        }
    }
    if(DEBUG)
    {
        for(type_LLU i=0;i<6;i++)
            printf("Paramaters[%u] = %.1f mm\n",i,Paramaters[i]*1e3);
    }
    double Leff_cp = 0.0;
    double Leff_nc = 0.0;
    switch(Row_table)
    {
        case(1):
        {
            Leff_cp = 2*CONST_M_PI*m;
            Leff_nc = _alpha * m;
            if(OUT)
            {
                printf("Bolt-row location: Bolt-row adjacent to a stiffener\n");
                printf(" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp*1e3);
                printf(" Leff_nc = alpha * m = %.2f * %.1f = %.1f mm\n",_alpha,m*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(2):
        {
            Leff_cp = 2*CONST_M_PI*m;
            Leff_nc = 4*m+1.25*e;
            if(OUT)
            {
                printf("Bolt-row location: Other inner bolt-row\n");
                printf(" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n",m*1e3,Leff_cp*1e3);
                printf(" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n",m*1e3,e*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(3):
        {
            double Leff_cp1 = 2*CONST_M_PI*m;
            double Leff_cp2 = CONST_M_PI*m+2*e1;
            Leff_cp = min(Leff_cp1,Leff_cp2);
            double Leff_nc1 = 4*m+1.25*e;
            double Leff_nc2 = e1+2*m+0.625*e;
            Leff_nc = min(Leff_nc1,Leff_nc2);
            if(OUT)
            {
                printf("Bolt-row location: Other end bolt-row\n");
                printf(" Leff_cp1 = 2*PI*m = 2*PI* %.1f = %.1f mm\n",m*1e3,Leff_cp1*1e3);
                printf(" Leff_cp2 = PI*m+2*e1 = PI* %.1f + 2* %.1f = %.1f mm\n",m*1e3,e1*1e3,Leff_cp2*1e3);
                printf(" Leff_cp = min(Leff_cp1,Leff_cp2) = min(%.1f,%.1f) = %.1f mm",
                     Leff_cp1*1e3,Leff_cp2*1e3,Leff_cp*1e3);
                printf(" Leff_nc1 = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n",m*1e3,e*1e3,Leff_nc1*1e3);
                printf(" Leff_nc2 = e1+2*m+0.625*e = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n",e1*1e3,m*1e3,e*1e3,Leff_nc2*1e3);
                printf(" Leff_nc = min(%.1f,%.1f) = %.1f mm\n",
                        Leff_nc1*1e3,Leff_nc2*1e3,Leff_nc*1e3);
            }
            break;
        }
        case(4):
        {
            double Leff_cp1 = 2*CONST_M_PI*m;
            double Leff_cp2 = CONST_M_PI*m+2*e1;
            Leff_cp = min(Leff_cp1,Leff_cp2);
            double Leff_nc  = e1+_alpha*m-(2*m+0.625*e);
            if(OUT)
            {
                printf("Bolt-row location: End bolt-row adjacent to a stiffener\n");
                printf(" Leff_cp1 = 2*PI*m = 2*PI* %.1f = %.1f mm\n",m*1e3,Leff_cp1*1e3);
                printf(" Leff_cp2 = PI*m+2*e1 = PI* %.1f + 2* %.1f = %.1f mm\n",m*1e3,e1*1e3,Leff_cp2*1e3);
                printf(" Leff_cp = min(Leff_cp1,Leff_cp2) = min(%.1f,%.1f) = %.1f mm",
                     Leff_cp1*1e3,Leff_cp2*1e3,Leff_cp*1e3);
                printf(" Leff_nc = e1+alpha*m-2*m+0.625*e = %.1f + %.3f*%.1f- 2* %.1f + 0,625* %.1f = %.1f mm\n",e1*1e3,_alpha,m*1e3,m*1e3,e*1e3,Leff_nc*1e3);
            }
            break;
        }
        default:
            print_name("Error in Leff_table6_6_individ");
            FATAL();
    }
    Leff1 = min(Leff_nc,Leff_cp);
    if(OUT)printf("Leff1 = min(Leff_nc,Leff_cp) = min(%.1f,%.1f)\n",Leff_nc*1e3,Leff_cp*1e3);
    Leff2 = Leff_nc;
    if(OUT)printf("Leff1 = %.1f mm\n",Leff1*1e3);
    if(OUT)printf("Leff2 = %.1f mm\n",Leff2*1e3);
}

