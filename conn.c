//
//        ***************--------
//        *             *       |
//        *   =======   *       |
//        *      |      *       |
//        *      |      *       |
//        *      |      *       Yplus
//        *      |      *       |
//        *      |      *       |
//        *   ======= -----------
//        *             *       Yminus
//        ***************--------
//        |             |
//        |-----X-------|
//

class EndPlate {
    double Py, Pu;

public:
    double X, Yplus, Yminus, THK;
    double Fy() { return Py; }
    double Fu() { return Pu; }
    char* Category;
    EndPlate(double thk, double x, double yplus, double yminus, const char* _Category)
    {
        X = x;
        Yplus = yplus;
        Yminus = yminus;
        THK = thk;
        Category = new char[10];
        strcpy(Category, _Category);
        SNiP(Category, THK, Py, Pu);
    };
    EndPlate(double thk, double x, double y, IS* I, char* _Category)
    {
        X = x;
        Yplus = I->h() + (y - I->h()) / 2.;
        Yminus = -(y - I->h()) / 2.;
        THK = thk;
        Category = _Category;
        SNiP(Category, THK, Py, Pu);
    };
    void Printf()
    {
        printf("End plate\n");
        if (true)
            printf("Size %.1f x (%.1f;%.1f) mm\n", X * 1e3, Yplus * 1e3, Yminus * 1e3);
        else
            printf("Size %.1f x %.1f mm\n", X * 1e3, (Yplus + Yminus) * 1e3);
        printf("Thickness = %.1fmm\n", THK * 1e3);
        printf("Py = %.1fMPa\n", Py * 1e-6);
        printf("Pu = %.1fMPa\n\n", Pu * 1e-6);
    };
};


// TYPE 1
//        ***************--------
//        *             *       |
//        *             *       A
//        *             *       |
//        ***************--------
//        |             |
//        |----- B -----|
//
// TYPE 2
//        ******-------------
//        *     \           |
//        *      \          A
//        *      *          |
//        ********-----------
//        |      |
//        |-- B -|
//

enum Gusset_Type { Gusset_Type1,
    Gusset_Type2 };
class Gusset {
    bool LIFE;
    double Py, Pu;
    Gusset_Type type;
    void IN(double thk, double a, double b, double _a_A, double _a_B, Gusset_Type t, const char* _Category)
    {
        LIFE = true;
        A = a;
        B = b;
        THK = thk;
        a_A = _a_A;
        a_B = _a_B;
        if (a_A < 1e-6)
            a_A = max(0.003, thk * 0.85 - 0.00025);
        if (a_B < 1e-6)
            a_B = max(0.003, thk * 0.85 - 0.00025);
        t = type;
        Category = new char[10];
        //        strncpy (Category,_Category,4);
        //        Category[4] = '\0';
        strcpy(Category, _Category);
        SNiP(Category, THK, Py, Pu);
    };

public:
    double A, B, THK;
    double a_A, a_B; //Weld
    double Fy() { return Py; }
    double Fu() { return Pu; }
    char* Category;
    bool Life() { return LIFE; }
    Gusset() { LIFE = false; }
    Gusset(double thk, double a, double b, double _a_A, double _a_B, Gusset_Type t, const char* _Category)
    {
        IN(thk, a, b, _a_A, _a_B, t, _Category);
    };
    Gusset(double thk, double a, double b, Gusset_Type t, const char* _Category)
    {
        IN(thk, a, b, 0, 0, t, _Category);
    }
    void Printf(bool OUT = true)
    {
        printf("Gusset\n");
        printf("Size %.1f x %.1f mm (%s)\n", A * 1e3, B * 1e3, Category);
        printf("Thickness = %.1fmm\n", THK * 1e3);
        printf("Py = %.1fMPa\n", Py * 1e-6);
        printf("Pu = %.1fMPa\n", Pu * 1e-6);
        printf("Weld of gusset:\n");
        printf("a_A = %.1fmm\n", a_A * 1e3);
        printf("a_B = %.1fmm\n", a_B * 1e3);
    };
    bool Check(bool OUT = true)
    {
        bool ex = true;
        if (OUT)
            printf("Checking of gusset:\n");
        if (THK < 1e-6)
            FATAL();
        if (!EN1993_1_8_punkt_4_5_2(a_A, OUT))
            ex = false;
        if (!SNiP_II_23_punkt12_8(a_A, THK, OUT))
            ex = false;
        if (!EN1993_1_8_punkt_4_5_1(a_A, A, OUT))
            ex = false;

        if (!EN1993_1_8_punkt_4_5_2(a_B, OUT))
            ex = false;
        if (!SNiP_II_23_punkt12_8(a_B, THK, OUT))
            ex = false;
        if (!EN1993_1_8_punkt_4_5_1(a_B, B, OUT))
            ex = false;
        if (!ex)
            FATAL();
        return ex;
    }
    double Calculate(double Fed, double Ved, bool OUT = true)
    {
        double factor = 0;
        Printf(OUT);
        Check(OUT);
        Plate AA(THK, A);
        Plate BB(THK, B);
        // check 1
        {
            double Npl_Rd;
            EN1993_1_1_Formula6_6(AA.A(), Py, 1.1, Npl_Rd, OUT);
            Printf_CALC(Fed / Npl_Rd, factor, OUT);
        }
        // check 2
        {
            double Npl_Rd;
            EN1993_1_1_Formula6_6(BB.A(), Py, 1.1, Npl_Rd, OUT);
            Printf_CALC(Fed / Npl_Rd, factor, OUT);
        }
        // check 3
        {
            double Vpl_Rd;
            EN1993_1_1_Formula6_18(Py, AA.A(), 1.1, Vpl_Rd, OUT);
            Printf_CALC(Ved / Vpl_Rd, factor, OUT);
        }
        // check 4
        {
            double Vpl_Rd;
            EN1993_1_1_Formula6_18(Py, BB.A(), 1.1, Vpl_Rd, OUT);
            Printf_CALC(Ved / Vpl_Rd, factor, OUT);
        }
        // Weld
        Plate a_AA(a_A, A);
        Plate a_BB(a_B, B);
        // check 5
        {
            if (OUT)
                printf("Calculation: Weld checking\n");
            double AreaW = a_AA.A() * 2;
            double Ww = a_AA.Wx() * 2;
            double Nw = fabs(Ved);
            double Mw = fabs(Fed) * a_BB.h() / 2.;
            double Vw = fabs(Fed);
            double Sigma_perp = Nw / AreaW / sqrt(2.) + Mw / Ww / sqrt(2.);
            double Tau_perp = Nw / AreaW / sqrt(2.) + Mw / Ww / sqrt(2.);
            double Tau_parral = Vw / AreaW;

            Printf_CALC(EN1993_1_8_Formula4_1(Sigma_perp, Tau_perp, Tau_parral,
                            Pu, Category, 1.25, OUT),
                factor, OUT);
        }
        // check 6
        {
            if (OUT)
                printf("Calculation: Weld checking\n");
            double AreaW = a_BB.A() * 2;
            double Ww = a_BB.Wx() * 2;
            double Nw = fabs(Fed);
            double Mw = 0;
            double Vw = fabs(Ved);
            double Sigma_perp = Nw / AreaW / sqrt(2.) + Mw / Ww / sqrt(2.);
            double Tau_perp = Nw / AreaW / sqrt(2.) + Mw / Ww / sqrt(2.);
            double Tau_parral = Vw / AreaW;

            Printf_CALC(EN1993_1_8_Formula4_1(Sigma_perp, Tau_perp, Tau_parral,
                            Pu, Category, 1.25, OUT),
                factor, OUT);
        }
        return factor;
    }
};

///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
double Calc_MomentConnection
     (
        IS*       I,
        EndPlate* E,
        Bolt* Bolts,
        double   af,
        double   aw,
        double gamma_M0,
        double gamma_M1,
        Array<STAAD_FORCE_FULL> arrSTAAD_FORCES,
        Gusset *UP_Hanged1,
        Gusset *UP_Hanged2,
        Gusset *DOWN_Hanged1,
        Gusset *DOWN_Hanged2,
        bool OUT = true
     )
{
    double factor = 0;
    Array<double> f;f.SetSize(arrSTAAD_FORCES.GetSize());
    for(type_LLU i=0;i<arrSTAAD_FORCES.GetSize();i++)
    {
        f.Set(i,0.0);
        arrSTAAD_FORCES.Get(i).Printf(NULL,true);
        f[i] = Calc_MomentConnection(I, E, Bolts, af, aw, gamma_M0, gamma_M1, arrSTAAD_FORCES.Get(i),
                                          UP_Hanged1, UP_Hanged2, DOWN_Hanged1, DOWN_Hanged2, OUT);
        Printf_CALC(f[i],factor,OUT);
    }
    if(OUT)
    {
        printf("Result of calculation for %u cases : %.1f %%\n",arrSTAAD_FORCES.GetSize(),factor*1e2);
        printf("Table. Factor value\n");
        type_LLU j=0;
        for(type_LLU i=0;i<f.GetSize();i++)
        {
            if(f[i]>0.75)
            {
                j++;
                printf("%u\t%2.1f%%\t",i,f[i]*1e2);
                arrSTAAD_FORCES.Get(i).Printf(NULL,true);
            }
        }
        printf("Number of cases with factor more 75%% is %u\n",j);
    }
    return factor;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
//double Calc_MomentConnection
//     (
//        IS*       I,
//        STAAD_FORCE_FULL Force,
//        bool THK_END_PLATE_ANALYSYS = false,
//        bool OUT = true
//     )
//{
//    bool DEBUG = false;//true;
//    if(DEBUG) printf("DEBUG OUTPUT\n");
//    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
//    {
//
//Calc_MomentConnection1:
//
//        if(DEBUG) printf("DEBUG: Diameter of bolts: %.0f\n",HB_Diameter_Bolt[i]*1000);
//        BOLT_CLASS bc = g5_6;
//        if(HB_Diameter_Bolt.Get(i)>0.016)
//            bc = g8_8;
//        if(HB_Diameter_Bolt.Get(i)>0.031)
//        {
//            if(DEBUG) print_name("NOT CONNECTION");
//            break;//return 10;
//        }
//
//
//        double A,M,L;ANUREV(HB_Diameter_Bolt.Get(i), A, M, L, false);
//        double X_bolts = max(I->b()/2.,I->tw()+2*M);
//        X_bolts = (((type_LLU)(X_bolts*1000./10.))+1)*10./1000.;
//        double Y_bolts = min(I->h()-2.*I->tf()-2*M,I->h()-2*1.2*(HB_Diameter_Bolt.Get(i)+0.003));
//        Y_bolts = (((type_LLU)(Y_bolts*1000./10.))-1)*10./1000.;
//        double X_plate = max(X_bolts+4.*(HB_Diameter_Bolt.Get(i)+0.003),I->b());
//        X_plate = (((type_LLU)(X_plate*1000./10.))+1)*10./1000.;
//        double Y_plate = I->h();
//
//        if((I->h()-Y_bolts)/2. < max(1.2*(HB_Diameter_Bolt.Get(i)+0.0029),I->tf()+M) ||
//           Y_bolts<2.200*(HB_Diameter_Bolt.Get(i)+0.003) ||
//           Y_bolts<A)
//        {
//            if(DEBUG)
//            {
//                print_name("NOT CONNECTION e2");
//                I->Printf();
//                printf("Y_bolts = %.1f mm\n",Y_bolts*1e3);
//                printf("M = %.1f mm\n",M*1e3);
//                printf("HB_Diameter_Bolt.Get(i) = %.1f mm\n",HB_Diameter_Bolt.Get(i)*1e3);
//            }
//            break;
//        }
//        if(Y_bolts+2*M+2*I->tf() >= I->h())
//        {
//            if(DEBUG) print_name("Don`t find good connection by Calc_MomentConnection");
//            i++; goto Calc_MomentConnection1;
//        }
//
//        Gusset gNULL;
//        Bolt      B(HB_Diameter_Bolt.Get(i),4,bc,STANDART_Eurocode,1.25);
//                  B.Rectangle(X_bolts,Y_bolts,2,2);
//                  B.Move     (0.00000,I->h()/2.00);
//
//        double thk_plate = THK_more(max(I->tf(),I->tw()));
//        if(!THK_END_PLATE_ANALYSYS)thk_plate = THK_more(HB_Diameter_Bolt.Get(i));
//        for(;thk_plate<THK_more(HB_Diameter_Bolt.Get(i))+0.0001;thk_plate+=0.002)
//        {
//            thk_plate = THK_more(thk_plate);
//            EndPlate  E(thk_plate,X_plate,Y_plate,0.0,"C235");
//            double factor = Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&gNULL,&gNULL,DEBUG);
//            if(factor < 1.00)
//            {
//                // CORRECT OUTPUT //
//                //if(DEBUG)Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&gNULL,&gNULL,true);
//                if(OUT)
//                {
//                    printf("Short output about moment connection:\n");
//                    printf("Section  : %s\n",I->name);
//                    printf("End plate: %.1f x %.1f x %.1f mm\n",E.THK*1e3,E.X*1e3,(E.Yplus-E.Yminus)*1e3);
//                    printf("Diameter of bolts: HM%.0f\n",B.Diameter_()*1e3);
//                    printf("Number of bolts = %2u\n",B.positionBolts.GetSize());
//                    for(type_LLU bb = 0;bb<B.positionBolts.GetSize();bb++)
//                        printf("%2u - [%.1f,%.1f] mm\n",bb,
//                                                        B.positionBolts.Get(bb).x*1e3,
//                                                        B.positionBolts.Get(bb).y*1e3);
//                    Printf(bc);
//                }
//                return factor;
//            }
//        }
//    }
//    return 10;
//}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
double Calc_MomentConnection_WithDown
     (
        IS*       I,
        STAAD_FORCE_FULL Force,
        bool THK_END_PLATE_ANALYSYS = false,
        bool OUT = true
     )
{
    bool DEBUG =  false;//true;//
    if(DEBUG) printf("DEBUG OUTPUT\n");
    // with Down_hanges
    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {

Calc_MomentConnection_WithDown1:

        if(DEBUG) printf("DEBUG: Diameter of bolts: %.0f\n",HB_Diameter_Bolt[i]*1000);
        BOLT_CLASS bc = g5_6;
        if(HB_Diameter_Bolt.Get(i)>0.016)
            bc = g8_8;
        if(HB_Diameter_Bolt.Get(i)>0.031)
        {
            if(DEBUG) print_name("NOT CONNECTION");
            break;//return 10;
        }

        double A,M,L;ANUREV(HB_Diameter_Bolt.Get(i), A, M, L, false);
        double X_bolts = max(I->b()/2.,I->tw()+2*M);
        X_bolts = (((type_LLU)(X_bolts*1000./10.))+1)*10./1000.;
        double Y_bolts = min(I->h()-2.*I->tf()-2*M,I->h()-2*1.2*(HB_Diameter_Bolt.Get(i)+0.003));
        Y_bolts = (((type_LLU)(Y_bolts*1000./10.))-1)*10./1000.;
        double X_plate = max(X_bolts+4.*(HB_Diameter_Bolt.Get(i)+0.003),I->b());
        X_plate = (((type_LLU)(X_plate*1000./10.))+1)*10./1000.;

        if((I->h()-Y_bolts)/2. < max(1.2*(HB_Diameter_Bolt.Get(i)+0.0029),I->tf()+M) ||
           Y_bolts<2.200*(HB_Diameter_Bolt.Get(i)+0.003) ||
           Y_bolts<A)
        {
            if(DEBUG)
            {
                print_name("NOT CONNECTION e2");
                I->Printf();
                printf("Y_bolts = %.1f mm\n",Y_bolts*1e3);
                printf("M = %.1f mm\n",M*1e3);
                printf("HB_Diameter_Bolt.Get(i) = %.1f mm\n",HB_Diameter_Bolt.Get(i)*1e3);
            }
            break;
        }
        if(Y_bolts+2*M+2*I->tf() >= I->h())
        {
            if(DEBUG) print_name("Don`t find good connection by Calc_MomentConnection_WithDown");
            i++; goto Calc_MomentConnection_WithDown1;//break;//return 10;
        }

        // vertical distance between bolts in Down_hanges to profile
        double L2 = 2.*tan(RADIANS(30))*(L/2./tan(RADIANS(30))-HB_Diameter_Bolt.Get(i)+(X_plate-X_bolts)/2.);
        L2 = max(L2,max(L2/2.,2.5*HB_Diameter_Bolt.Get(i)-(I->h()-Y_bolts)/2.)+max(2.*HB_Diameter_Bolt.Get(i),L2/2.));
        L2 = (((type_LLU)(L2*1000./5.))+1)*5./1000.;

        double Y_plate = I->h() + L2 + I->tf();
        Y_plate = (((type_LLU)(Y_plate*1000./10.))+1)*10./1000.;



        Gusset gNULL;
        Gusset g1(THK_more(I->tw()),Y_plate-THK_more(I->tf())-I->h(),2*(Y_plate-THK_more(I->tf())),Gusset_Type2,"C235");
        Gusset g2(THK_more(I->tf()),X_plate                  ,2*(Y_plate-THK_more(I->tf())),Gusset_Type1,"C235");

        double thk_plate = THK_more(max(I->tf(),I->tw()));
        if(!THK_END_PLATE_ANALYSYS)thk_plate = THK_more(HB_Diameter_Bolt.Get(i));
        for(;thk_plate<THK_more(HB_Diameter_Bolt.Get(i))+0.0001;thk_plate+=0.002)
        {
            double y_p;
            Bolt      B(HB_Diameter_Bolt.Get(i),6,bc,STANDART_Eurocode,1.25);
                      y_p = min(4*thk_plate+0.040-0.001,L2/2.);
                      B.SetXY(0,-X_bolts/2.,y_p);
                      B.SetXY(1,+X_bolts/2.,y_p);
                      y_p = (I->h()-Y_bolts)/2.+(Y_plate-I->h());
                      B.SetXY(2,-X_bolts/2.,y_p);
                      B.SetXY(3,+X_bolts/2.,y_p);
                      y_p = Y_plate-I->h()/2.+Y_bolts/2.;//(I->h()-Y_bolts)/2.+y_p;
                      B.SetXY(4,-X_bolts/2.,y_p);
                      B.SetXY(5,+X_bolts/2.,y_p);

            thk_plate = THK_more(thk_plate);
            EndPlate  E(thk_plate,X_plate,I->h(),I->h()-Y_plate,"C235");
            double factor = Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&g1,&g2,DEBUG);
            if(factor < 1.00)
            {
                // CORRECT OUTPUT //
                //if(DEBUG)Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&gNULL,&gNULL,true);
                if(OUT)
                {
                    printf("Short output about moment connection with down hanges:\n");
                    printf("Section  : %s\n",I->name);
                    printf("End plate: %.1f x %.1f x %.1f mm\n",E.THK*1e3,E.X*1e3,(E.Yplus-E.Yminus)*1e3);
                    printf("Diameter of bolts: HM%.0f\n",B.Diameter_()*1e3);
                    printf("Number of bolts = %2u\n",B.positionBolts.GetSize());
                    for(type_LLU bb = 0;bb<B.positionBolts.GetSize();bb++)
                        printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                        B.positionBolts.Get(bb).x*1e3,
                                                        B.positionBolts.Get(bb).y*1e3);
                    Printf(bc);
                    g1.Printf(true);
                    g2.Printf(true);
                }
                return factor;
            }
        }
    }

    return 10;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
double Calc_FinPlate_Angle
    (
        Angle a,
        double angle,
        double ForceKnee,
        bool OUT = true
    )
{

    bool DEBUG = false;//true;
    if(DEBUG)printf("\n\nForceKnee = %.1f N\n\n",ForceKnee);
    double Py,Pu;
    SNiP("C235",a.thk(),Py,Pu);

    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {
        double Do = HB_Diameter_Bolt.Get(i)+0.003;
        BOLT_CLASS bc = g5_6;
        if(HB_Diameter_Bolt.Get(i)>0.016)
            bc = g8_8;
        if(HB_Diameter_Bolt.Get(i)>0.025)
        {
            print_name("NOT CONNECTION");
            return 10;
        }
        double e1 = (double)((type_LLU)(2*Do*1000/5.)+1)*5/1000.;
        double p1 = (double)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
        double e3 = e1+0.025/cos(angle)+a.b()/2.*tan(angle);
               e3 = (double)((type_LLU)(e3*1000/5.)+1)*5/1000.;
        for(type_LLU n=2;n<4;n++) // number of bolts
        {
            Bolt Bo(HB_Diameter_Bolt.Get(i),1*n,bc,STANDART_Eurocode,1.25);
            if(n != 1)
            {
                Bo.Rectangle(0.000,(n-1)*p1,1,n);
                Bo.Move     (0.000,e3+(n-1)*p1/2.);
            }
            else
            {
               Bo.SetXY(0,0,e3);
            }
            if(DEBUG) Bo.Printf();
            double ff = 0;
            // check plate //
            double WitdhPlate = a.b()+2*0.005;
            Printf_CALC(
            Calc_FinPlate ( WitdhPlate, e1+(n-1)*p1+e3, WitdhPlate/2., 0,
                            THK_more(a.thk()), Pu, Py,
                            &Bo,
                            0, 1.0, 1.1,
                            0 , -ForceKnee, 0,
                            0 ,  0 , 0 ,
                            true , DEBUG),ff, DEBUG);
            double Nu_Rd;
            if(n == 1)
            {
                EN1993_1_8_Formula3_11(a.b()/2., Do, a.thk(), Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            else if(n == 2)
            {
                double Anet = a.A() - Do*a.thk();
                EN1993_1_8_Formula3_12(p1, Do, Anet, Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            else
            {
                double Anet = a.A() - Do*a.thk();
                EN1993_1_8_Formula3_13( p1, Do, Anet, Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            Printf_CALC(fabs(ForceKnee)/Nu_Rd,ff,DEBUG);
            if(ff < 1.00)
            {
                printf("Short output about connection:\n");
                Bo.Printf();
                printf("Size of Plate(%.1f x %.1f)mm\n",WitdhPlate*1e3,(e1+(n-1)*p1+e3)*1e3);
                printf("e1 = %.1f mm\n",e1*1e3);
                printf("p1 = %.1f mm\n",p1*1e3);
                printf("Number of bolts = %2u\n",n);
                return ff;
            }
        }
    }
    return 10;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
double Calc_FinPlate_Angle
    (
        Angle_LD a,
        double angle,
        double ForceKnee,
        bool OUT = true
    )
{

    bool DEBUG = false;//true;
    if(DEBUG)printf("\n\nForceKnee = %.1f N\n\n",ForceKnee);
    double Py,Pu;
    SNiP("C235",a.thk(),Py,Pu);

    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {
        double Do = HB_Diameter_Bolt.Get(i)+0.003;
        BOLT_CLASS bc = g5_6;
        if(HB_Diameter_Bolt.Get(i)>0.016)
            bc = g8_8;
        if(HB_Diameter_Bolt.Get(i)>0.025)
        {
            print_name("NOT CONNECTION");
            return 10;
        }
        double e1 = (double)((type_LLU)(2*Do*1000/5.)+1)*5/1000.;
        double p1 = (double)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
        double e3 = e1+0.025/cos(angle)+(a.THK_PLATE()+2*a.b())/2.*tan(angle);
               e3 = (double)((type_LLU)(e3*1000/5.)+1)*5/1000.;
        for(type_LLU n=2;n<3;n++) // number of bolts
        {
            Bolt Bo(HB_Diameter_Bolt.Get(i),2*n,bc,STANDART_Eurocode,1.25);
            if(n != 1)
            {
                Bo.Rectangle(p1,(n-1)*p1,2,n);
                Bo.Move     (0.000,e3+(n-1)*p1/2.);
            }
            else
            {
               Bo.SetXY(0,-a.b()-a.THK_PLATE()/2.,e3);
               Bo.SetXY(1,+a.b()+a.THK_PLATE()/2.,e3);
            }
            if(DEBUG) Bo.Printf();
            double ff = 0;
            // check plate //
            double WitdhPlate = a.THK_PLATE()+2*a.b()+2*0.005;
            Printf_CALC(
            Calc_FinPlate ( WitdhPlate, e1+(n-1)*p1+e3, WitdhPlate/2., 0,
                            THK_more(a.thk()), Pu, Py,
                            &Bo,
                            0, 1.0, 1.1,
                            0 , -ForceKnee, 0,
                            0 ,  0 , 0 ,
                            true , DEBUG),ff, DEBUG);
            double Nu_Rd;
            if(n == 1)
            {
                EN1993_1_8_Formula3_11(a.b()/2., Do, a.thk(), Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            else if(n == 2)
            {
                double Anet = a.A() - Do*a.thk();
                EN1993_1_8_Formula3_12(p1, Do, Anet, Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            else
            {
                double Anet = a.A() - Do*a.thk();
                EN1993_1_8_Formula3_13( p1, Do, Anet, Pu, Bo.GetGammaM2(), Nu_Rd, DEBUG);
            }
            Printf_CALC(fabs(ForceKnee/2.)/Nu_Rd,ff,DEBUG);
            if(ff < 1.00)
            {
                if(OUT)
                {
                    printf("Short output about connection:\n");
                    printf("Diameter of bolts: HM%2u\n",unsigned(Bo.Diameter_()*1e3));
                    printf("Number of bolts = %2u\n",n);
                    for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
                        printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                       Bo.positionBolts.Get(bb).x*1e3,
                                                       Bo.positionBolts.Get(bb).y*1e3);
                    Printf(bc);
                    printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",THK_more(a.thk())*1e3,WitdhPlate*1e3,(e1*2+0.010)*1e3);
                    printf("e1 = %.1f mm\n",e1*1e3);
                    printf("p1 = %.1f mm\n",p1*1e3);
                }
                return ff;
            }
        }
    }
    return 10;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
//double Calc_FinPlate_Is_row_one
//    (
//        IS is,
//        STAAD_FORCE_FULL Force,
//        bool OUT = true
//    )
//{
//
//    bool DEBUG = false;//true;//
//    if(DEBUG)printf("Start of Calc_FinPlate_Is_row_one\n");
//    double Py,Pu;
//    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
//
//    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
//    {
//        double Do = HB_Diameter_Bolt.Get(i)+0.003;
//        BOLT_CLASS bc = g5_6;
//        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
//        {
//            if(DEBUG)printf("0\n");
//            if(i_bc) bc = g5_6;
//            else     bc = g8_8;
//            if(HB_Diameter_Bolt.Get(i) == 0.016 )bc = g5_6;
//            double e1 = (double)((type_LLU)(2*Do*1000/5.)+1)*5/1000.;
//            double p1 = (double)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
//            double WitdhPlate = is.h()-2*is.tf()-2*is.r()-2*0.003;
//            double ThkPlate   = THK_more(is.tw());
//            type_LLU n = ((double)((WitdhPlate-2.*e1)/p1));// number of bolts
//            WitdhPlate = 2.*e1+p1*(n-1u);
//            if(n<2 || WitdhPlate > (is.h()-2*is.tf()-2*is.r()-2*0.003))
//            {
//                if(DEBUG)
//                {
//                    print_name("This is not good connection(Calc_FinPlate_Is_row_one)");
//                    printf("n = %u\n",n);
//                    printf("e1 = %.1f mm\n",e1*1e3);
//                    printf("p1 = %.1f mm\n",p1*1e3);
//                    printf("WitdhPlate = %.1f\t> %.1f\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r()-2*0.003)*1e3);
//                }
//                return 10;
//            }
//            Bolt Bo(HB_Diameter_Bolt.Get(i),n,bc,STANDART_Eurocode,1.25);
//            Bo.Rectangle(p1*(n-1),0,n,1);
//            Bo.Move(0,0.010+e1);
//            if(DEBUG) Bo.Printf();
//            double ff = 0;
//            // check plate //
//            Printf_CALC(
//            Calc_FinPlate ( WitdhPlate, e1*2+0.010, WitdhPlate/2., 0,
//                            ThkPlate, Pu, Py,
//                            &Bo,
//                            0, 1.0, 1.1,
//                            Force.V() , Force.N() , 0 ,
//                            0 , 0 , Force.V()*(e1+0.010)+Force.M() ,
//                            true , DEBUG),ff, DEBUG);
//            Printf_CALC(
//            Calc_FinPlate ( is.h(), e1*2+0.010, WitdhPlate/2., 0,
//                            is.tw(), Pu, Py,
//                            &Bo,
//                            0, 1.0, 1.1,
//                            Force.V() , Force.N() , 0 ,
//                            0 , 0 , Force.M() ,
//                            true , DEBUG),ff, DEBUG);
//            if(ff < 1.00)
//            {
//                if(OUT)
//                {
//                    printf("Short output about connection(Calc_FinPlate_Is_row_one):\n");
//                    printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
//                    printf("Number of bolts = %2u\n",n);
//                    for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
//                        printf("%2u - [%.1f,%.1f] mm\n",bb,
//                                                       Bo.positionBolts.Get(bb).x*1e3,
//                                                       Bo.positionBolts.Get(bb).y*1e3);
//                    Printf(bc);
//                    printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+0.010)*1e3);
//                    printf("e1 = %.1f mm\n",e1*1e3);
//                    printf("p1 = %.1f mm\n",p1*1e3);
//                }
//                return ff;
//            }
//        }
//    }
//    return 10;
//}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
//double Calc_FinPlate_Is_row_two
//    (
//        IS is,
//        STAAD_FORCE_FULL Force,
//        bool OUT = true
//    )
//{
//    bool DEBUG = false;//true;
//    double Py,Pu;
//    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
//
//    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
//    {
//        double Do = HB_Diameter_Bolt.Get(i)+0.003;
//        BOLT_CLASS bc = g5_6;
//        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
//        {
//            if(i_bc) bc = g5_6;
//            else     bc = g8_8;
//            if(HB_Diameter_Bolt.Get(i) == 0.016 )bc = g5_6;
//            double e1 = (double)((type_LLU)(2*Do*1000/5.)+1)*5/1000.;
//            double p1 = (double)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
//            double WitdhPlate = is.h()-2*is.tf()-2*is.r()-2*0.003;
//            double ThkPlate   = THK_more(is.tw());
//            type_LLU n = ((double)((WitdhPlate-2.*e1)/p1));// number of bolts
//            WitdhPlate = 2.*e1+p1*(n-1u);
//            if(n<2 || WitdhPlate > (is.h()-2*is.tf()-2*is.r()-2*0.003))
//            {
//                if(DEBUG)
//                {
//                    print_name("This is not good connection(Calc_FinPlate_Is_row_two)");
//                    printf("n = %u\n",n);
//                    printf("e1 = %.1f mm\n",e1*1e3);
//                    printf("p1 = %.1f mm\n",p1*1e3);
//                    printf("WitdhPlate = %.1f\t> %.1f\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r()-2*0.003)*1e3);
//                }
//                return 10;
//            }
//            Bolt Bo(HB_Diameter_Bolt.Get(i),n*2,bc,STANDART_Eurocode,1.25);
//            Bo.Rectangle(p1*(n-1),p1,n,2);
//            Bo.Move(0,0.010+e1+p1/2.);
//            if(DEBUG) Bo.Printf();
//            double ff = 0;
//            // check plate //
//            Printf_CALC(
//            Calc_FinPlate ( WitdhPlate, e1*2+p1+0.010, WitdhPlate/2., 0,
//                            ThkPlate, Pu, Py,
//                            &Bo,
//                            0, 1.0, 1.1,
//                            Force.V() , Force.N() , 0 ,
//                            0 , 0 , Force.V()*(e1+p1/2.+0.010)+Force.M() ,
//                            true , DEBUG),ff, DEBUG);
//            Printf_CALC(
//            Calc_FinPlate ( is.h(), e1*2+p1+0.010, WitdhPlate/2., 0,
//                            is.tw(), Pu, Py,
//                            &Bo,
//                            0, 1.0, 1.1,
//                            Force.V() , Force.N() , 0 ,
//                            0 , 0 , Force.M() ,
//                            true , DEBUG),ff, DEBUG);
//            if(ff < 1.00)
//            {
//                if(OUT)
//                {
//                    printf("Short output about connection(Calc_FinPlate_Is_row_two):\n");
//                    printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
//                    printf("Number of bolts = %2u\n",n);
//                    for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
//                        printf("%2u - [%.1f,%.1f] mm\n",bb,
//                                                       Bo.positionBolts.Get(bb).x*1e3,
//                                                       Bo.positionBolts.Get(bb).y*1e3);
//                    Printf(bc);
//                    printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+0.010)*1e3);
//                    printf("e1 = %.1f mm\n",e1*1e3);
//                    printf("p1 = %.1f mm\n",p1*1e3);
//                }
//                return ff;
//            }
//        }
//    }
//    return 10;
//}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
//double Calc_FinPlate_Is_column
//    (
//        IS is,
//        STAAD_FORCE_FULL Force,
//        bool OUT = true
//    )
//{
//
//    bool DEBUG = false;//true;
//    double Py,Pu;
//    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
//
//    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
//    {
//        double Do = HB_Diameter_Bolt.Get(i)+0.003;
//        BOLT_CLASS bc = g5_6;
//        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
//        {
//            if(i_bc) bc = g5_6;
//            else     bc = g8_8;
//            if(HB_Diameter_Bolt.Get(i) == 0.016 )bc = g5_6;
//            double e1 = (double)((type_LLU)(2*Do*1000/5.)+1)*5/1000.;
//            double p1 = (double)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
//            double WitdhPlate = is.h()-2*is.tf()-2*is.r()-2*0.003;
//            double ThkPlate   = THK_more(is.tw());
//            WitdhPlate = 2.*e1;
//            if(WitdhPlate > (is.h()-2*is.tf()-2*is.r()-2*0.003))
//            {
//                if(DEBUG)
//                {
//                    print_name("This is not good connection(Calc_FinPlate_Is_column)");
//                    printf("e1 = %.1f mm\n",e1*1e3);
//                    printf("p1 = %.1f mm\n",p1*1e3);
//                    printf("WitdhPlate = %.1f\t> %.1f\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r()-2*0.003)*1e3);
//                }
//                return 10;
//            }
//            for(type_LLU n=2;n<4;n++)
//            {
//                Bolt Bo(HB_Diameter_Bolt.Get(i),n,bc,STANDART_Eurocode,1.25);
//                Bo.Rectangle(0,p1*(n-1),1,n);
//                Bo.Move(0,0.010+e1+p1*(n-1)/2.);
//                if(DEBUG) Bo.Printf();
//                double ff = 0;
//                // check plate //
//                Printf_CALC(
//                Calc_FinPlate ( WitdhPlate,0.010+e1*2+p1*(n-1), WitdhPlate/2., 0,
//                                ThkPlate, Pu, Py,
//                                &Bo,
//                                0, 1.0, 1.1,
//                                Force.V() , Force.N() , 0 ,
//                                0 , 0 , Force.V()*(e1+p1*(n-1)/2.+0.010)+Force.M() ,
//                                true , DEBUG),ff, DEBUG);
//                Printf_CALC(
//                Calc_FinPlate ( is.h(), 0.010+e1*2+p1*(n-1), WitdhPlate/2., 0,
//                                is.tw(), Pu, Py,
//                                &Bo,
//                                0, 1.0, 1.1,
//                                Force.V() , Force.N() , 0 ,
//                                0 , 0 , Force.M() ,
//                                true , DEBUG),ff, DEBUG);
//                if(ff < 1.00)
//                {
//                    if(OUT)
//                    {
//                        printf("Short output about connection(Calc_FinPlate_Is_column):\n");
//                        printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
//                        printf("Number of bolts = %2u\n",n);
//                        for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
//                            printf("%2u - [%.1f,%.1f] mm\n",bb,
//                                                           Bo.positionBolts.Get(bb).x*1e3,
//                                                           Bo.positionBolts.Get(bb).y*1e3);
//                        Printf(bc);
//                        printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+0.010)*1e3);
//                        printf("e1 = %.1f mm\n",e1*1e3);
//                        printf("p1 = %.1f mm\n",p1*1e3);
//                    }
//                    return ff;
//                }
//            }
//        }
//    }
//    return 10;
//}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
enum ConnectionAutoEnum
{
    ConnectionAuto_Good , // good solution <100% //
    ConnectionAuto_Error, // error in input data //
    ConnectionAuto_Bad    // bad  solution >100% //
};
void Printf(ConnectionAutoEnum cae)
{
    switch(cae)
    {
        case(ConnectionAuto_Good ): printf("Good  output\n"); break;
        case(ConnectionAuto_Error): printf("Error output\n"); break;
        case(ConnectionAuto_Bad  ): printf("Bad   output\n"); break;
        default:
            print_name("Error in Printf(ConnectionAutoEnum cae)");
            FATAL();
    }
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
ConnectionAutoEnum Calc_FinPlate_Is_column
    (
        IS is,
        STAAD_FORCE_FULL Force,
        double &factor,
        bool OUT = true
    )
{
    bool DEBUG = false;//true;//
    double Py,Pu;
    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
    float DiatanceBetweenStructureAndIS = 0.015;

    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {
        double Do = HB_Diameter_Bolt.Get(i)+0.003;
        BOLT_CLASS bc;
        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
        {
            if(i_bc == 0) bc = g5_6;
            else          bc = g8_8;
            if(HB_Diameter_Bolt.Get(i) <= 0.016 )bc = g5_6;
            float e1 = (float)((type_LLU)(1.5*Do*1000/5.)+1)*5/1000.;
            float p1 = (float)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
            for(type_LLU n=2;n<5;n++)
            {
                if(DEBUG){line();printf("NUMBER OF BOLT = %u\n",n);line();}
                float WitdhPlate = min(is.h()-2*is.tf()-2*is.r()-2*0.003,2*e1+0.006);
                WitdhPlate = (float)((type_LLU)(WitdhPlate*1000/5.)-1)*5/1000.;
                if(DEBUG)
                {
                    printf("e1 = %.1f mm\n",e1*1e3);
                    printf("p1 = %.1f mm\n",p1*1e3);
                    printf("WitdhPlate = [%.1f\t; %.1f]\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r())*1e3);
                }
                if(WitdhPlate+1e-6<2*e1)
                    {if(DEBUG)printf("BREAK");break;}//return ConnectionAuto_Error;
                float ThkPlate   = THK_more(((float)n-0.75)*is.tw()+0.001);
                Bolt Bo(HB_Diameter_Bolt.Get(i),n,bc,STANDART_Eurocode,1.25);
                Bo.Rectangle(0,p1*(n-1),1,n);
                Bo.Move(0,DiatanceBetweenStructureAndIS+e1+p1*(n-1)/2.);
                if(DEBUG) Bo.Printf();
                factor = 0;
                // check plate //
                Printf_CALC(
                Calc_FinPlate ( WitdhPlate,DiatanceBetweenStructureAndIS+e1*2+p1*(n-1), WitdhPlate/2., 0,
                                ThkPlate, Pu, Py,
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.V()*(e1+p1*(n-1)/2.+DiatanceBetweenStructureAndIS)+Force.M() ,
                                true , DEBUG),factor, DEBUG);
                Printf_CALC(
                Calc_FinPlate ( is.h() ,DiatanceBetweenStructureAndIS+e1*2+p1*(n-1), is.h(), 0,
                                is.tw(), Pu, Py,//min(is.h(),4*is.tw()+0.040)
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.M() ,
                                true , DEBUG),factor, DEBUG);
                //check profile
                if(factor < 1.00)
                {
                    if(OUT)
                    {
                        printf("Short output about connection(Calc_FinPlate_Is_column):\n");
                        printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
                        printf("Number of bolts = %2u\n",n);
                        for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
                            printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                           Bo.positionBolts.Get(bb).x*1e3,
                                                           Bo.positionBolts.Get(bb).y*1e3);
                        Printf(bc);
                        printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+p1*(n-1)+DiatanceBetweenStructureAndIS)*1e3);
                        printf("e1 = %.1f mm\n",e1*1e3);
                        printf("p1 = %.1f mm\n",p1*1e3);
                    }
                    return ConnectionAuto_Good;
                }
            }
        }
    }
    return ConnectionAuto_Bad;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
ConnectionAutoEnum Calc_FinPlate_Is_row_one
    (
        IS is,
        STAAD_FORCE_FULL Force,
        double &factor,
        bool OUT = true
    )
{
    bool DEBUG = false;//true;//
    double Py,Pu;
    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
    float DiatanceBetweenStructureAndIS = 0.015;
    if(HB_Diameter_Bolt.GetSize() < 1)
    {
        print_name("Please add HB_Diameter_Bolt in Calc_FinPlate_Is_row_one");
        printf("Size of HB_Diameter_Bolt = %u\n",HB_Diameter_Bolt.GetSize());
        Printf(HB_Diameter_Bolt);
        FATAL();
    }

    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {
        double Do = HB_Diameter_Bolt.Get(i)+0.003;
        BOLT_CLASS bc;
        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
        {
            if(DEBUG){line();printf("HB_Diameter_Bolt.Get(i) = %f\n",HB_Diameter_Bolt.Get(i)*1e3);line();}
            if(i_bc == 0) bc = g5_6;
            else          bc = g8_8;
            if(HB_Diameter_Bolt.Get(i) <= 0.016 )bc = g5_6;
            float e1 = (float)((type_LLU)(1.5*Do*1000/5.)+1)*5/1000.;
            float p1 = (float)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
            for(type_LLU n=2;n<5;n++)
            {
                if(DEBUG){line();printf("NUMBER OF BOLT = %u\n",n);line();}
                float WitdhPlate = min(is.h()-2*is.tf()-2*is.r()-2*0.003,2*e1+(n-1)*p1+0.006);
                WitdhPlate = (float)((type_LLU)(WitdhPlate*1000/5.)-1)*5/1000.;
                if(DEBUG)
                {
                    printf("e1 = %.1f mm\n",e1*1e3);
                    printf("p1 = %.1f mm\n",p1*1e3);
                    printf("WitdhPlate = [%.1f\t; %.1f]\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r())*1e3);
                }
                if(WitdhPlate+1e-6<2*e1+(n-1)*p1)
                    {if(DEBUG)printf("BREAK");break;}//return ConnectionAuto_Error;
                float ThkPlate   = THK_more(is.tw());//(((float)n-0.75)*is.tw()+0.001);
                Bolt Bo(HB_Diameter_Bolt.Get(i),1*n,bc,STANDART_Eurocode,1.25);
                Bo.Rectangle(p1*(n-1),p1,n,1);
                Bo.Move(0,DiatanceBetweenStructureAndIS+e1);
                if(DEBUG) Bo.Printf();
                factor = 0;
                // check plate //
                Printf_CALC(
                Calc_FinPlate ( WitdhPlate,DiatanceBetweenStructureAndIS+e1*2, WitdhPlate/2., 0,
                                ThkPlate, Pu, Py,
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.V()*(e1+DiatanceBetweenStructureAndIS)+Force.M() ,
                                true , DEBUG),factor, DEBUG);
                Printf_CALC(
                Calc_FinPlate ( is.h() ,DiatanceBetweenStructureAndIS+e1*2, is.h()/2., 0,
                                is.tw(), Pu, Py,//min(is.h(),4*is.tw()+0.040)
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.M() ,
                                true , DEBUG),factor, DEBUG);
                //check profile
                if(factor < 1.00)
                {
                    if(OUT)
                    {
                        printf("Short output about connection(Calc_FinPlate_Is_row_one):\n");
                        printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
                        printf("Number of bolts = %2u\n",n);
                        for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
                            printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                           Bo.positionBolts.Get(bb).x*1e3,
                                                           Bo.positionBolts.Get(bb).y*1e3);
                        Printf(bc);
                        printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+p1*(n-1)+DiatanceBetweenStructureAndIS)*1e3);
                        printf("e1 = %.1f mm\n",e1*1e3);
                        printf("p1 = %.1f mm\n",p1*1e3);
                    }
                    return ConnectionAuto_Good;
                }
            }
        }
    }
    return ConnectionAuto_Bad;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
ConnectionAutoEnum Calc_FinPlate_Is_row_two
    (
        IS is,
        STAAD_FORCE_FULL Force,
        double &factor,
        bool OUT = true
    )
{
    bool DEBUG = false;//true;//
    double Py,Pu;
    SNiP("C235",max(is.tw(),is.tf()),Py,Pu);
    float DiatanceBetweenStructureAndIS = 0.015;

    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {
        double Do = HB_Diameter_Bolt.Get(i)+0.003;
        BOLT_CLASS bc;
        for(type_LLU i_bc = 0;i_bc < 2;i_bc++)
        {
            if(DEBUG){line();printf("HB_Diameter_Bolt.Get(i) = %f\n",HB_Diameter_Bolt.Get(i)*1e3);line();}
            if(i_bc == 0) bc = g5_6;
            else          bc = g8_8;
            if(HB_Diameter_Bolt.Get(i) <= 0.016 )bc = g5_6;
            float e1 = (float)((type_LLU)(1.5*Do*1000/5.)+1)*5/1000.;
            float p1 = (float)((type_LLU)(2.5*Do*1000/5.)+1)*5/1000.;
            for(type_LLU n=2;n<5;n++)
            {
                if(DEBUG){line();printf("NUMBER OF BOLT = %u\n",n);line();}
                float WitdhPlate = min(is.h()-2*is.tf()-2*is.r()-2*0.003,2*e1+(n-1)*p1+0.006);
                WitdhPlate = (float)((type_LLU)(WitdhPlate*1000/5.)-1)*5/1000.;
                if(DEBUG)
                {
                    printf("e1 = %.1f mm\n",e1*1e3);
                    printf("p1 = %.1f mm\n",p1*1e3);
                    printf("WitdhPlate = [%.1f\t; %.1f]\n",WitdhPlate*1e3,(is.h()-2*is.tf()-2*is.r())*1e3);
                }
                if(WitdhPlate+1e-6<2*e1+(n-1)*p1)
                    {if(DEBUG)printf("BREAK");break;}//return ConnectionAuto_Error;
                float ThkPlate   = THK_more(1.9*is.tw());//(((float)n-0.75)*is.tw()+0.001);
                Bolt Bo(HB_Diameter_Bolt.Get(i),2*n,bc,STANDART_Eurocode,1.25);
                Bo.Rectangle(p1*(n-1),p1,n,2);
                Bo.Move(0,DiatanceBetweenStructureAndIS+e1+p1/2.);
                if(DEBUG) Bo.Printf();
                factor = 0;
                // check plate //
                Printf_CALC(
                Calc_FinPlate ( WitdhPlate,DiatanceBetweenStructureAndIS+e1*2+p1, WitdhPlate/2., 0,
                                ThkPlate, Pu, Py,
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.V()*(e1+p1+DiatanceBetweenStructureAndIS)+Force.M() ,
                                true , DEBUG),factor, DEBUG);
                Printf_CALC(
                Calc_FinPlate ( is.h() ,DiatanceBetweenStructureAndIS+e1*2+p1, is.h()/2., 0,
                                is.tw(), Pu, Py,//min(is.h(),4*is.tw()+0.040)
                                &Bo,
                                0, 1.0, 1.1,
                                Force.V() , Force.N() , 0 ,
                                0 , 0 , Force.M() ,
                                true , DEBUG),factor, DEBUG);
                //check profile
                if(factor < 1.00)
                {
                    if(OUT)
                    {
                        printf("Short output about connection(Calc_FinPlate_Is_row_two):\n");
                        printf("Diameter of bolts: HM%2u\n",(type_LLU)(Bo.Diameter_()*1e3));
                        printf("Number of bolts = %2u\n",n);
                        for(type_LLU bb = 0;bb<Bo.positionBolts.GetSize();bb++)
                            printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                           Bo.positionBolts.Get(bb).x*1e3,
                                                           Bo.positionBolts.Get(bb).y*1e3);
                        Printf(bc);
                        printf("Size of Plate(%.1f x %.1f x %.1f)mm\n",ThkPlate*1e3,WitdhPlate*1e3,(e1*2+p1*(n-1)+DiatanceBetweenStructureAndIS)*1e3);
                        printf("e1 = %.1f mm\n",e1*1e3);
                        printf("p1 = %.1f mm\n",p1*1e3);
                    }
                    return ConnectionAuto_Good;
                }
            }
        }
    }
    return ConnectionAuto_Bad;
}
///////////////////////////////////////////
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
//                                       //
///////////////////////////////////////////
double Calc_MomentConnection
     (
        IS*       I,
        STAAD_FORCE_FULL Force,
        bool THK_END_PLATE_ANALYSYS = false,
        bool OUT = true
     )
{
    bool DEBUG = true;//false;//
    if(DEBUG) printf("DEBUG OUTPUT\n");
    for(type_LLU i=0;i<HB_Diameter_Bolt.GetSize();i++)
    {

Calc_MomentConnection1:

        if(DEBUG) printf("DEBUG: Diameter of bolts: %.0f\n",HB_Diameter_Bolt[i]*1000);
        BOLT_CLASS bc = g5_6;
        if(HB_Diameter_Bolt.Get(i)>0.016)
            bc = g8_8;
        if(HB_Diameter_Bolt.Get(i)>0.031)
        {
            if(DEBUG) print_name("NOT CONNECTION");
            break;//return 10;
        }


        double A,M,L;ANUREV(HB_Diameter_Bolt.Get(i), A, M, L, DEBUG);
        double X_bolts = max(I->b()/2.,I->tw()+2*M);
        X_bolts = (((type_LLU)(X_bolts*1000./10.))+1)*10./1000.;
        double Y_bolts = min(I->h()-2.*I->tf()-2*M,I->h()-2*1.2*(HB_Diameter_Bolt.Get(i)+0.003));
        Y_bolts = (((type_LLU)(Y_bolts*1000./10.))-1)*10./1000.;
        double X_plate = max(X_bolts+4.*(HB_Diameter_Bolt.Get(i)+0.003),I->b());
        X_plate = (((type_LLU)(X_plate*1000./10.))+1)*10./1000.;
        double Y_plate = I->h();

        if((I->h()-Y_bolts)/2. < max(1.2*(HB_Diameter_Bolt.Get(i)+0.0029),I->tf()+M) ||
           Y_bolts<2.200*(HB_Diameter_Bolt.Get(i)+0.003) ||
           Y_bolts<A)
        {
            if(DEBUG)
            {
                print_name("NOT CONNECTION e2");
                I->Printf();
                printf("Y_bolts = %.1f mm\n",Y_bolts*1e3);
                printf("M = %.1f mm\n",M*1e3);
                printf("HB_Diameter_Bolt.Get(i) = %.1f mm\n",HB_Diameter_Bolt.Get(i)*1e3);
            }
            break;
        }
        if(Y_bolts+2*M+2*I->tf() >= I->h())
        {
            if(DEBUG) print_name("Don`t find good connection by Calc_MomentConnection");
            i++; goto Calc_MomentConnection1;
        }

        Gusset gNULL;
        Bolt      B(HB_Diameter_Bolt.Get(i),4,bc,STANDART_Eurocode,1.25);
                  B.Rectangle(X_bolts,Y_bolts,2,2);
                  B.Move     (0.00000,I->h()/2.00);
        if(DEBUG) B.Printf();

        double thk_plate = THK_more(max(I->tf(),I->tw()));
        if(!THK_END_PLATE_ANALYSYS)thk_plate = THK_more(HB_Diameter_Bolt.Get(i));
        for(;thk_plate<THK_more(HB_Diameter_Bolt.Get(i))+0.0001;thk_plate+=0.002)
        {
            thk_plate = THK_more(thk_plate);
            if(DEBUG)printf("thk_plate = %.2f mm\n",thk_plate*1e3);
            if(DEBUG)printf("X_plate   = %.2f mm\n",X_plate*1e3);
            if(DEBUG)printf("Y_plate   = %.2f mm\n",Y_plate*1e3);
            EndPlate  E(thk_plate,X_plate,Y_plate,0.0,"C235");
            if(DEBUG) E.Printf();
            double factor = Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&gNULL,&gNULL,DEBUG);
            if(factor < 1.00)
            {
                // CORRECT OUTPUT //
                //if(DEBUG)Calc_MomentConnection (I, &E, &B, 0,0, 1.0,1.1, Force,&gNULL,&gNULL,&gNULL,&gNULL,true);
                if(OUT)
                {
                    printf("Short output about moment connection:\n");
                    printf("Section  : %s\n",I->name);
                    printf("End plate: %.1f x %.1f x %.1f mm\n",E.THK*1e3,E.X*1e3,(E.Yplus-E.Yminus)*1e3);
                    printf("Diameter of bolts: HM%.0f\n",B.Diameter_()*1e3);
                    printf("Number of bolts = %2u\n",B.positionBolts.GetSize());
                    for(type_LLU bb = 0;bb<B.positionBolts.GetSize();bb++)
                        printf("%2u - [%.1f,%.1f] mm\n",bb,
                                                        B.positionBolts.Get(bb).x*1e3,
                                                        B.positionBolts.Get(bb).y*1e3);
                    Printf(bc);
                }
                return factor;
            }
        }
    }
    return 10;
}

