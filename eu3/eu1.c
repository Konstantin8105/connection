void EN1993_1_1_Formula6_18(double Fy, double A, double gamma_M0,
    double& Vpl_Rd, bool OUT = true)
{
    ///
    /// Shear
    ///
    Vpl_Rd = (1 / sqrt(3)) * Fy * A / gamma_M0;
    if (OUT) {
        printf("Calculation: Shear Yielding Verification of the Plate\n");
        printf("Shear Yielding Verification of the Plate\n");
        printf("Ved <= Vpl_Rd, see formula (6.17) EN1993-1-1\n");
        printf("Vpl_Rd = Fy*Av/gamma_M0/sqrt(3)\n");
        printf("Av = %.1f mm2\n", A * 1e6);
        printf("Vpl_Rd = %.1f * %.1f /%.2f/%.2f\n", Fy * 1e-6, A * 1e6, gamma_M0, sqrt(3));
        printf("Vpl_Rd = %.1f kN\n", Vpl_Rd * 1e-3);
    }
}

void EN1993_1_1_Formula6_6(double A, double Fy, double gamma_M0, double& Npl_Rd, bool OUT = true)
{
    Npl_Rd = A * Fy / gamma_M0;
    if (OUT) {
        printf("Calculate Npl_Rd by formula (6.6) EN1993-1-1\n");
        printf("Npl_Rd = A*Fy/gamma_M0\n");
        printf("Npl_Rd = %.1f*%.1f/%.3f\n", A * 1e6, Fy * 1e-6, gamma_M0);
        printf("Npl_Rd = %.1f kN\n", Npl_Rd * 1e-3);
    }
}

double EN1993_1_1_Formula6_12(double Med, double Mc_Rd, bool OUT = true)
{
    double factor = Med / Mc_Rd;
    if (OUT) {
        printf("Calculate by formula (6.12) EN1993-1-1\n");
        printf("Med/Mc_Rd <= 1\n");
        printf("Med   = %.1f kN*m\n", Med * 1e-3);
        printf("Mc_Rd = %.1f kN*m\n", Mc_Rd * 1e-3);
    }
    return factor;
}

void EN1993_1_1_Formula6_13(double Wpl, double Fy, double gamma_M0, double& Mpl_Rd, bool OUT = true)
{
    Mpl_Rd = Wpl * Fy / gamma_M0;
    if (OUT) {
        printf("Calculate Mpl_Rd by formula (6.13) EN1993-1-1\n");
        printf("Mpl_Rd = W*Fy/gamma_M0\n");
        printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wpl * 1e6, Fy * 1e-6, gamma_M0);
        printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd * 1e-3);
    }
}

void EN1993_1_1_Formula6_14(double Wei, double Fy, double gamma_M0, double& Mpl_Rd, bool OUT = true)
{
    Mpl_Rd = Wei * Fy / gamma_M0;
    if (OUT) {
        printf("Calculate Mpl_Rd by formula (6.14) EN1993-1-1\n");
        printf("Mpl_Rd = W*Fy/gamma_M0\n");
        printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wei * 1e6, Fy * 1e-6, gamma_M0);
        printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd * 1e-3);
    }
}

void EN1993_1_1_Formula6_15(double Wei, double Fy, double gamma_M0, double& Mpl_Rd, bool OUT = true)
{
    Mpl_Rd = Wei * Fy / gamma_M0;
    if (OUT) {
        printf("Calculate Mpl_Rd by formula (6.15) EN1993-1-1\n");
        printf("Mpl_Rd = W*Fy/gamma_M0\n");
        printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wei * 1e6, Fy * 1e-6, gamma_M0);
        printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd * 1e-3);
    }
}

enum CrossSection { CrossSection_Class1,
    CrossSection_Class2,
    CrossSection_Class3,
    CrossSection_Class4,
};
CrossSection EN1993_1_1_Table_5_2_tube(double d, double t, double Fy, bool OUT = true)
{
    CrossSection CS;
    if (Fy < 100e6 || t > d) {
        print_name("Error in EN1993_1_1_Table_5_2_tube");
        FATAL();
    }
    double e = sqrt(235.e6 / Fy);
    if (d / t <= 50 * pow(e, 2.))
        CS = CrossSection_Class1;
    else if (d / t <= 70 * pow(e, 2.))
        CS = CrossSection_Class2;
    else if (d / t <= 90 * pow(e, 2.))
        CS = CrossSection_Class3;
    else {
        print_name("Error in EN1993_1_1_Table_5_2_tube. See EN1993-1-6");
        FATAL();
    }
    if (OUT) {
        printf("Tube %.1fx%.1f mm - by table 5.2 EN1993-1-1 class ", d * 1e3, t * 1e3);
        if (CS == CrossSection_Class1)
            printf("1\n");
        else if (CS == CrossSection_Class2)
            printf("2\n");
        else if (CS == CrossSection_Class3)
            printf("3\n");
    }
    return CS;
}
CrossSection EN1993_1_1_Table_5_2_Isection(double h, double b, double tw, double tf, double Fy, bool OUT = true)
{
    CrossSection CS1, CS2;
    if (Fy < 100e6 || tf < tw || b < tf || b < tw) {
        print_name("Error in EN1993_1_1_Table_5_2_Isection");
        FATAL();
    }
    double e = sqrt(235.e6 / Fy);
    double c, t;
    c = h - 2. * tf;
    t = tw;
    if (c / t <= 33. * e)
        CS1 = CrossSection_Class1;
    else if (c / t <= 38. * e)
        CS1 = CrossSection_Class2;
    else if (c / t <= 42. * e)
        CS1 = CrossSection_Class3;
    else {
        print_name("Error in EN1993_1_1_Table_5_2_Isection. See EN1993-1-6");
        FATAL();
    }
    c = b / 2. - 2. * tw;
    t = tf;
    if (c / t <= 09. * e)
        CS2 = CrossSection_Class1;
    else if (c / t <= 10. * e)
        CS2 = CrossSection_Class2;
    else if (c / t <= 14. * e)
        CS2 = CrossSection_Class3;
    else {
        print_name("Error in EN1993_1_1_Table_5_2_Isection. See EN1993-1-6");
        FATAL();
    }
    CrossSection CS = CS1;
    if (CS > CS2)
        CS = CS2;
    if (OUT) {
        printf("I-section - by table 5.2 EN1993-1-1 class ");
        if (CS == CrossSection_Class1)
            printf("1\n");
        else if (CS == CrossSection_Class2)
            printf("2\n");
        else if (CS == CrossSection_Class3)
            printf("3\n");
    }
    return CS;
}

// Buckling curve
enum BUCKLING_CURVE { BUCKLING_CURVE_Ao,
    BUCKLING_CURVE_A,
    BUCKLING_CURVE_B,
    BUCKLING_CURVE_C,
    BUCKLING_CURVE_D };

void Printf(BUCKLING_CURVE Bc)
{
    printf("Buckling curve: ");
    switch (Bc) {
    case (BUCKLING_CURVE_Ao):
        printf("ao");
        break;
    case (BUCKLING_CURVE_A):
        printf("a ");
        break;
    case (BUCKLING_CURVE_B):
        printf("b ");
        break;
    case (BUCKLING_CURVE_C):
        printf("c ");
        break;
    case (BUCKLING_CURVE_D):
        printf("d ");
        break;
    default:
        print_name("ERROR in EN1991_1_1_Table6_1\n");
        FATAL();
    }
    printf("\n");
}

double EN1991_1_1_Table6_1_Alpha(BUCKLING_CURVE Bc, bool OUT = true)
{
    double alpha = -1;
    switch (Bc) {
    case (BUCKLING_CURVE_Ao):
        alpha = 0.13;
        break;
    case (BUCKLING_CURVE_A):
        alpha = 0.21;
        break;
    case (BUCKLING_CURVE_B):
        alpha = 0.34;
        break;
    case (BUCKLING_CURVE_C):
        alpha = 0.49;
        break;
    case (BUCKLING_CURVE_D):
        alpha = 0.76;
        break;
    default:
        print_name("ERROR in EN1991_1_1_Table6_1\n");
        FATAL();
    }
    if (OUT) {
        Printf(Bc);
        printf("Imperfection factor alpha = %.2f by Table 6.1 EN1993-1-1\n", alpha);
    }
    return alpha;
}

void EN1993_1_1_table_6_2_Isection(double tf, char* Category, BUCKLING_CURVE& YY, BUCKLING_CURVE& ZZ, bool OUT = true)
{
    if (tf > 0.100 || tf < 0.00) {
        print_name("Please check - stupid tf in EN1993_1_1_table_6_2_Isection");
        FATAL();
    }
    if (strcmp(Category, "C235") == 0 || strcmp(Category, "C275") == 0 || strcmp(Category, "C355") == 0 || strcmp(Category, "C420") == 0) {
        if (tf < 0.040) {
            YY = BUCKLING_CURVE_B;
            ZZ = BUCKLING_CURVE_C;
        }
        else {
            YY = BUCKLING_CURVE_C;
            ZZ = BUCKLING_CURVE_D;
        }
    }
    else {
        print_name("Please add this material to EN1993_1_1_table_6_2_Isection");
        FATAL();
    }
    if (OUT) {
        printf("I-section by table 6.2:\n");
        printf("for axes Y-Y :");
        Printf(YY);
        printf("\n");
        printf("for axes Z-Z :");
        Printf(ZZ);
        printf("\n");
    }
}

void EN1993_1_1_Formula6_47(double Xi, double A, double Fy, double gamma_M1, CrossSection beam, double& Nb_Rd, bool OUT = true)
{
    if (beam != CrossSection_Class1 && beam != CrossSection_Class2 && beam != CrossSection_Class3) {
        printf("Please use another formula. EN1993_1_1_Formula6_47");
        FATAL();
    }
    Nb_Rd = Xi * A * Fy / gamma_M1;
    if (OUT) {
        printf("Calculate Nb_Rd by formula (6.47) EN1993-1-1\n");
        printf("Nb_Rd = Xi*A*Fy/gamma_M1\n");
        printf("Nb_Rd = %.1f*%.1f*%.1f/%.3f\n", Xi, A * 1e6, Fy * 1e-6, gamma_M1);
        printf("Nb_Rd = %.1f kN\n", Nb_Rd * 1e-3);
    }
}

double EN1993_1_1_Formula6_49_Xi(double Lambda_, BUCKLING_CURVE Bc, bool OUT = true)
{
    double Fi = 0.5 * (1 + EN1991_1_1_Table6_1_Alpha(Bc, OUT) * (Lambda_ - 0.2) + pow(Lambda_, 2.));
    double Xi = 1. / (Fi + sqrt(Fi * Fi - Lambda_ * Lambda_));
    if (Xi > 1.0)
        Xi = 1.0;
    if (OUT) {
        printf("Fi = %f\n", Fi);
        printf("Xi = %f\n", Xi);
    }
    return Xi;
}

void EN1993_1_1_Punkt6_3_1_2_Lambda_c1(double A, double Fy, double Ncr, double& Lambda_, bool OUT = true)
{
    Lambda_ = sqrt(A * Fy / Ncr);
    if (OUT) {
        print_name("ewrwerwe add good view");
        printf("Lambda_ = %.2f\n", Lambda_);
    }
}
double EN1993_1_1_Punkt6_3_1_2_Lambda_c2(double Lcr, double i, double Lambda1, bool OUT = true)
{
    double Lambda_ = Lcr / i / Lambda1;
    if (OUT) {
        printf("Lambda_ = %.2f\n", Lambda_);
    }
    return Lambda_;
}

double EN1993_1_1_Punkt6_50_Lambda1(double E, double Fy, bool OUT = true)
{
    double Lambda1 = CONST_M_PI * sqrt(E / Fy);
    if (OUT) {
        printf("Lambda1 = PI*sqrt(E/Fy) = PI*sqrt(%.2f/%.2f) = %f\n", E, Fy, Lambda1);
    }
    return Lambda1;
}

void EN1993_1_1_Ncr(double ELASITY, double J_moment_inetria, double L, double& Ncr, bool OUT = true)
{
    Ncr = pow(CONST_M_PI / L, 2.0) * ELASITY * J_moment_inetria;
    if (OUT) {
        printf("Ncr = (PI/L)^2*E*J\n");
        printf("Ncr = (PI/%.1f)^2*%.1f*%.1f\n", L * 1e3, ELASITY * 1e-6, J_moment_inetria * 1e8);
        printf("Ncr = %.1f kN\n", Ncr * 1e-3);
    }
}

