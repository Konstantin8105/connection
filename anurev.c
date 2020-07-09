
///
/// FROM BOOK Anurev
///

void ANUREV(double DiameterOfBolt, double& A, double& M, double& L, bool OUT = true)
{
    if (DiameterOfBolt == 0.010) {
        A = 0.034;
        M = 0.017;
        L = 0.052;
    }
    else if (DiameterOfBolt == 0.012) {
        A = 0.036;
        M = 0.019;
        L = 0.060;
    }
    else if (DiameterOfBolt == 0.016) {
        A = 0.048;
        M = 0.025;
        L = 0.078;
    }
    else if (DiameterOfBolt == 0.020) {
        A = 0.058;
        M = 0.030;
        L = 0.098;
    }
    else if (DiameterOfBolt == 0.024) {
        A = 0.068;
        M = 0.036;
        L = 0.110;
    }
    else if (DiameterOfBolt == 0.030) {
        A = 0.090;
        M = 0.045;
        L = 0.140;
    }
    else if (DiameterOfBolt == 0.036) {
        A = 0.105;
        M = 0.052;
        L = 0.160;
    }
    else if (DiameterOfBolt == 0.048) {
        A = 0.140;
        M = 0.070;
        L = 0.210;
    }
    else {
        print_name("Please add information about this bolts from Anuriev");
        FATAL();
    }
    if (OUT) {
        printf("Distances for bolt HM%.0f by Anurev:\n", DiameterOfBolt * 1e3);
        printf("A = %.1f mm\tM = %.1f mm\tL = %.1f mm\n", A * 1e3, M * 1e3, L * 1e3);
    }
}

bool ANUREV_MinDistance(double MinDistance, double DiameterOfBolt, bool OUT = true)
{
    bool ex = false;
    double A, M, L;
    ANUREV(DiameterOfBolt, A, M, L, OUT);
    if (MinDistance >= A - 1e-6)
        ex = true;
    if (OUT) {
        printf("Check minimal distance by book Anurev:");
        if (ex)
            printf("OK\n");
        else {
            printf("None\n");
            printf("Diameter Of Bolt = %.1f mm\n", DiameterOfBolt * 1e3);
            printf("A = %.1f mm\n", A * 1e3);
            FATAL();
        }
    }
    return ex;
}

void ANUREV_Schema2(double P, double l, double J, double& Re, double& M, double& Vmax, double& Qmax, bool OUT = true)
{
    Re = P;
    M = P * l;
    Vmax = P * pow(l, 3.) / (3. * CONST_ELACITY * J);
    Qmax = P * pow(l, 2.) / (2. * CONST_ELACITY * J);
    if (OUT) {
        printf("Beam by schema 2 book Anurev\n");
        printf("P = %.1f N\n", P);
        printf("L = %.1f m\n", l);
        printf("J = %.1f cm4\n", J * 1.e8);
        printf("P = %.1f N   - shear force\n", P);
        printf("R = %.1f N   - reaction shear force\n", P / 2.);
        printf("M = %.1f N*m - moment\n", M);
        printf("Q = %.3f deg - angle of rotation\n", GRADIANS(Qmax));
        printf("v = %.2f mm - deflection\n\n", Vmax * 1.e3);
    }
}

void ANUREV_Schema5(double P, double l, double J, double& Re, double& M, double& Vmax, bool OUT = true)
{
    Re = P / 2.;
    M = P * l / 4.;
    Vmax = 1. / 48. * P * pow(l, 3.) / (CONST_ELACITY * J);
    if (OUT) {
        printf("Beam by schema 5 book Anurev\n");
        printf("P = %.1f N\n", P);
        printf("L = %.1f m\n", l);
        printf("J = %.1f cm4\n", J * 1.e8);
        printf("R = %.1f N   - reaction shear force\n", Re);
        printf("M = %.1f N*m - moment\n", M);
        printf("v = %.2f mm - deflection\n\n", Vmax * 1.e3);
    }
}

void ANUREV_Schema6(double q, double l, double J, double& Re, double& M, double& Vmax, bool OUT = true)
{
    Re = q * l / 2.;
    M = 1. / 8. * q * pow(l, 2.);
    Vmax = 5. / 384. * q * pow(l, 4.) / (CONST_ELACITY * J);
    if (OUT) {
        printf("Beam by schema 6 book Anurev\n");
        printf("q = %.1f N/m\n", q);
        printf("L = %.1f m\n", l);
        printf("J = %.1f cm4\n", J * 1.e8);
        printf("R = %.1f N   - reaction shear force\n", Re);
        printf("M = %.1f N*m - moment\n", M);
        printf("v = %.2f mm - deflection\n\n", Vmax * 1.e3);
    }
}

void ANUREV_Schema11(double P, double l, double J, double& Re, double& M, double& Vmax, double& Qmax, bool OUT = true)
{
    Re = max(11. / 16., 5. / 16.) * P;
    M = 3. / 16. * P * l;
    Vmax = 0.0093 * P * pow(l, 3.) / (CONST_ELACITY * J);
    Qmax = P * pow(l, 2.) / (32. * CONST_ELACITY * J);
    if (OUT) {
        printf("Beam by schema 2 book Anurev\n");
        printf("P = %.1f N\n", P);
        printf("L = %.1f m\n", l);
        printf("J = %.1f cm4\n", J * 1.e8);
        printf("P = %.1f N   - shear force\n", P);
        printf("R = %.1f N   - reaction shear force\n", P / 2.);
        printf("M = %.1f N*m - moment\n", M);
        printf("Q = %.3f deg - angle of rotation\n", GRADIANS(Qmax));
        printf("v = %.2f mm - deflection\n\n", Vmax * 1.e3);
    }
}

void ANUREV_Schema15(double q, double l, double J, double& Re, double& M, double& Vmax, bool OUT = true)
{
    Re = q * l / 2.;
    M = 1. / 12. * q * pow(l, 2.);
    Vmax = 1. / 384. * q * pow(l, 4.) / (CONST_ELACITY * J);
    if (OUT) {
        printf("Beam by schema 15 book Anurev\n");
        printf("q = %.1f N/m\n", q);
        printf("L = %.1f m\n", l);
        printf("J = %.1f cm4\n", J * 1.e8);
        printf("R = %.1f N   - reaction shear force\n", Re);
        printf("M = %.1f N*m - moment\n", M);
        printf("v = %.2f mm - deflection\n\n", Vmax * 1.e3);
    }
}

