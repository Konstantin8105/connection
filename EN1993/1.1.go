package eu311

import (
	"fmt"
	"math"
)

type CrossSection int

const (
	CrossSection_Class1 CrossSection = 0
	CrossSection_Class2              = 1
	CrossSection_Class3              = 2
	CrossSection_Class4              = 3
)

type BUCKLING_CURVE int

const (
	BUCKLING_CURVE_Ao BUCKLING_CURVE = 0
	BUCKLING_CURVE_A                 = 1
	BUCKLING_CURVE_B                 = 2
	BUCKLING_CURVE_C                 = 3
	BUCKLING_CURVE_D                 = 4
)

func Formula6_18(Fy float64, A float64, gamma_M0 float64, Vpl_Rd float64, out bool) {
	Vpl_Rd = 1 / math.Sqrt(3) * Fy * A / gamma_M0
	if out {
		fmt.Printf("Calculation: Shear Yielding Verification of the Plate\n")
		fmt.Printf("Shear Yielding Verification of the Plate\n")
		fmt.Printf("Ved <= Vpl_Rd, see formula (6.17) EN1993-1-1\n")
		fmt.Printf("Vpl_Rd = Fy*Av/gamma_M0/sqrt(3)\n")
		fmt.Printf("Av = %.1f mm2\n", A*1e+06)
		fmt.Printf("Vpl_Rd = %.1f * %.1f /%.2f/%.2f\n", Fy*1e-06, A*1e+06, gamma_M0, math.Sqrt(3))
		fmt.Printf("Vpl_Rd = %.1f kN\n", Vpl_Rd*0.001)
	}
}

func Formula6_6(A float64, Fy float64, gamma_M0 float64, Npl_Rd float64, out bool) {
	Npl_Rd = A * Fy / gamma_M0
	if out {
		fmt.Printf("Calculate Npl_Rd by formula (6.6) EN1993-1-1\n")
		fmt.Printf("Npl_Rd = A*Fy/gamma_M0\n")
		fmt.Printf("Npl_Rd = %.1f*%.1f/%.3f\n", A*1e+06, Fy*1e-06, gamma_M0)
		fmt.Printf("Npl_Rd = %.1f kN\n", Npl_Rd*0.001)
	}
}

func Formula6_12(Med float64, Mc_Rd float64, out bool) float64 {
	var factor float64 = Med / Mc_Rd
	if out {
		fmt.Printf("Calculate by formula (6.12) EN1993-1-1\n")
		fmt.Printf("Med/Mc_Rd <= 1\n")
		fmt.Printf("Med   = %.1f kN*m\n", Med*0.001)
		fmt.Printf("Mc_Rd = %.1f kN*m\n", Mc_Rd*0.001)
	}
	return factor
}

func Formula6_13(Wpl float64, Fy float64, gamma_M0 float64, Mpl_Rd float64, out bool) {
	Mpl_Rd = Wpl * Fy / gamma_M0
	if out {
		fmt.Printf("Calculate Mpl_Rd by formula (6.13) EN1993-1-1\n")
		fmt.Printf("Mpl_Rd = W*Fy/gamma_M0\n")
		fmt.Printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wpl*1e+06, Fy*1e-06, gamma_M0)
		fmt.Printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	}
}

func Formula6_14(Wei float64, Fy float64, gamma_M0 float64, Mpl_Rd float64, out bool) {
	Mpl_Rd = Wei * Fy / gamma_M0
	if out {
		fmt.Printf("Calculate Mpl_Rd by formula (6.14) EN1993-1-1\n")
		fmt.Printf("Mpl_Rd = W*Fy/gamma_M0\n")
		fmt.Printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wei*1e+06, Fy*1e-06, gamma_M0)
		fmt.Printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	}
}

func Formula6_15(Wei float64, Fy float64, gamma_M0 float64, Mpl_Rd float64, out bool) {
	Mpl_Rd = Wei * Fy / gamma_M0
	if out {
		fmt.Printf("Calculate Mpl_Rd by formula (6.15) EN1993-1-1\n")
		fmt.Printf("Mpl_Rd = W*Fy/gamma_M0\n")
		fmt.Printf("Mpl_Rd = %.1f*%.1f/%.3f\n", Wei*1e+06, Fy*1e-06, gamma_M0)
		fmt.Printf("Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	}
}

func Table_5_2_tube(d float64, t float64, Fy float64, out bool) CrossSection {
	var CS CrossSection
	if Fy < 1e+08 || t > d {
		fmt.Printf("Error in Table_5_2_tube")
		panic("sss")
	}
	var e float64 = math.Sqrt(2.35e+08 / Fy)
	if d/t <= 50*math.Pow(e, 2) {
		CS = CrossSection_Class1
	} else if d/t <= 70*math.Pow(e, 2) {
		CS = CrossSection_Class2
	} else if d/t <= 90*math.Pow(e, 2) {
		CS = CrossSection_Class3
	} else {
		fmt.Printf("Error in Table_5_2_tube. See EN1993-1-6")
		panic("sss")
	}
	if out {
		fmt.Printf("Tube %.1fx%.1f mm - by table 5.2 EN1993-1-1 class ", d*0, t*0)
		if (CS) == (CrossSection_Class1) {
			fmt.Printf("1\n")
		} else if (CS) == (CrossSection_Class2) {
			fmt.Printf("2\n")
		} else if (CS) == (CrossSection_Class3) {
			fmt.Printf("3\n")
		}
	}
	return CS
}

func Table_5_2_Isection(h float64, b float64, tw float64, tf float64, Fy float64, out bool) CrossSection {
	var CS1 CrossSection
	var CS2 CrossSection
	if Fy < 1e+08 || tf < tw || b < tf || b < tw {
		fmt.Printf("Error in Table_5_2_Isection")
		panic("sss")
	}
	var e float64 = math.Sqrt(2.35e+08 / Fy)
	var c float64
	var t float64
	c = h - 2*tf
	t = tw
	if c/t <= 33*e {
		CS1 = CrossSection_Class1
	} else if c/t <= 38*e {
		CS1 = CrossSection_Class2
	} else if c/t <= 42*e {
		CS1 = CrossSection_Class3
	} else {
		fmt.Printf("Error in Table_5_2_Isection. See EN1993-1-6")
		panic("sss")
	}
	c = b/2 - 2*tw
	t = tf
	if c/t <= 9*e {
		CS2 = CrossSection_Class1
	} else if c/t <= 10*e {
		CS2 = CrossSection_Class2
	} else if c/t <= 14*e {
		CS2 = CrossSection_Class3
	} else {
		fmt.Printf("Error in Table_5_2_Isection. See EN1993-1-6")
		panic("sss")
	}
	var CS CrossSection = CS1
	if (CS) > (CS2) {
		CS = CS2
	}
	if out {
		fmt.Printf("I-section - by table 5.2 EN1993-1-1 class ")
		if (CS) == (CrossSection_Class1) {
			fmt.Printf("1\n")
		} else if (CS) == (CrossSection_Class2) {
			fmt.Printf("2\n")
		} else if (CS) == (CrossSection_Class3) {
			fmt.Printf("3\n")
		}
	}
	return CS
}

func Printf(Bc BUCKLING_CURVE) {
	// Buckling curve
	fmt.Printf("Buckling curve: ")
	switch Bc {
	case (BUCKLING_CURVE_Ao):
		fmt.Printf("ao")
	case (BUCKLING_CURVE_A):
		fmt.Printf("a ")
	case (BUCKLING_CURVE_B):
		fmt.Printf("b ")
	case (BUCKLING_CURVE_C):
		fmt.Printf("c ")
	case (BUCKLING_CURVE_D):
		fmt.Printf("d ")
	default:
		fmt.Printf("ERROR in EN1991_1_1_Table6_1\n")
		panic("sss")
	}
	fmt.Printf("\n")
}

func Table6_1_Alpha(Bc BUCKLING_CURVE, out bool) float64 {
	var alpha float64 = float64(-1)
	switch Bc {
	case (BUCKLING_CURVE_Ao):
		alpha = 0.13
	case (BUCKLING_CURVE_A):
		alpha = 0.21
	case (BUCKLING_CURVE_B):
		alpha = 0.34
	case (BUCKLING_CURVE_C):
		alpha = 0.49
	case (BUCKLING_CURVE_D):
		alpha = 0.76
	default:
		fmt.Printf("ERROR in EN1991_1_1_Table6_1\n")
		panic("sss")
	}
	if out {
		Printf(Bc)
		fmt.Printf("Imperfection factor alpha = %.2f by Table 6.1 EN1993-1-1\n", alpha)
	}
	return alpha
}

func table_6_2_Isection(tf float64, Category string, YY BUCKLING_CURVE, ZZ BUCKLING_CURVE, out bool) {
	if tf > 0.1 || tf < 0 {
		fmt.Printf("Please check - stupid tf in table_6_2_Isection")
		panic("sss")
	}
	if Category == "C235" || Category == "C275" || Category == "C355" || Category == "C420" {
		if tf < 0.04 {
			YY = BUCKLING_CURVE_B
			ZZ = BUCKLING_CURVE_C
		} else {
			YY = BUCKLING_CURVE_C
			ZZ = BUCKLING_CURVE_D
		}
	} else {
		fmt.Printf("Please add this material to table_6_2_Isection")
		panic("sss")
	}
	if out {
		fmt.Printf("I-section by table 6.2:\n")
		fmt.Printf("for axes Y-Y :")
		Printf(YY)
		fmt.Printf("\n")
		fmt.Printf("for axes Z-Z :")
		Printf(ZZ)
		fmt.Printf("\n")
	}
}

func Formula6_47(Xi float64, A float64, Fy float64, gamma_M1 float64, beam CrossSection, Nb_Rd float64, out bool) {
	if (beam) != (CrossSection_Class1) && (beam) != (CrossSection_Class2) && (beam) != (CrossSection_Class3) {
		fmt.Printf("Please use another formula. Formula6_47")
		panic("sss")
	}
	Nb_Rd = Xi * A * Fy / gamma_M1
	if out {
		fmt.Printf("Calculate Nb_Rd by formula (6.47) EN1993-1-1\n")
		fmt.Printf("Nb_Rd = Xi*A*Fy/gamma_M1\n")
		fmt.Printf("Nb_Rd = %.1f*%.1f*%.1f/%.3f\n", Xi, A*1e+06, Fy*1e-06, gamma_M1)
		fmt.Printf("Nb_Rd = %.1f kN\n", Nb_Rd*0.001)
	}
}

func Formula6_49_Xi(Lambda_ float64, Bc BUCKLING_CURVE, out bool) float64 {
	var Fi float64 = 0.5 * (1 + Table6_1_Alpha(Bc, out)*(Lambda_-0.2) + math.Pow(Lambda_, 2))
	var Xi float64 = 1 / (Fi + math.Sqrt(Fi*Fi-Lambda_*Lambda_))
	if Xi > 1 {
		Xi = 1
	}
	if out {
		fmt.Printf("Fi = %f\n", Fi)
		fmt.Printf("Xi = %f\n", Xi)
	}
	return Xi
}

func Par6_3_1_2_Lambda_c1(A float64, Fy float64, Ncr float64, Lambda_ float64, out bool) {
	Lambda_ = math.Sqrt(A * Fy / Ncr)
	if out {
		fmt.Printf("ewrwerwe add good view")
		fmt.Printf("Lambda_ = %.2f\n", Lambda_)
	}
}

func Par6_3_1_2_Lambda_c2(Lcr float64, i float64, Lambda1 float64, out bool) float64 {
	var Lambda_ float64 = Lcr / i / Lambda1
	if out {
		fmt.Printf("Lambda_ = %.2f\n", Lambda_)
	}
	return Lambda_
}

func Par6_50_Lambda1(E float64, Fy float64, out bool) float64 {
	var Lambda1 float64 = math.Pi * math.Sqrt(E/Fy)
	if out {
		fmt.Printf("Lambda1 = PI*sqrt(E/Fy) = PI*sqrt(%.2f/%.2f) = %f\n", E, Fy, Lambda1)
	}
	return Lambda1
}

func Ncr(ELASITY float64, J_moment_inetria float64, L float64, Ncr float64, out bool) {
	Ncr = math.Pow(math.Pi/L, 2) * ELASITY * J_moment_inetria
	if out {
		fmt.Printf("Ncr = (PI/L)^2*E*J\n")
		fmt.Printf("Ncr = (PI/%.1f)^2*%.1f*%.1f\n", L*0, ELASITY*1e-06, J_moment_inetria*1e+08)
		fmt.Printf("Ncr = %.1f kN\n", Ncr*0.001)
	}
}
