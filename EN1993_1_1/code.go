package en1993_1_1

import (
	"fmt"
	"io"
	"math"

	"github.com/Konstantin8105/connection/unit"
)

// constants by par.3.2.6
const (
	E  = 210000.0 // Modulus of elasticity. N/sq.mm
	G  = 81000.0  // Shear modulus. N/sq.mm
	Mν = 0.3      // Poissin`s ratio in elastic stage
	Mα = 12e-6    // Coefficient of linear thermal expansion for T <= 100deg.C
)

type Steel string

const (
	S235EN10025_2 Steel = "S235 EN10025-2"
	S275EN10025_2 Steel = "S275 EN10025-2"
	S355EN10025_2 Steel = "S355 EN10025-2"
	S450EN10025_2 Steel = "S470 EN10025-2"

	S275N_NLEN10025_3 Steel = "S275 N/NL EN10025-3"
	S355N_NLEN10025_3 Steel = "S355 N/NL EN10025-3"
	S420N_NLEN10025_3 Steel = "S420 N/NL EN10025-3"
	S460N_NLEN10025_3 Steel = "S460 N/NL EN10025-3"

	S275M_MLEN10025_4 Steel = "S275 M/ML EN10025-4"
	S355M_MLEN10025_4 Steel = "S355 M/ML EN10025-4"
	S420M_MLEN10025_4 Steel = "S420 M/ML EN10025-4"
	S460M_MLEN10025_4 Steel = "S460 M/ML EN10025-4"

	S235WEN10025_5 Steel = "S235 W EN10025-5"
	S355WEN10025_5 Steel = "S335 W EN10025-5"

	S460Q_QL_QL1WEN10025_6 Steel = "S460 Q/QL/QL1 EN10025-6"

	S235HEN10210_1 Steel = "S235 H EN10210-1"
	S275HEN10210_1 Steel = "S275 H EN10210-1"
	S355HEN10210_1 Steel = "S355 H EN10210-1"
)

func (s Steel) String() string {
	return string(s)
}

// TODO: add unit everythere

// Strength by table 3.1.
// Units:
//	t  - nominal thickness, mm
//	Fy - nominal value of yield strength, Pa
//	Fu - nominal value of ultimate tensile strength, Pa
func (s Steel) Strength(w io.Writer, t unit.Meter) (Fy, Fu unit.Pa) {
	t *= 1000.0        // convert from meter to mm
	var fy, fu float64 // unit: MPa
	for _, data := range []struct {
		steel    Steel
		fy1, fu1 float64 // MPa
		fy2, fu2 float64 // MPa
	}{
		// data of materials
		{S235EN10025_2, 235.0, 360.0, 215.0, 360.0},
		{S275EN10025_2, 275.0, 430.0, 255.0, 410.0},
		{S355EN10025_2, 355.0, 510.0, 335.0, 470.0},
		{S450EN10025_2, 440.0, 550.0, 410.0, 550.0},

		{S275N_NLEN10025_3, 275.0, 390.0, 255.0, 370.0},
		{S355N_NLEN10025_3, 355.0, 490.0, 335.0, 470.0},
		{S420N_NLEN10025_3, 420.0, 520.0, 390.0, 520.0},
		{S460N_NLEN10025_3, 460.0, 540.0, 430.0, 540.0},

		{S275M_MLEN10025_4, 275.0, 370.0, 255.0, 360.0},
		{S355M_MLEN10025_4, 355.0, 470.0, 335.0, 450.0},
		{S420M_MLEN10025_4, 420.0, 520.0, 390.0, 500.0},
		{S460M_MLEN10025_4, 460.0, 540.0, 430.0, 530.0},

		{S235WEN10025_5, 235.0, 360.0, 215.0, 340.0},
		{S355WEN10025_5, 355.0, 510.0, 335.0, 490.0},

		{S460Q_QL_QL1WEN10025_6, 460.0, 570.0, 440.0, 550.0},

		{S235HEN10210_1, 235.0, 360.0, 215.0, 340.0},
		{S275HEN10210_1, 275.0, 430.0, 255.0, 410.0},
		{S355HEN10210_1, 355.0, 510.0, 335.0, 490.0},
	} {
		if data.steel != s {
			continue
		}
		if t <= 40 {
			fy, fu = data.fy1, data.fu1
			break
		}
		if 40 < t && t <= 80 {
			fy, fu = data.fy2, data.fu2
			break
		}
		panic(fmt.Errorf("Not implemented material strength for %s with thk=%.2f mm", s, t))
	}
	// convert from MPa to Pa
	Fy = unit.Pa(fy * 1.0e6)
	Fu = unit.Pa(fu * 1.0e6)

	fmt.Fprintf(w, "For steel %s with thickness %.2 mm:\nFy = %.1f MPa\nFu = %.1f MPa\n",
		s, t.ToMM(), Fy.ToMPa(), Fu.ToMPa)
	return
}

// Part11 - for store factors par.6.1 by Eurocode 3: Design of steel structures.
// Part 1-1: General rules and rules for building
type Part11 struct {
	γm0 float64 // resistance of cross-sections whever the class
	γm1 float64 // resistance of members to instability assessed by member checks
	γm2 float64 // resistance of cross-sections in tension to fracture
}

// New - all factors by note 2B par.6.1
func New() *Part11 {
	return &Part11{
		γm0: 1.00,
		γm1: 1.00,
		γm2: 1.25,
	}
}

// Generate with custom factors
func Generate(γm0, γm1, γm2 float64) *Part11 {
	return &Part11{
		γm0: γm0,
		γm1: γm1,
		γm2: γm2,
	}
}

func (p Part11) Formula6_6(w io.Writer, A, Fy float64) (Npl_Rd float64) {
	Npl_Rd = A * Fy / p.γm0

	fmt.Fprintf(w, "Calculate Npl_Rd by formula (6.6) EN1993-1-1\n")
	fmt.Fprintf(w, "Npl_Rd = A*Fy/γm0\n")
	fmt.Fprintf(w, "Npl_Rd = %.1f*%.1f/%.3f\n", A*1e+06, Fy*1e-06, p.γm0)
	fmt.Fprintf(w, "Npl_Rd = %.1f kN\n", Npl_Rd*0.001)
	return
}

func (p Part11) Formula6_12(w io.Writer, Med, Mc_Rd float64) (ratio float64) {
	ratio = Med / Mc_Rd

	fmt.Fprintf(w, "Calculate by formula (6.12) EN1993-1-1\n")
	fmt.Fprintf(w, "Med   = %.1f kN*m\n", Med*0.001)
	fmt.Fprintf(w, "Mc_Rd = %.1f kN*m\n", Mc_Rd*0.001)
	fmt.Fprintf(w, "Med/Mc_Rd = %.2f <= 1\n", ratio)
	return
}

func (p Part11) Formula6_13(w io.Writer, Wpl, Fy, Mpl_Rd float64) {
	Mpl_Rd = Wpl * Fy / p.γm0

	fmt.Fprintf(w, "Calculate Mpl_Rd by formula (6.13) EN1993-1-1\n")
	fmt.Fprintf(w, "Mpl_Rd = W*Fy/γm0\n")
	fmt.Fprintf(w, "Mpl_Rd = %.1f*%.1f/%.3f\n", Wpl*1e+06, Fy*1e-06, p.γm0)
	fmt.Fprintf(w, "Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	return
}

func (p Part11) Formula6_14(w io.Writer, Wei, Fy float64) (Mpl_Rd float64) {
	Mpl_Rd = Wei * Fy / p.γm0

	fmt.Fprintf(w, "Calculate Mpl_Rd by formula (6.14) EN1993-1-1\n")
	fmt.Fprintf(w, "Mpl_Rd = W*Fy/γm0\n")
	fmt.Fprintf(w, "Mpl_Rd = %.1f*%.1f/%.3f\n", Wei*1e+06, Fy*1e-06, p.γm0)
	fmt.Fprintf(w, "Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	return
}

func (p Part11) Formula6_15(w io.Writer, Wei, Fy float64) (Mpl_Rd float64) {
	Mpl_Rd = Wei * Fy / p.γm0

	fmt.Fprintf(w, "Calculate Mpl_Rd by formula (6.15) EN1993-1-1\n")
	fmt.Fprintf(w, "Mpl_Rd = W*Fy/γm0\n")
	fmt.Fprintf(w, "Mpl_Rd = %.1f*%.1f/%.3f\n", Wei*1e+06, Fy*1e-06, p.γm0)
	fmt.Fprintf(w, "Mpl_Rd = %.1f kN*m\n", Mpl_Rd*0.001)
	return
}

func (p Part11) Formula6_18(w io.Writer, Fy, Av float64) (Vpl_Rd float64) {
	Vpl_Rd = Fy / math.Sqrt(3) * Av / p.γm0

	fmt.Fprintf(w, "Calculation: Shear Yielding Verification of the Plate\n")
	fmt.Fprintf(w, "Shear Yielding Verification of the Plate\n")
	fmt.Fprintf(w, "Ved <= Vpl_Rd, see formula (6.17) EN1993-1-1\n")
	fmt.Fprintf(w, "Vpl_Rd = Fy*Av/γm0/sqrt(3)\n")
	fmt.Fprintf(w, "Av = %.1f mm2\n", Av*1e+06)
	fmt.Fprintf(w, "Vpl_Rd = %.1f * %.1f /%.2f/%.2f\n", Fy*1e-06, Av*1e+06, p.γm0, math.Sqrt(3))
	fmt.Fprintf(w, "Vpl_Rd = %.1f kN\n", Vpl_Rd*0.001)
	return
}

type Class string

const (
	Class1 Class = "Class 1"
	Class2       = "Class 2"
	Class3       = "Class 3"
	Class4       = "Class 4"
)

func (c Class) String() string {
	return string(c)
}

type SectionClassType string

const (
	Internal SectionClassType = "Internal compression parts"
	Outstand                  = "Outstand flanges"
	Angles                    = "Angles"
	Tubes                     = "Tubular sections"
)

// compression section part
func Table5_2(s SectionClassType, c, t float64, Fy float64) (class Class, err error) {
	e := math.Sqrt(235.0e6 / Fy)
	ratio := c / t
	for _, p := range []struct {
		s        SectionClassType
		ratioMax float64
		class    Class
	}{
		{Internal, 33.0 * e, Class1},
		{Internal, 38.0 * e, Class2},
		{Internal, 42.0 * e, Class3},
		{Outstand, 9.0 * e, Class1},
		{Outstand, 10.0 * e, Class2},
		{Outstand, 14.0 * e, Class3},
		{Tubes, 50.0 * e * e, Class1},
		{Tubes, 70.0 * e * e, Class2},
		{Tubes, 90.0 * e * e, Class3},
	} {
		if s == p.s && ratio < p.ratioMax {
			return p.class, nil
		}
	}
	if s == Angles && ratio < 15*e && (c+c)/2/t <= 11.5*e {
		return Class3, nil
	}

	err = fmt.Errorf("not acceptable parameters for %s: c=%.2e, f=%.2e, Fy=%.2e", s, c, t, Fy)
	return
}

// func Table_5_2_Isection(w io.Writer, h, b, tw, tf, Fy float64) Class {
// 	var CS1 CrossSection
// 	var CS2 CrossSection
// 	if Fy < 1e+08 || tf < tw || b < tf || b < tw {
// 		fmt.Fprintf(w, "Error in Table_5_2_Isection")
// 		panic("sss")
// 	}
// 	var e float64 = math.Sqrt(2.35e+08 / Fy)
// 	var c float64
// 	var t float64
// 	c = h - 2*tf
// 	t = tw
// 	if c/t <= 33*e {
// 		CS1 = CrossSection_Class1
// 	} else if c/t <= 38*e {
// 		CS1 = CrossSection_Class2
// 	} else if c/t <= 42*e {
// 		CS1 = CrossSection_Class3
// 	} else {
// 		fmt.Fprintf(w, "Error in Table_5_2_Isection. See EN1993-1-6")
// 		panic("sss")
// 	}
// 	c = b/2 - 2*tw
// 	t = tf
// 	if c/t <= 9*e {
// 		CS2 = CrossSection_Class1
// 	} else if c/t <= 10*e {
// 		CS2 = CrossSection_Class2
// 	} else if c/t <= 14*e {
// 		CS2 = CrossSection_Class3
// 	} else {
// 		fmt.Fprintf(w, "Error in Table_5_2_Isection. See EN1993-1-6")
// 		panic("sss")
// 	}
// 	var CS CrossSection = CS1
// 	if (CS) > (CS2) {
// 		CS = CS2
// 	}
// 	if out {
// 		fmt.Fprintf(w, "I-section - by table 5.2 EN1993-1-1 class ")
// 		if (CS) == (CrossSection_Class1) {
// 			fmt.Fprintf(w, "1\n")
// 		} else if (CS) == (CrossSection_Class2) {
// 			fmt.Fprintf(w, "2\n")
// 		} else if (CS) == (CrossSection_Class3) {
// 			fmt.Fprintf(w, "3\n")
// 		}
// 	}
// 	return CS
// }
//
// func Table_5_2_tube(w io.Writer, d float64, t float64, Fy float64, out bool) CrossSection {
// 	var CS CrossSection
// 	if Fy < 1e+08 || t > d {
// 		fmt.Fprintf(w, "Error in Table_5_2_tube")
// 		panic("sss")
// 	}
// 	var e float64 = math.Sqrt(2.35e+08 / Fy)
// 	if d/t <= 50*math.Pow(e, 2) {
// 		CS = CrossSection_Class1
// 	} else if d/t <= 70*math.Pow(e, 2) {
// 		CS = CrossSection_Class2
// 	} else if d/t <= 90*math.Pow(e, 2) {
// 		CS = CrossSection_Class3
// 	} else {
// 		fmt.Fprintf(w, "Error in Table_5_2_tube. See EN1993-1-6")
// 		panic("sss")
// 	}
// 	if out {
// 		fmt.Fprintf(w, "Tube %.1fx%.1f mm - by table 5.2 EN1993-1-1 class ", d*0, t*0)
// 		if (CS) == (CrossSection_Class1) {
// 			fmt.Fprintf(w, "1\n")
// 		} else if (CS) == (CrossSection_Class2) {
// 			fmt.Fprintf(w, "2\n")
// 		} else if (CS) == (CrossSection_Class3) {
// 			fmt.Fprintf(w, "3\n")
// 		}
// 	}
// 	return CS
// }

type BucklingCurve string

const (
	BcAo BucklingCurve = "Buckling curve Ao"
	BcA                = "Buckling curve A"
	BcB                = "Buckling curve B"
	BcC                = "Buckling curve C"
	BcD                = "Buckling curve D"
)

func (b BucklingCurve) String() string {
	return string(b)
}

func Table6_1(w io.Writer, b BucklingCurve) (α float64) {
	switch b {
	case BcAo:
		α = 0.13
	case BcA:
		α = 0.21
	case BcB:
		α = 0.34
	case BcC:
		α = 0.49
	case BcD:
		α = 0.76
	}

	fmt.Fprintf(w, "%s\n", b)
	fmt.Fprintf(w, "Imperfection factor alpha = %.2f by Table 6.1 EN1993-1-1\n", α)
	return
}

// func table_6_2_Isection(w io.Writer, tf float64, Category string, YY BUCKLING_CURVE, ZZ BUCKLING_CURVE, out bool) {
// 	if tf > 0.1 || tf < 0 {
// 		fmt.Fprintf(w, "Please check - stupid tf in table_6_2_Isection")
// 		panic("sss")
// 	}
// 	if Category == "C235" || Category == "C275" || Category == "C355" || Category == "C420" {
// 		if tf < 0.04 {
// 			YY = BUCKLING_CURVE_B
// 			ZZ = BUCKLING_CURVE_C
// 		} else {
// 			YY = BUCKLING_CURVE_C
// 			ZZ = BUCKLING_CURVE_D
// 		}
// 	} else {
// 		fmt.Fprintf(w, "Please add this material to table_6_2_Isection")
// 		panic("sss")
// 	}
// 	if out {
// 		fmt.Fprintf(w, "I-section by table 6.2:\n")
// 		fmt.Fprintf(w, "for axes Y-Y :")
// 		Fprintf(w, YY)
// 		fmt.Fprintf(w, "\n")
// 		fmt.Fprintf(w, "for axes Z-Z :")
// 		Fprintf(w, ZZ)
// 		fmt.Fprintf(w, "\n")
// 	}
// }

func Formula6_47(w io.Writer, Xi float64, A float64, Fy float64, gamma_M1 float64, beam CrossSection, Nb_Rd float64, out bool) {
	if (beam) != (CrossSection_Class1) && (beam) != (CrossSection_Class2) && (beam) != (CrossSection_Class3) {
		fmt.Fprintf(w, "Please use another formula. Formula6_47")
		panic("sss")
	}
	Nb_Rd = Xi * A * Fy / gamma_M1
	if out {
		fmt.Fprintf(w, "Calculate Nb_Rd by formula (6.47) EN1993-1-1\n")
		fmt.Fprintf(w, "Nb_Rd = Xi*A*Fy/gamma_M1\n")
		fmt.Fprintf(w, "Nb_Rd = %.1f*%.1f*%.1f/%.3f\n", Xi, A*1e+06, Fy*1e-06, gamma_M1)
		fmt.Fprintf(w, "Nb_Rd = %.1f kN\n", Nb_Rd*0.001)
	}
}

func Formula6_49_Xi(w io.Writer, Lambda_ float64, Bc BucklingCurve, out bool) float64 {
	var Fi float64 = 0.5 * (1 + Table6_1(w,Bc)*(Lambda_-0.2) + math.Pow(Lambda_, 2))
	var Xi float64 = 1 / (Fi + math.Sqrt(Fi*Fi-Lambda_*Lambda_))
	if Xi > 1 {
		Xi = 1
	}
	if out {
		fmt.Fprintf(w, "Fi = %f\n", Fi)
		fmt.Fprintf(w, "Xi = %f\n", Xi)
	}
	return Xi
}

func Par6_3_1_2_Lambda_c1(w io.Writer, A float64, Fy float64, Ncr float64, Lambda_ float64, out bool) {
	Lambda_ = math.Sqrt(A * Fy / Ncr)
	if out {
		fmt.Fprintf(w, "ewrwerwe add good view")
		fmt.Fprintf(w, "Lambda_ = %.2f\n", Lambda_)
	}
}

func Par6_3_1_2_Lambda_c2(w io.Writer, Lcr float64, i float64, Lambda1 float64, out bool) float64 {
	var Lambda_ float64 = Lcr / i / Lambda1
	if out {
		fmt.Fprintf(w, "Lambda_ = %.2f\n", Lambda_)
	}
	return Lambda_
}

func Par6_50_Lambda1(w io.Writer, E float64, Fy float64, out bool) float64 {
	var Lambda1 float64 = math.Pi * math.Sqrt(E/Fy)
	if out {
		fmt.Fprintf(w, "Lambda1 = PI*sqrt(E/Fy) = PI*sqrt(%.2f/%.2f) = %f\n", E, Fy, Lambda1)
	}
	return Lambda1
}

func Ncr(w io.Writer, ELASITY float64, J_moment_inetria float64, L float64, Ncr float64, out bool) {
	Ncr = math.Pow(math.Pi/L, 2) * ELASITY * J_moment_inetria
	if out {
		fmt.Fprintf(w, "Ncr = (PI/L)^2*E*J\n")
		fmt.Fprintf(w, "Ncr = (PI/%.1f)^2*%.1f*%.1f\n", L*0, ELASITY*1e-06, J_moment_inetria*1e+08)
		fmt.Fprintf(w, "Ncr = %.1f kN\n", Ncr*0.001)
	}
}
