package eu3

import (
	"fmt"
	"math"
)

var DEBUG int = 1

type CODE_NORM = int

const (
	Code_SNiP     CODE_NORM = 0
	Code_Eurocode           = 1
)

var CODE CODE_NORM

type STANDART = int

const (
	STANDART_SNiP     STANDART = 0
	STANDART_Eurocode          = 1
)

type STEEL = int

const (
	C235 STEEL = 0
)


func Printf_CALC(a float64, b float64, g int) {
}

func Graph(z1 []float64, z2 []float64, a float64, a2 float64) float64 {
	return 0
}

func Strcmp(s1, s2 string) int {
	return 0
}

func Betta_W_table4_1(Category string) float64 {
	if Strcmp(Category, ("C235\x00")) == 0 {
		return 0.8
	} else {
		fmt.Printf(("Error in Betta_W_table4_1\x00"))
		panic("ddd")
	}
	return 1e+30
}

func SNiP(Category string, thk float64, PY float64, PU float64) {
	if uint(int((CODE))) == uint(int((Code_SNiP))) {
		if Strcmp(Category, ("C235\x00")) == 0 {
			if 0.002 <= thk && thk <= 0.02 {
				//(Category == "C235")
				PY = 2.3e+08
				PU = 3.5e+08
			} else if 0.02 <= thk && thk <= 0.04 {
				PY = 2.2e+08
				PU = 3.5e+08
			} else if 0.04 <= thk && thk <= 0.1 {
				PY = 2.1e+08
				PU = 3.5e+08
			} else if 0.1 <= thk {
				PY = 1.9e+08
				PU = 3.5e+08
			} else {
				fmt.Printf(("You can`t use C235\x00"))
				fmt.Printf(("Category = %s\n\x00"), Category)
				fmt.Printf(("thk = %f\n\x00"), thk)
				panic("ddd")
			}
		}
		if Strcmp(Category, ("C245\x00")) == 0 {
			if 0.002 <= thk && thk <= 0.02 {
				//(Category == "C245")
				PY = 2.4e+08
				PU = 3.6e+08
			} else {
				fmt.Printf(("You can`t use C245\x00"))
				panic("ddd")
			}
		}
	} else if uint(int((CODE))) == uint(int((Code_Eurocode))) {
		if Strcmp(Category, ("C235\x00")) == 0 {
			if thk <= 0.04 {
				//(Category == "C235")
				PY = 2.35e+08
				PU = 3.6e+08
			} else {
				PY = 2.15e+08
				PU = 3.4e+08
			}
		} else if Strcmp(Category, ("C275\x00")) == 0 {
			if thk <= 0.04 {
				//(Category == "C275")
				PY = 2.75e+08
				PU = 4.3e+08
			} else {
				PY = 2.55e+08
				PU = 4.1e+08
			}
		} else {
			fmt.Printf(("You can`t use this with Eurocode\x00"))
			panic("ddd")
		}
	} else {
		fmt.Printf(("You can`t use this: What is this code\x00"))
		panic("ddd")
	}
}

func STEEL_STRESS(Steels STEEL, thk float64, Py float64, Pu float64, STD STANDART, OUT int) {
	if uint(int((STD))) == uint(int((STANDART_SNiP))) {
		switch uint(int((Steels))) {
		case uint(int((C235))):
			if 0.002 <= thk && thk <= 0.02 {
				Py = 2.3e+08
				Pu = 3.5e+08
			} else if 0.02 <= thk && thk <= 0.04 {
				Py = 2.2e+08
				Pu = 3.5e+08
			} else if 0.04 <= thk && thk <= 0.1 {
				Py = 2.1e+08
				Pu = 3.5e+08
			} else if 0.1 <= thk {
				Py = 1.9e+08
				Pu = 3.5e+08
			} else {
				fmt.Printf(("You can`t use C235\x00"))
				panic("ddd")
			}
		default:
			fmt.Printf(("ERROR: enum STEEL_STRESS\x00"))
			panic("ddd")
		}
	} else if uint(int((STD))) == uint(int((STANDART_SNiP))) {
		switch uint(int((Steels))) {
		case uint(int((C235))):
			if thk <= 0.04 {
				Py = 2.35e+08
				Pu = 3.6e+08
			} else {
				Py = 2.15e+08
				Pu = 3.4e+08
			}
			fallthrough
		default:
			fmt.Printf(("ERROR: enum STEEL_STRESS\x00"))
			panic("ddd")
		}
	} else {
		panic("ddd")
	}
	if OUT != 0 {
		fmt.Printf(("Py = %5.1 MPa\n\x00"), Py*1e-06)
		fmt.Printf(("Pu = %5.1 MPa\n\x00"), Pu*1e-06)
	}
}

// Picture_6_1_APLHA - transpiled function from  GOPATH/src/github.com/Konstantin8105/connection/eu3/eu8.c:322
func Picture_6_1_APLHA(lambda1 float64, lambda2 float64, OUT int) float64 {
	if DEBUG != 0 {
		///
		/// BOLT CONNECTION: COEFFICIENT FROM PICTURE 6.1
		///
		fmt.Printf(("lambda1=%.2f,lambda2=%.2f\n\x00"), lambda1, lambda2)
	}
	if lambda1 < 0 || lambda2 < 0 {
		fmt.Printf("Stupid check\n")
		fmt.Printf(("lambda1=%+.2f,lambda2=%+.2f\n\x00"), lambda1, lambda2)
	}
	var ALPHA float64 = 1e+30
	var TABLE [][]float64 = [][]float64{{8, 0.9625419, -3.6454601, 5.0308278, 3.0586037}, {7, 0.9636739, -2.6572564, 3.9919931, 2.9602062}, {2 * math.Pi, 0.9726341, -1.9932994, 3.2341635, 2.7782066}, {6, 0.9496498, -1.480309, 2.6267325, 3.1078904}, {5.5, 0.9785034, -1.1933579, 2.3203117, 2.9238812}, {5, 0.9637229, -0.7743167, 1.9115241, 2.9607118}, {4.75, 0.9946669, -0.7516043, 1.9434759, 1.9056765}, {4.5, 1.2055126, -3.8888717, 7.1149599, 0.9088292}, {4.45, 1.2265766, -3.2360189, 6.2854489, 1.0927973}}
	var LAMBDA1 []float64 = make([]float64, 9)
	{
		var i int
		for ; i < 9; i++ {
			var K float64 = TABLE[i][1]
			var A float64 = TABLE[i][2]
			var B float64 = TABLE[i][3]
			var C float64 = TABLE[i][4]
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
			//    if(lambda1 <= F1) value1 = 2* math.Pi;
			////    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
			//    if(lambda1 >= F2) value1 =  4.45 ;
			////    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
			//    if(F1 < lambda1 && lambda1 < F2)
			//    {
			//        if(lambda2 >= 0.45)
			//            value1 = math.Min(F3,2*math.Pi);
			////        printf("ALPHA3   = %.4e\n",ALPHA  );//DEBUG
			//        if((0.2768*lambda1+0.14) <= lambda2 && lambda2 < 0.45)
			//            value1 = math.Min(F4,2*math.Pi);
			////        printf("ALPHA4   = %.4e\n",ALPHA  );//DEBUG
			//        if((1.2971*lambda1-0.7782) <= lambda2 && lambda2 < (0.2768*lambda1+0.14))
			//            value1 = math.Min(F5,2*math.Pi);
			////        printf("ALPHA5   = %.4e\n",ALPHA  );//DEBUG
			//        if(lambda2 < (1.2971*lambda1-0.7782))
			//            value1 = math.Min(F6,2*math.Pi);
			////        printf("ALPHA6   = %.4e\n",ALPHA  );//DEBUG
			//    }
			//    printf("value1 = %f\n",value1);
			/////////////////////////////////
			/// TABLE
			/////////////////////////////////
			//    printf("DEBUG 05\n");//DEBUG
			LAMBDA1[i] = K + A*lambda2/math.Pow(1+math.Pow(B*lambda2, C), 1/C)
			if DEBUG != 0 {
				fmt.Printf(("[%+.5e, %+.5e, %+.5e, %+.5e] = %+.5e\n\x00"), K, A, B, C, LAMBDA1[i])
			}
		}
	}
	{
		var x []float64 = make([]float64, 9)
		var y []float64 = make([]float64, 9)
		{
			var i int
			for ; i < 9; i++ {
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
				//Graph(Array<double> x, Array<double> y, double pX, double pY)
				//     Array<double> x,y;
				x[i] = TABLE[0][5*i]
				y[i] = LAMBDA1[i]
			}
		}
		Graph(y, x, lambda1, ALPHA)
	}
	if OUT != 0 || DEBUG != 0 {
		//         x.Delete();
		//         y.Delete();
		//    printf("F1      = %.4e\n",F1);     //DEBUG
		//    printf("F2      = %.4e\n",F2);     //DEBUG
		//    printf("F3      = %.4e\n",F3);     //DEBUG
		//    printf("lambda1 = %.4e\n",lambda1);//DEBUG
		//    printf("lambda2 = %.4e\n",lambda2);//DEBUG
		//    printf("ALPHA   = %.4e\n",ALPHA  );//DEBUG
		fmt.Printf(("Alpha(lambda1=%.2f,lambda2=%.2f) = %.2f by Figure 6.11 EN1993-1-8\n\x00"), lambda1, lambda2, ALPHA)
	}
	return ALPHA
}

func Formula3_2(Fu float64, D float64, Thk float64, gammaM2 float64, Fb_Rd float64, OUT int) {
	///
	/// Design bearing resistance
	///
	Fb_Rd = 1.5 * Fu * D * Thk / gammaM2
	if OUT != 0 {
		fmt.Printf("Design bearing resistance by Formula(3.2) EN1993-1-8\n")
		fmt.Printf("Fb_Rd = 1.5*Fu*d*t/gammaM2\n")
		fmt.Printf(("Fb_Rd = 1.5*%.1f*%.1f*%.1f/%.2f\n\x00"), Fu*1e-06, D*1000, Thk*1000, gammaM2)
		fmt.Printf(("Fb_Rd = %.1f kN\n\x00"), Fb_Rd*0.001)
	}
}

func Formula3_9(Fu float64, Ant float64, gammaM2 float64, Fy float64, Anv float64, gammaM0 float64, Veff_t_Rd float64, OUT int) {
	///
	/// Block Tearing Verification
	///
	Veff_t_Rd = Fu*Ant/gammaM2 + 1/math.Sqrt(3)*Fy*Anv/gammaM0
	if OUT != 0 {
		fmt.Printf("Block Tearing Verification\n")
		fmt.Printf(("Ant = %.1f mm2\n\x00"), Ant*1e+06)
		fmt.Printf(("Anv = %.1f mm2\n\x00"), Anv*1e+06)
		fmt.Printf("Veff = Fu*Ant/gamma_M2 + (1/sqrt(3))*Fy*Anv/gamma_M0\n")
		fmt.Printf(("Veff = %.1f *%.1f /%.2f + 1/%.2f *%.1f *%.1f /%.2f \n\x00"), Fu*1e-06, Ant*1e+06, gammaM2, math.Sqrt(3), Fy*1e-06, Anv*1e+06, gammaM0)
		fmt.Printf(("Veff_Rd = %.1f kN\n\x00"), Veff_t_Rd*0.001)
	}
}

func Formula4_1(Sigma_perp float64, Tau_perp float64, Tau_parr float64, Fu float64, Category string, gamma_M2 float64, OUT int) float64 {
	var factor float64
	var SS float64 = math.Sqrt(math.Pow(Sigma_perp, 2) + 3*(math.Pow(Tau_perp, 2)+math.Pow(Tau_parr, 2)))
	var SS_left float64 = Fu / Betta_W_table4_1(Category) / gamma_M2
	var SS_left2 float64 = 0.9 * Fu / gamma_M2
	if Fu < 1e+06 || gamma_M2 < 0.99 {
		///
		/// Weld Verification
		///
		fmt.Printf(("STUPID in Formula4_1\x00"))
		panic("ddd")
	}
	if OUT != 0 {
		fmt.Printf("Weld Verification\n")
		fmt.Printf(("Sigma_perp = %.2f MPa\n\x00"), Sigma_perp*1e-06)
		fmt.Printf(("Tau_perp   = %.2f MPa\n\x00"), Tau_perp*1e-06)
		fmt.Printf(("Tau_parr   = %.2f MPa\n\x00"), Tau_parr*1e-06)
		fmt.Printf("sqrt(Sigma_perp^2+3*(Tau_perp^2+Tau_parr^2)) <= Fu/Betta_W_table4_1/gamma_M2\n")
		fmt.Printf(("sqrt(%.2f^2+3*(%.2f^2+%.2f^2)) <= %.2f/%.2f/%.2f\n\x00"), Sigma_perp*1e-06, Tau_perp*1e-06, Tau_parr*1e-06, Fu*1e-06, Betta_W_table4_1(Category), gamma_M2)
		fmt.Printf(("%.2f MPa<= %.2f MPa\n\x00"), SS*1e-06, SS_left*1e-06)
	}
	Printf_CALC(SS/SS_left, factor, OUT)
	if OUT != 0 {
		fmt.Printf("Sigma_perp <= 0.9*Fu/gamma_M2\n")
		fmt.Printf(("%.2f <= 0.9*%.2f/%.2f\n\x00"), Sigma_perp*1e-06, Fu*1e-06, gamma_M2)
		fmt.Printf(("%.2f MPa <= %.2f MPa\n\x00"), Sigma_perp*1e-06, SS_left2*1e-06)
	}
	Printf_CALC(Sigma_perp/SS_left2, factor, OUT)
	return factor
}

func Table6_2_M_pl_1_RD(l_eff_1 float64, tf float64, Fy float64, gamma_M0 float64, OUT int) float64 {
	var M_pl_1_RD float64 = 0.25 * l_eff_1 * math.Pow(tf, 2) * Fy / gamma_M0
	if OUT != 0 {
		///
		/// TABLE 6.2 EN1993-1-8
		///
		fmt.Printf("M_pl_1_Rd = 0.25*l_eff_1*tf^2*Fy/gamma_M0\n")
		fmt.Printf(("M_pl_1_Rd = 0.25*%.1f*%.1f^2*%.1f/%.2f\n\x00"), l_eff_1*1000, tf*1000, Fy*1e-06, gamma_M0)
		fmt.Printf(("M_pl_1_Rd = %.1f kN*m\n\x00"), M_pl_1_RD*0.001)
	}
	return M_pl_1_RD
}

func Table6_2_M_pl_2_RD(l_eff_2 float64, tf float64, Fy float64, gamma_M0 float64, OUT int) float64 {
	var M_pl_2_RD float64 = 0.25 * l_eff_2 * math.Pow(tf, 2) * Fy / gamma_M0
	if OUT != 0 {
		fmt.Printf("M_pl_2_RD = 0.25*l_eff_2*tf^2*Fy/gamma_M0\n")
		fmt.Printf(("M_pl_2_RD = 0.25*%.1f*%.1f^2*%.1f/%.2f\n\x00"), l_eff_2*1000, tf*1000, Fy*1e-06, gamma_M0)
		fmt.Printf(("M_pl_2_RD = %.1f kN*m\n\x00"), M_pl_2_RD*0.001)
	}
	return M_pl_2_RD
}

func Table6_2_M_bp_RD(l_eff_1 float64, tbp float64, Fy_bp float64, gamma_M0 float64, OUT int) float64 {
	var M_bp_RD float64 = 0.25 * l_eff_1 * math.Pow(tbp, 2) * Fy_bp / gamma_M0
	if OUT != 0 {
		fmt.Printf("M_bp_RD = 0.25*l_eff_1*tbp^2*Fy_bp/gamma_M0\n")
		fmt.Printf(("M_bp_RD = 0.25*%.1f*%.1f^2*%.1f/%.2f\n\x00"), l_eff_1*1000, tbp*1000, Fy_bp*1e-06, gamma_M0)
		fmt.Printf(("M_bp_RD = %.1f kN*m\n\x00"), M_bp_RD*0.001)
	}
	return M_bp_RD
}

func Table6_2_F_t1_Rd(M_pl_1_RD float64, m float64, OUT int) float64 {
	var F_t1_Rd float64 = 4 * M_pl_1_RD / m
	if OUT != 0 {
		fmt.Printf("F_t1_Rd = 4*M_pl_1_RD/m\n")
		fmt.Printf(("F_t1_Rd = 4*%.1f/%.1f\n\x00"), M_pl_1_RD*0.001, m*1000)
		fmt.Printf(("F_t1_Rd = %.1f kN\n\x00"), F_t1_Rd*0.001)
	}
	return F_t1_Rd
}

func Table6_2_F_t2_Rd(M_pl_2_RD float64, n float64, F_t_Rd float64, m float64, number_bolts int, OUT int) float64 {
	var F_t2_Rd float64 = (2*M_pl_2_RD + n*float64(number_bolts)*F_t_Rd) / (m + n)
	if OUT != 0 {
		fmt.Printf("F_t2_Rd = (2*M_pl_2_RD+n*number_bolts*F_t_Rd)/(m+n)\n")
		fmt.Printf(("F_t2_Rd = (2*%.1f+%.1f*%d*%.1f)/(%.1f+%.1f)\n\x00"), M_pl_2_RD*0.001, n*1000, number_bolts, F_t_Rd*0.001, m*1000, n*1000)
		fmt.Printf(("F_t2_Rd = %.1f kN\n\x00"), F_t2_Rd*0.001)
	}
	return F_t2_Rd
}

func Table6_2_F_t3_Rd(F_t_Rd float64, number_bolts int, OUT int) float64 {
	var F_t3_Rd float64 = float64(number_bolts) * F_t_Rd
	if OUT != 0 {
		fmt.Printf("F_t3_Rd = number_bolts*F_t_Rd\n")
		fmt.Printf(("F_t3_Rd = %d*%.1f\n\x00"), number_bolts, F_t_Rd*0.001)
		fmt.Printf(("F_t3_Rd = %.1f kN\n\x00"), F_t3_Rd*0.001)
	}
	return F_t3_Rd
}

func Table6_2_F_t_Rd(F_t1_Rd float64, F_t2_Rd float64, F_t3_Rd float64, OUT int) float64 {
	var F_t_Rd float64 = math.Min(F_t1_Rd, math.Min(F_t2_Rd, F_t3_Rd))
	if OUT != 0 {
		fmt.Printf("F_t_Rd  = math.Min(F_t1_Rd,F_t2_Rd,F_t3_Rd)\n")
		fmt.Printf(("F_t_Rd  = math.Min(%.1f,%.1f,%.1f)\n\x00"), F_t1_Rd*0.001, F_t2_Rd*0.001, F_t3_Rd*0.001)
		fmt.Printf(("F_t_Rd  = %.1f kN\n\x00"), F_t_Rd*0.001)
	}
	return F_t_Rd
}

func Formula6_22_F_t_wb_Rd(Beff_t_wb float64, twb float64, Fy_wb float64, gamma_M0 float64, OUT int) float64 {
	var F_t_wb_Rd float64 = Beff_t_wb * twb * Fy_wb / gamma_M0
	if OUT != 0 {
		fmt.Printf("F_t_wb_Rd = Beff_t_wb * twb * Fy_wb/gamma_M0, see formula(6.22) EN1993-1-8\n")
		fmt.Printf(("F_t_wb_Rd = %.1f * %.1f * %.1f/%.2f\n\x00"), Beff_t_wb*1000, twb*1000, Fy_wb*1e-06, gamma_M0)
		fmt.Printf(("F_t_wb_Rd = %.1f kN\n\x00"), F_t_wb_Rd*0.001)
	}
	return F_t_wb_Rd
}

func Table_6_6_Leff_group(Leff1 float64, Leff2 float64, Paramaters []float64, Row_table int, OUT int) {
	var m float64 = Paramaters[0]
	var p float64 = Paramaters[1]
	var _alpha float64 = Paramaters[2]
	var e float64 = Paramaters[3]
	{
		var i int
		for ; i < 4; i++ {
			if Paramaters[i] < 0 {
				///
				/// Table 6.6
				///
				fmt.Printf("Stupid check in Table_6_6_Leff_group\n")
				{
					var j int
					for ; j < 4; j++ {
						fmt.Printf(("Paramaters[%d] = %.1f mm\x00"), j, Paramaters[j]*1000)
					}
				}
				panic("ddd")
			}
		}
	}
	var Leff_cp float64
	var Leff_nc float64
	switch Row_table {
	case (1):
		fmt.Printf("Bolt-row location: Bolt-row outside tension flange of beam\n")
		panic("ddd")
	case (2):
		Leff_cp = math.Pi*m + p
		Leff_nc = 0.5*p + _alpha*m - (2*m + 0.625*e)
		if OUT != 0 {
			fmt.Printf("Bolt-row location: First bolt-row below tension flange of beam\n")
			fmt.Printf(" Leff_cp = PI*m+p\n")
			fmt.Printf((" Leff_cp = PI*%.1f+%.1f\n\x00"), m*1000, p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e)\n")
			fmt.Printf((" Leff_nc = 0.5*%.1f+%.2f*%.1f-(2*%.1f+0.625*%.1f)\n\x00"), p*1000, _alpha, m*1000, m*1000, e*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	case (3):
		Leff_cp = 2 * p
		Leff_nc = p
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other inner bolt-row\n")
			fmt.Printf(" Leff_cp = 2*p\n")
			fmt.Printf((" Leff_cp = 2*%.1f\n\x00"), p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = p\n")
			fmt.Printf((" Leff_nc = %.1f\n\x00"), p*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	case (4):
		Leff_cp = math.Pi*m + p
		Leff_nc = 2*m + 0.625*e + 0.5*p
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other end bolt-row\n")
			fmt.Printf(" Leff_cp = PI*m+p\n")
			fmt.Printf((" Leff_cp = PI*%.1f+%.1f\n\x00"), m*1000, p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = 2*m+0.625*e+0.5*p\n")
			fmt.Printf((" Leff_nc = 2*%.1f+0.625*%.1f+0.5*%.1f\n\x00"), m*1000, e*1000, p*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	default:
		fmt.Printf(("Error in Leff_table6_6_group\x00"))
		panic("ddd")
	}
	Leff1 = math.Min(Leff_nc, Leff_cp)
	if OUT != 0 {
		fmt.Printf(("  Leff1 = math.Min(Leff_nc,Leff_cp) = math.Min(%.1f,%.1f)\n\x00"), Leff_nc*1000, Leff_cp*1000)
	}
	Leff2 = Leff_nc
	if OUT != 0 {
		fmt.Printf(("  Leff1 = %.1f mm\n\x00"), Leff1*1000)
	}
	if OUT != 0 {
		fmt.Printf(("  Leff2 = %.1f mm\n\x00"), Leff2*1000)
	}
}

func Table_6_6_Leff_individ(Leff1 float64, Leff2 float64, Paramaters []float64, Row_table int, OUT int) {
	var DEBUG int
	var mx float64 = Paramaters[0]
	var w float64 = Paramaters[1]
	var e float64 = Paramaters[2]
	var ex float64 = Paramaters[3]
	var bp float64 = Paramaters[4]
	var _alpha float64 = Paramaters[5]
	var m float64 = Paramaters[6]
	{
		var i int
		for ; i < 7; i++ {
			if Paramaters[i] < 0 {
				fmt.Printf("Stupid check in Table_6_6_Leff_individ\n")
				{
					var j int
					for ; j < 7; j++ {
						fmt.Printf(("Paramaters[%d] = %.1f mm\n\x00"), j, Paramaters[j]*1000)
					}
				}
				panic("ddd")
			}
		}
	}
	if DEBUG != 0 {
		{
			var i int
			for ; i < 7; i++ {
				fmt.Printf(("Paramaters[%d] = %.1f mm\n\x00"), i, Paramaters[i]*1000)
			}
		}
	}
	var Leff_cp float64
	var Leff_nc float64
	switch Row_table {
	case (1):
		var Leff_cp1 float64 = 2 * math.Pi * mx
		var Leff_cp2 float64 = math.Pi*mx + w
		var Leff_cp3 float64 = math.Pi*mx + 2*e
		Leff_cp = math.Min(Leff_cp1, math.Min(Leff_cp2, Leff_cp3))
		var Leff_nc1 float64 = 4*mx + 1.25*ex
		var Leff_nc2 float64 = e + 2*mx + 0.625*ex
		var Leff_nc3 float64 = 0.5 * bp
		var Leff_nc4 float64 = 0.5*w + 2*mx + 0.625*ex
		Leff_nc = math.Min(Leff_nc1, math.Min(Leff_nc2, math.Min(Leff_nc3, Leff_nc4)))
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Bolt-row outside tension flange of beam\n")
			fmt.Printf((" Leff_cp1 = 2*PI*mx = 2*PI* %.1f = %.1f mm\n\x00"), mx*1000, Leff_cp1*1000)
			fmt.Printf((" Leff_cp2 = PI*mx+w = PI* %.1f + %.1f = %.1f mm\n\x00"), mx*1000, w*1000, Leff_cp2*1000)
			fmt.Printf((" Leff_cp3 = PI*mx+2*e = PI* %.1f + 2* %.1f = %.1f mm\n\x00"), mx*1000, e*1000, Leff_cp3*1000)
			fmt.Printf((" Leff_cp = math.Min(Leff_cp1,Leff_cp2,Leff_cp3) = math.Min(%.1f,%.1f,%.1f) = %.1f mm\x00"), Leff_cp1*1000, Leff_cp2*1000, Leff_cp3*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc1 = 4*mx+1.25*ex = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), mx*1000, ex*1000, Leff_nc1*1000)
			fmt.Printf((" Leff_nc2 = e+2*mx+0.625*ex = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n\x00"), e*1000, mx*1000, ex*1000, Leff_nc2*1000)
			fmt.Printf((" Leff_nc3 = 0.5*bp = 0.5 * %.1f = %.1f mm\n\x00"), bp*1000, Leff_nc3*1000)
			fmt.Printf((" Leff_nc4 = 0.5*w+2*mx+0.625*ex = 0.5* %.1f +2* %.1f +0.625* %.1f = %.1f mm\n\x00"), w*1000, mx*1000, ex*1000, Leff_nc4*1000)
			fmt.Printf(" Leff_nc = math.Min(Leff_nc1,Leff_nc2,Leff_nc3,Leff_nc4)\n")
			fmt.Printf((" Leff_nc = math.Min(%.1f,%.1f,%.1f,%.1f) = %.1f mm\n\x00"), Leff_nc1*1000, Leff_nc2*1000, Leff_nc3*1000, Leff_nc4*1000, Leff_nc*1000)
		}
	case (2):
		Leff_cp = 2 * math.Pi * m
		Leff_nc = _alpha * m
		if OUT != 0 {
			fmt.Printf("Bolt-row location: First bolt-row below tension flange of beam\n")
			fmt.Printf((" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = alpha * m = %.2f * %.1f = %.1f mm\n\x00"), _alpha, m*1000, Leff_nc*1000)
		}
	case (3):
		Leff_cp = 2 * math.Pi * m
		Leff_nc = 4*m + 1.25*e
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other inner bolt-row\n")
			fmt.Printf((" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), m*1000, e*1000, Leff_nc*1000)
		}
	case (4):
		Leff_cp = 2 * math.Pi * m
		Leff_nc = 4*m + 1.25*e
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other end bolt-row\n")
			fmt.Printf((" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), m*1000, e*1000, Leff_nc*1000)
		}
	case (5):
		var Leff_cp1 float64 = 2 * math.Pi * mx
		var Leff_cp2 float64 = 2*m + 0.625*e + ex
		var Leff_cp3 float64 = 2 * math.Pi * m
		Leff_cp = math.Min(Leff_cp1, math.Min(Leff_cp2, Leff_cp3))
		var Leff_nc1 float64 = _alpha*mx - (2*mx + 0.625*e) + ex
		var Leff_nc2 float64 = 4*m + 1.25*e
		var Leff_nc3 float64 = e + 2*mx + 0.625*ex
		var Leff_nc4 float64 = 4*mx + 1.25*ex
		Leff_nc = math.Min(Leff_nc1, math.Min(Leff_nc2, math.Min(Leff_nc3, Leff_nc4)))
		if OUT != 0 {
			fmt.Printf("Bolt-row location: bolt-row with gusset between bolts\n")
			fmt.Printf((" Leff_cp1 = 2*PI*mx = 2*PI* %.1f = %.1f mm\n\x00"), mx*1000, Leff_cp1*1000)
			fmt.Printf((" Leff_cp2 = 2*m+0.625*e+ex = 2*%.1f+0.625*%.1f+%.1f = %.1f mm\n\x00"), m*1000, e*1000, ex*1000, Leff_cp2*1000)
			fmt.Printf((" Leff_cp3 = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp3*1000)
			fmt.Printf((" Leff_nc1 = alpha*mx-(2*mx+0.625*e)+ex = %.1f*%.1f-(2*%.1f+0.625*%.1f)+%.1f = %.1f mm\n\x00"), _alpha*1000, mx*1000, mx*1000, e*1000, ex*1000, Leff_nc1*1000)
			fmt.Printf((" Leff_nc2 = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), m*1000, e*1000, Leff_nc2*1000)
			fmt.Printf((" Leff_nc3 = e+2*mx+0.625*ex = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n\x00"), e*1000, mx*1000, ex*1000, Leff_nc3*1000)
			fmt.Printf((" Leff_nc4 = 4*mx+1.25*ex = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), mx*1000, ex*1000, Leff_nc4*1000)
		}
	default:
		fmt.Printf(("Error in Leff_table6_6_individ\x00"))
		panic("ddd")
	}
	Leff1 = math.Min(Leff_nc, Leff_cp)
	if OUT != 0 {
		fmt.Printf(("Leff1 = math.Min(Leff_nc,Leff_cp) = math.Min(%.1f,%.1f)\n\x00"), Leff_nc*1000, Leff_cp*1000)
	}
	Leff2 = Leff_nc
	if OUT != 0 {
		fmt.Printf(("Leff1 = %.1f mm\n\x00"), Leff1*1000)
	}
	if OUT != 0 {
		fmt.Printf(("Leff2 = %.1f mm\n\x00"), Leff2*1000)
	}
}

func Table_3_3_e1(Do float64, e1 float64, THK float64, OUT int) int {
	var ex int
	if 1.2*Do <= e1 && e1 <= 4*THK+0.04 {
		ex = 1
	}
	if OUT != 0 {
		fmt.Printf("Check end distanse e1 (on load direction) by table 3.3 EN1993-1-8:\n")
		fmt.Printf("1.2 * Do <= e1\n")
		fmt.Printf(("1.2 * %.1f <= %.1f\n\x00"), Do*1000, e1*1000)
		fmt.Printf(("%.1f <= %.1f\n\x00"), 1.2*Do*1000, e1*1000)
		fmt.Printf("e1 <= 4*THK+0.040\n")
		fmt.Printf(("%.1f <= 4*%.1f + 0.040\n\x00"), e1*1000, THK*1000)
		fmt.Printf(("%.1f <= %.1f \x00"), e1*1000, (4*THK+0.04)*1000)
		if ex != 0 {
			fmt.Printf("- OK\n")
		} else {
			fmt.Printf("- None\n")
			panic("ddd")
		}
	}
	return ex
}

func Table_3_3_e2(Do float64, e2 float64, THK float64, OUT int) int {
	var ex int
	if 1.2*Do <= e2 {
		//&& e2 <= 4*THK+0.040)
		ex = 1
	}
	if OUT != 0 {
		fmt.Printf("Check end distanse e2 (perpendicular on load direction) by table 3.3 EN1993-1-8:\n")
		fmt.Printf("1.2 * Do <= e2\n")
		fmt.Printf(("1.2 * %.1f <= %.1f\n\x00"), Do*1000, e2*1000)
		fmt.Printf(("%.1f <= %.1f \n\x00"), 1.2*Do*1000, e2*1000)
		fmt.Printf("e2 <= 4*THK+0.040\n")
		fmt.Printf(("%.1f <= 4*%.1f+0.040\n\x00"), e2*1000, THK*1000)
		fmt.Printf(("%.1f <= %.1f\n\x00"), e2*1000, (4*THK+0.04)*1000)
		if ex != 0 {
			fmt.Printf("- OK\n")
		} else {
			//printf("- None\n");
			fmt.Printf("1.2*Do <= e2 - ")
			if 1.2*Do <= e2 {
				fmt.Printf("OK\n")
			} else {
				fmt.Printf("None\n")
			}
			fmt.Printf("e2 <= 4*THK+0.040 - ")
			if e2 <= 4*THK+0.04 {
				fmt.Printf("OK\n")
			} else {
				fmt.Printf("None\n")
			}
			panic("ddd")
		}
	}
	return ex
}

func Table_3_3_p1(Do float64, p1 float64, THK float64, OUT int) int {
	var ex int
	if 2.2*Do <= p1 && p1 <= math.Min(14*THK, 0.2) {
		ex = 1
	}
	if OUT != 0 {
		fmt.Printf("Check spacing p1 by table 3.3 EN1993-1-8:\n")
		fmt.Printf("2.2 * Do <= p1\n")
		fmt.Printf(("2.2 * %.1f <= %.1f\n\x00"), Do*1000, p1*1000)
		fmt.Printf(("%.1f <= %.1f \x00"), 2.2*Do*1000, p1*1000)
		if ex != 0 {
			fmt.Printf("- OK\n")
		} else {
			fmt.Printf("- None\n")
			panic("ddd")
		}
	}
	return ex
}

func Par_4_5_2(a float64, OUT int) int {
	var ex int
	if a >= 0.003 {
		ex = 1
	}
	if OUT != 0 {
		fmt.Printf("Check math.Minimum weld thickness by p.4.5.2(2) EN1993-1-8: ")
		if ex != 0 {
			fmt.Printf("OK\n")
		} else {
			fmt.Printf("None\n")
			panic("ddd")
		}
	}
	return ex
}

func Par_4_5_1(a float64, Leff float64, OUT int) int {
	var ex int
	if Leff >= math.Min(6*a, 0.03) {
		ex = 1
	}
	if OUT != 0 {
		fmt.Printf("Check lenght of weld  by p.4.5.1 EN1993-1-8: ")
		if ex != 0 {
			fmt.Printf("OK\n")
		} else {
			fmt.Printf("None\n")
			panic("ddd")
		}
	}
	return ex
}

func Formula3_11(e2 float64, Do float64, t float64, Fu float64, gamma_M2 float64, Nu_Rd float64, OUT int) {
	Nu_Rd = 2 * (e2 - 0.5*Do) * t * Fu / gamma_M2
	if OUT != 0 {
		fmt.Printf("Calculate by Formula 3.11 EN1993-1-8\n")
		fmt.Printf("Tension force by a single row of one bolt in one leg\n")
		fmt.Printf("Nu_Rd = 2.0*(e2-0.5*Do)*t*Fu/gamma_M2\n")
		fmt.Printf(("Nu_Rd = 2.0*(%.1f-0.5*%.1f)*%.1f*%.1f/%.2f\n\x00"), e2*1000, Do*1000, t*1000, Fu*1e-06, gamma_M2)
		fmt.Printf(("Nu_Rd = %.1f kN\n\x00"), Nu_Rd*0.001)
	}
}

func Formula3_12(p1 float64, Do float64, Anet float64, Fu float64, gamma_M2 float64, Nu_Rd float64, OUT int) {
	var x []float64 = []float64{2.5 * Do, 5 * Do}
	var y []float64 = []float64{0.4, 0.7}
	var betta2 float64
	// Array<double> x,y;
	//     x.AddList(2,
	//     y.AddList(2,);
	Graph(x, y, p1, betta2)
	Nu_Rd = betta2 * Anet * Fu / gamma_M2
	if OUT != 0 {
		fmt.Printf("Calculate by Formula 3.12 EN1993-1-8\n")
		fmt.Printf("Tension force by a single row of two bolts in one leg\n")
		fmt.Printf("Nu_Rd = betta2*Anet*Fu/gamma_M2\n")
		fmt.Printf(("Nu_Rd = %.3f*%.1f*%.1f/%.2f\n\x00"), betta2, Anet*1e+06, Fu*1e-06, gamma_M2)
		fmt.Printf(("Nu_Rd = %.1f kN\n\x00"), Nu_Rd*0.001)
	}
}

func Formula3_13(p1 float64, Do float64, Anet float64, Fu float64, gamma_M2 float64, Nu_Rd float64, OUT int) {
	var x []float64 = []float64{2.5 * Do, 5 * Do}
	var y []float64 = []float64{0.5, 0.7}
	var betta3 float64
	Graph(x, y, p1, betta3)
	Nu_Rd = betta3 * Anet * Fu / gamma_M2
	if OUT != 0 {
		fmt.Printf("Calculate by Formula 3.13 EN1993-1-8\n")
		fmt.Printf("Tension force by a single row of 3 or more bolts in one leg\n")
		fmt.Printf("Nu_Rd = betta2*Anet*Fu/gamma_M2\n")
		fmt.Printf(("Nu_Rd = %.3f*%.1f*%.1f/%.2f\n\x00"), betta3, Anet*1e+06, Fu*1e-06, gamma_M2)
		fmt.Printf(("Nu_Rd = %.1f kN\n\x00"), Nu_Rd*0.001)
	}
}

func Table_6_5_Leff_group(Leff1 float64, Leff2 float64, Paramaters []float64, Row_table int, OUT int) {
	var m float64 = Paramaters[0]
	var p float64 = Paramaters[1]
	var _alpha float64 = Paramaters[2]
	var e float64 = Paramaters[3]
	var e1 float64 = Paramaters[4]
	{
		var i int
		for ; i < 5; i++ {
			if Paramaters[i] < 0 {
				/// Table 6.5
				///
				fmt.Printf("Stupid check in Table_6_5_Leff_group\n")
				{
					var j int
					for ; j < 5; j++ {
						fmt.Printf(("Paramaters[%d] = %.1f mm\x00"), j, Paramaters[j]*1000)
					}
				}
				panic("ddd")
			}
		}
	}
	var Leff_cp float64
	var Leff_nc float64
	switch Row_table {
	case (1):
		Leff_cp = math.Pi*m + p
		Leff_nc = 0.5*p + _alpha*m - (2*m + 0.625*e)
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Bolt-row adjacent to a stiffener\n")
			fmt.Printf(" Leff_cp = PI*m+p\n")
			fmt.Printf((" Leff_cp = PI*%.1f+%.1f\n\x00"), m*1000, p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = 0.5*p+_alpha*m-(2*m+0.625*e)\n")
			fmt.Printf((" Leff_nc = 0.5*%.1f+%.2f*%.1f-(2*%.1f+0.625*%.1f)\n\x00"), p*1000, _alpha, m*1000, m*1000, e*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	case (2):
		Leff_cp = 2 * p
		Leff_nc = p
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other inner bolt-row\n")
			fmt.Printf(" Leff_cp = 2*p\n")
			fmt.Printf((" Leff_cp = 2*%.1f\n\x00"), p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = p\n")
			fmt.Printf((" Leff_nc = %.1f\n\x00"), p*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	case (3):
		Leff_cp = math.Min(math.Pi*m+p, 2*e1+p)
		Leff_nc = math.Min(2*m+0.625*e+0.5*p, e1+0.5*p)
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other end bolt-row\n")
			fmt.Printf(" Leff_cp = math.Min(PI*m+p,2*e1+p)\n")
			fmt.Printf((" Leff_cp = math.Min(PI*%.1f+%.1f,2*%.1f+%.1f)\n\x00"), m*1000, p*1000, e1*1000, p*1000)
			fmt.Printf((" Leff_cp = %.1f mm\n\x00"), Leff_cp*1000)
			fmt.Printf(" Leff_nc = math.Min(2*m+0.625*e+0.5*p,e1+0.5*p)\n")
			fmt.Printf((" Leff_nc = math.Min(2*%.1f+0.625*%.1f+0.5*%.1f,%.1f+0.5*%.1f)\n\x00"), m*1000, e*1000, p*1000, e1*1000, p*1000)
			fmt.Printf((" Leff_nc = %.1f mm\n\x00"), Leff_nc*1000)
		}
	case (4):
		fmt.Printf(("NOT RELEVANT\x00"))
		panic("ddd")
	default:
		fmt.Printf(("Error in Leff_table6_5_group\x00"))
		panic("ddd")
	}
	Leff1 = math.Min(Leff_nc, Leff_cp)
	if OUT != 0 {
		fmt.Printf(("Leff1 = math.Min(Leff_nc,Leff_cp) = math.Min(%.1f,%.1f)\n\x00"), Leff_nc*1000, Leff_cp*1000)
	}
	Leff2 = Leff_nc
	if OUT != 0 {
		fmt.Printf(("Leff1 = %.1f mm\n\x00"), Leff1*1000)
	}
	if OUT != 0 {
		fmt.Printf(("Leff2 = %.1f mm\n\x00"), Leff2*1000)
	}
}

func Table_6_5_Leff_individ(Leff1 float64, Leff2 float64, Paramaters []float64, Row_table int, OUT int) {
	var DEBUG int
	var e float64 = Paramaters[0]
	var e1 float64 = Paramaters[1]
	var _alpha float64 = Paramaters[2]
	var m float64 = Paramaters[3]
	{
		var i int
		for ; i < 6; i++ {
			if Paramaters[i] < 0 {
				fmt.Printf("Stupid check in Table_6_5_Leff_individ\n")
				{
					var j int
					for ; j < 6; j++ {
						fmt.Printf(("Paramaters[%d] = %.1f mm\n\x00"), j, Paramaters[j]*1000)
					}
				}
				panic("ddd")
			}
		}
	}
	if DEBUG != 0 {
		{
			var i int
			for ; i < 6; i++ {
				fmt.Printf(("Paramaters[%d] = %.1f mm\n\x00"), i, Paramaters[i]*1000)
			}
		}
	}
	var Leff_cp float64
	var Leff_nc float64
	switch Row_table {
	case (1):
		Leff_cp = 2 * math.Pi * m
		Leff_nc = _alpha * m
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Bolt-row adjacent to a stiffener\n")
			fmt.Printf((" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = alpha * m = %.2f * %.1f = %.1f mm\n\x00"), _alpha, m*1000, Leff_nc*1000)
		}
	case (2):
		Leff_cp = 2 * math.Pi * m
		Leff_nc = 4*m + 1.25*e
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other inner bolt-row\n")
			fmt.Printf((" Leff_cp = 2*PI*m = 2*PI * %.1f = %.1f mm\n\x00"), m*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), m*1000, e*1000, Leff_nc*1000)
		}
	case (3):
		var Leff_cp1 float64 = 2 * math.Pi * m
		var Leff_cp2 float64 = math.Pi*m + 2*e1
		Leff_cp = math.Min(Leff_cp1, Leff_cp2)
		var Leff_nc1 float64 = 4*m + 1.25*e
		var Leff_nc2 float64 = e1 + 2*m + 0.625*e
		Leff_nc = math.Min(Leff_nc1, Leff_nc2)
		if OUT != 0 {
			fmt.Printf("Bolt-row location: Other end bolt-row\n")
			fmt.Printf((" Leff_cp1 = 2*PI*m = 2*PI* %.1f = %.1f mm\n\x00"), m*1000, Leff_cp1*1000)
			fmt.Printf((" Leff_cp2 = PI*m+2*e1 = PI* %.1f + 2* %.1f = %.1f mm\n\x00"), m*1000, e1*1000, Leff_cp2*1000)
			fmt.Printf((" Leff_cp = math.Min(Leff_cp1,Leff_cp2) = math.Min(%.1f,%.1f) = %.1f mm\x00"), Leff_cp1*1000, Leff_cp2*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc1 = 4*m+1.25*e = 4*%.1f + 1.25*%.1f = %.1f mm\n\x00"), m*1000, e*1000, Leff_nc1*1000)
			fmt.Printf((" Leff_nc2 = e1+2*m+0.625*e = %.1f + 2* %.1f + 0,625* %.1f = %.1f mm\n\x00"), e1*1000, m*1000, e*1000, Leff_nc2*1000)
			fmt.Printf((" Leff_nc = math.Min(%.1f,%.1f) = %.1f mm\n\x00"), Leff_nc1*1000, Leff_nc2*1000, Leff_nc*1000)
		}
	case (4):
		var Leff_cp1 float64 = 2 * math.Pi * m
		var Leff_cp2 float64 = math.Pi*m + 2*e1
		Leff_cp = math.Min(Leff_cp1, Leff_cp2)
		var Leff_nc float64 = e1 + _alpha*m - (2*m + 0.625*e)
		if OUT != 0 {
			fmt.Printf("Bolt-row location: End bolt-row adjacent to a stiffener\n")
			fmt.Printf((" Leff_cp1 = 2*PI*m = 2*PI* %.1f = %.1f mm\n\x00"), m*1000, Leff_cp1*1000)
			fmt.Printf((" Leff_cp2 = PI*m+2*e1 = PI* %.1f + 2* %.1f = %.1f mm\n\x00"), m*1000, e1*1000, Leff_cp2*1000)
			fmt.Printf((" Leff_cp = math.Min(Leff_cp1,Leff_cp2) = math.Min(%.1f,%.1f) = %.1f mm\x00"), Leff_cp1*1000, Leff_cp2*1000, Leff_cp*1000)
			fmt.Printf((" Leff_nc = e1+alpha*m-2*m+0.625*e = %.1f + %.3f*%.1f- 2* %.1f + 0,625* %.1f = %.1f mm\n\x00"), e1*1000, _alpha, m*1000, m*1000, e*1000, Leff_nc*1000)
		}
	default:
		fmt.Printf(("Error in Leff_table6_6_individ\x00"))
		panic("ddd")
	}
	Leff1 = math.Min(Leff_nc, Leff_cp)
	if OUT != 0 {
		fmt.Printf(("Leff1 = math.Min(Leff_nc,Leff_cp) = math.Min(%.1f,%.1f)\n\x00"), Leff_nc*1000, Leff_cp*1000)
	}
	Leff2 = Leff_nc
	if OUT != 0 {
		fmt.Printf(("Leff1 = %.1f mm\n\x00"), Leff1*1000)
	}
	if OUT != 0 {
		fmt.Printf(("Leff2 = %.1f mm\n\x00"), Leff2*1000)
	}
}
