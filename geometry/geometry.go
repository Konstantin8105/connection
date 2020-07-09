package geometry

type XY struct {
	X, Y float64
}

func NewPlate(widthX, widthY float64) []XY {
	return []XY{
		XY{+widthX / 2.0, +widthY / 2.0},
		XY{-widthX / 2.0, +widthY / 2.0},
		XY{+widthX / 2.0, -widthY / 2.0},
		XY{-widthX / 2.0, -widthY / 2.0},
	}
}
