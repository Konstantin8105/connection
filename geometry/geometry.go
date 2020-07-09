package geometry

import "github.com/Konstantin8105/section"

type XY struct {
	X, Y float64
}

func NewPlate(widthX, widthY float64) []XY {
	return []XY{
		XY{X: +widthX / 2.0, Y: +widthY / 2.0},
		XY{X: -widthX / 2.0, Y: +widthY / 2.0},
		XY{X: +widthX / 2.0, Y: -widthY / 2.0},
		XY{X: -widthX / 2.0, Y: -widthY / 2.0},
	}
}

func NewISection(name string) []XY {
	for _, s := range section.Isections {
		if name != s.Name {
			continue
		}
		return []XY{
			//
			XY{X: 0.000, Y: s.H / 2.0},
			//
			XY{X: s.B / 2.0, Y: s.H / 2.0},
			XY{X: s.B / 2.0, Y: s.H/2.0 - s.Tf},
			XY{X: s.Tw / 2.0, Y: s.H/2.0 - s.Tf},
			XY{X: s.Tw / 2.0, Y: 0.000},
			XY{X: s.Tw / 2.0, Y: -s.H/2.0 + s.Tf},
			XY{X: s.B / 2.0, Y: -s.H/2.0 + s.Tf},
			XY{X: s.B / 2.0, Y: -s.H / 2.0},
			//
			XY{X: 0.000, Y: -s.H / 2.0},
			//
			XY{X: -s.B / 2.0, Y: -s.H / 2.0},
			XY{X: -s.B / 2.0, Y: -s.H/2.0 + s.Tf},
			XY{X: -s.Tw / 2.0, Y: -s.H/2.0 + s.Tf},
			XY{X: -s.Tw / 2.0, Y: 0.000},
			XY{X: -s.Tw / 2.0, Y: s.H/2.0 - s.Tf},
			XY{X: -s.B / 2.0, Y: s.H/2.0 - s.Tf},
			XY{X: -s.B / 2.0, Y: s.H / 2.0},
		}
	}

	panic("not found section : " + name)
}
