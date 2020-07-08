package baseplate_test

import (
	"os"

	bolt "github.com/Konstantin8105/Eurocode3.Bolt"
	"github.com/Konstantin8105/connection/baseplate"
	"github.com/Konstantin8105/connection/eu2"
	"github.com/Konstantin8105/connection/eu3"
	"github.com/Konstantin8105/connection/geometry"
	"github.com/Konstantin8105/connection/load"
)

func Example() {
	bp := baseplate.BasePlate{
		B: bolt.New(bolt.D24, bolt.G4p6),
		BPos: []geometry.XY{
			geometry.XY{X: +0.160, Y: +0.160},
			geometry.XY{X: -0.160, Y: +0.160},
			geometry.XY{X: +0.160, Y: -0.160},
			geometry.XY{X: -0.160, Y: -0.160},
		},

		PlatePos: []geometry.XY{
			geometry.XY{X: +0.170, Y: +0.170},
			geometry.XY{X: -0.170, Y: +0.170},
			geometry.XY{X: +0.170, Y: -0.170},
			geometry.XY{X: -0.170, Y: -0.170},
		},
		PlateThk: 0.030,

		SectionPos: []XY{}, // TODO: section.Generate("HEB200"),
		SectionProp: eu3.S235EN10025_2,

		ColumnPos: []geometry.XY{
			geometry.XY{X: +1.600, Y: +1.600},
			geometry.XY{X: -1.600, Y: +1.600},
			geometry.XY{X: +1.600, Y: -1.600},
			geometry.XY{X: -1.600, Y: -1.600},
		},
		ColumnProp: eu2.GradeC12_15,
	}

	bp.Calculate(os.Stdout, load.Load{N: -1000.0e3, V: 0.0, M: 100.0e3})
	// Output:
}
