package baseplate

import (
	"io"

	bolt "github.com/Konstantin8105/Eurocode3.Bolt"
	"github.com/Konstantin8105/connection/eu2"
	"github.com/Konstantin8105/connection/eu3"
	"github.com/Konstantin8105/connection/geometry"
	"github.com/Konstantin8105/connection/load"
)

type Washer struct {
	Pos  []geometry.XY
	Thk  float64
	Prop eu3.Steel
}

// BasePlate is general struct for base plate calculation without
// specific of sections
type BasePlate struct {
	B    bolt.Bolt     // bolts
	BPos []geometry.XY // bolt position

	PlatePos []geometry.XY // points of plate
	PlateThk float64       // thickness of plate

	SectionPos  []geometry.XY // section coordinates with stiffiners
	SectionProp eu3.Steel     // section steel

	Washers []Washer

	ColumnPos  []geometry.XY // concrete colimn size
	ColumnProp eu2.Concrete  // concrete column property
}

func (bp BasePlate) Check() error {
	return nil
}

func (bp BasePlate) Calculate(out io.Writer, forces load.Load) (err error) {
	if err = bp.Check(); err != nil {
		return
	}

	return nil
}
