package unit

type Meter float64

func (m Meter) ToMM() float64 {
	return float64(m) * 1000.0
}

type Pa float64

func (p Pa) ToMPa() float64 {
	return float64(p) * 1.0e-6
}

type Newton float64

func (n Newton) ToKN() float64 {
	return float64(n) * 1e-3
}
