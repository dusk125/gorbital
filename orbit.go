package gorbital

import (
	"math"
	"time"
)

/*
This is a Golang port of the Javascript Orrery from http://www.planetaryorbits.com/tutorial-javascript-orbit-simulation.html
*/

type Params struct {
	SemiMajorAxis          float64 // Semi-Major Axis (a) [AU]
	Eccentricity           float64 // Eccentricity (e) [radians]
	Inclination            float64 // Inclination (I) [degrees]
	MeanLongitude          float64 // Mean Longitude (L) [degrees]
	LongitudePerihelion    float64 // Longitude of Perihelion (w) [degrees]
	LongitudeAscendingNode float64 // Longitude of Ascending Node (W) [degrees]
	DSemiMajorAxis         float64 // Per century change in Semi-Major Axis (a) [AU/century]
	DEccentricity          float64 // Per century change in Eccentricity (e) [radians/century]
	DInclination           float64 // Per century change in Inclination (I) [degrees/century]
	DMeanLongitude         float64 // Per century change in Mean Longitude (L) [degrees/century]
	DLongPerihelion        float64 // Per century change in Longitude of Perihelion (w) [degrees/century]
	DLongAscendingNode     float64 // Per century change in Longitude of Ascending Node (W) [degrees/century]
}

const (
	secondsPerYear      = 86400.0
	julianUnixEpoch     = 2440587.5
	JulianEpochJ2000    = 2451545.0
	JulianCenturyInDays = 36525
)

// Julian converts a Golang time.Time structure to a Julian timestamp
func Julian(t time.Time) float64 {
	return float64(t.Unix()/secondsPerYear) + julianUnixEpoch
}

func JulianCenturiesSinceJ2000(t time.Time) float64 {
	return (Julian(t) - JulianEpochJ2000) / JulianCenturyInDays
}

// Converts radians to degrees
func degrees(radians float64) float64 {
	return radians * 180 / math.Pi
}

// Converts degrees to radians
func radians(degrees float64) float64 {
	return degrees * math.Pi / 180
}

// EstimateEccentricAnomaly calculates the eccentric anomaly (E) in degrees
func EstimateEccentricAnomaly(eccentricity, meanAnomaly float64) float64 {
	meanAnomaly = meanAnomaly / 360
	meanAnomaly = 2 * math.Pi * (meanAnomaly - math.Floor(meanAnomaly))

	eccentricAnomaly := math.Pi
	if eccentricity < 0.8 {
		eccentricAnomaly = meanAnomaly
	}

	F := eccentricAnomaly - eccentricity*math.Sin(meanAnomaly) - meanAnomaly

	for i := 0; math.Abs(F) > 10e-6 && i < 30; i++ {
		eccentricAnomaly = eccentricAnomaly - F/(1-eccentricity*math.Cos(eccentricAnomaly))
		F = eccentricAnomaly - eccentricity*math.Sin(eccentricAnomaly) - meanAnomaly
	}

	return degrees(eccentricAnomaly)
}

// Decompose calculates the 3D-cartisians coordinates for the given orbital parameters at a given time
func (p Params) Decompose(t time.Time) (x, y, z float64) {
	T := JulianCenturiesSinceJ2000(t)
	p = Params{
		SemiMajorAxis:          p.SemiMajorAxis + p.DSemiMajorAxis*T,
		Eccentricity:           p.Eccentricity + p.DEccentricity*T,
		Inclination:            math.Mod(p.Inclination+p.DInclination*T, 360),
		MeanLongitude:          math.Mod(p.MeanLongitude+p.DMeanLongitude*T, 360),
		LongitudePerihelion:    math.Mod(p.LongitudePerihelion+p.DLongPerihelion*T, 360),
		LongitudeAscendingNode: math.Mod(p.LongitudeAscendingNode+p.DLongAscendingNode*T, 360),
	}

	if p.LongitudePerihelion < 0 {
		p.LongitudePerihelion = 360 + p.LongitudePerihelion
	}
	if p.MeanLongitude < 0 {
		p.MeanLongitude = 360 + p.MeanLongitude
	}

	meanAnomaly := p.MeanLongitude - p.LongitudePerihelion
	if meanAnomaly < 0 {
		meanAnomaly = 360 + meanAnomaly
	}

	eccentricAnomaly := EstimateEccentricAnomaly(p.Eccentricity, meanAnomaly)

	argTrueAnomaly := (math.Sqrt((1 + p.Eccentricity) / (1 - p.Eccentricity))) * (math.Tan(radians(eccentricAnomaly) / 2))
	var n float64
	if argTrueAnomaly < 0 {
		n = 2 * (degrees(math.Atan(argTrueAnomaly)) + 180)
	} else {
		n = 2 * degrees(math.Atan(argTrueAnomaly))
	}

	radius := p.SemiMajorAxis * (1 - (p.Eccentricity * math.Cos(radians(eccentricAnomaly))))

	radLAN := radians(p.LongitudeAscendingNode)
	radIncl := radians(p.Inclination)
	long := radians(n + p.LongitudePerihelion - p.LongitudeAscendingNode)
	x = radius * (math.Cos(radLAN)*math.Cos(long) - math.Sin(radLAN)*math.Sin(long)*math.Cos(radIncl))
	y = radius * (math.Sin(radLAN)*math.Cos(long) + math.Cos(radLAN)*math.Sin(long)*math.Cos(radIncl))
	z = radius * (math.Sin(long) * math.Sin(radIncl))

	return
}
