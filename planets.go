package gorbital

// All elements are at J2000 and rates are per century past J2000
var (
	Mercury = Params{
		SemiMajorAxis:          0.38709927,
		Eccentricity:           0.20563593,
		Inclination:            7.00497902,
		MeanLongitude:          252.25032350,
		LongitudePerihelion:    77.45779628,
		LongitudeAscendingNode: 48.33076593,
		DSemiMajorAxis:         0.00000037,
		DEccentricity:          0.00001906,
		DInclination:           -0.00594749,
		DMeanLongitude:         149472.67411175,
		DLongPerihelion:        0.16047689,
		DLongAscendingNode:     -0.1253408,
	}
	Venus = Params{
		SemiMajorAxis:          0.72333566,
		Eccentricity:           0.00677672,
		Inclination:            3.39467605,
		MeanLongitude:          181.97909950,
		LongitudePerihelion:    131.60246718,
		LongitudeAscendingNode: 76.67984255,
		DSemiMajorAxis:         0.00000390,
		DEccentricity:          -0.00004107,
		DInclination:           -0.00078890,
		DMeanLongitude:         58517.81538729,
		DLongPerihelion:        0.00268329,
		DLongAscendingNode:     -0.27769418,
	}
	Earth = Params{
		SemiMajorAxis:          1.00000261,
		Eccentricity:           0.01671123,
		Inclination:            -0.00001531,
		MeanLongitude:          100.46457166,
		LongitudePerihelion:    102.93768193,
		LongitudeAscendingNode: 0.0,
		DSemiMajorAxis:         0.00000562,
		DEccentricity:          -0.00004392,
		DInclination:           -0.01294668,
		DMeanLongitude:         35999.37244981,
		DLongPerihelion:        0.32327364,
		DLongAscendingNode:     0.0,
	}
	Mars = Params{
		SemiMajorAxis:          1.52371034,
		Eccentricity:           0.09339410,
		Inclination:            1.84969142,
		MeanLongitude:          -4.55343205,
		LongitudePerihelion:    -23.94362959,
		LongitudeAscendingNode: 49.55953891,
		DSemiMajorAxis:         0.00001847,
		DEccentricity:          0.00007882,
		DInclination:           -0.00813131,
		DMeanLongitude:         19140.30268499,
		DLongPerihelion:        0.44441088,
		DLongAscendingNode:     -0.29257343,
	}
	Jupiter = Params{
		SemiMajorAxis:          5.20288700,
		Eccentricity:           0.04838624,
		Inclination:            1.30439695,
		MeanLongitude:          34.39644051,
		LongitudePerihelion:    14.72847983,
		LongitudeAscendingNode: 100.47390909,
		DSemiMajorAxis:         -0.00011607,
		DEccentricity:          -0.00013253,
		DInclination:           -0.00183714,
		DMeanLongitude:         3034.74612775,
		DLongPerihelion:        0.21252668,
		DLongAscendingNode:     0.20469106,
	}
	Saturn = Params{
		SemiMajorAxis:          19.18916464,
		Eccentricity:           0.04725744,
		Inclination:            0.77263783,
		MeanLongitude:          313.23810451,
		LongitudePerihelion:    170.95427630,
		LongitudeAscendingNode: 74.01692503,
		DSemiMajorAxis:         -0.00125060,
		DEccentricity:          -0.00050991,
		DInclination:           0.00193609,
		DMeanLongitude:         1222.49362201,
		DLongPerihelion:        -0.41897216,
		DLongAscendingNode:     -0.28867794,
	}
	Uranus = Params{
		SemiMajorAxis:          19.18916464,
		Eccentricity:           0.04725744,
		Inclination:            0.77263783,
		MeanLongitude:          313.23810451,
		LongitudePerihelion:    170.95427630,
		LongitudeAscendingNode: 74.01692503,
		DSemiMajorAxis:         -0.00196176,
		DEccentricity:          -0.00004397,
		DInclination:           -0.00242939,
		DMeanLongitude:         428.48202785,
		DLongPerihelion:        0.40805281,
		DLongAscendingNode:     0.04240589,
	}
	Neptune = Params{
		SemiMajorAxis:          30.06992276,
		Eccentricity:           0.00859048,
		Inclination:            1.77004347,
		MeanLongitude:          -55.12002969,
		LongitudePerihelion:    44.96476227,
		LongitudeAscendingNode: 131.78422574,
		DSemiMajorAxis:         0.00026291,
		DEccentricity:          0.00005105,
		DInclination:           0.00035372,
		DMeanLongitude:         218.45945325,
		DLongPerihelion:        -0.32241464,
		DLongAscendingNode:     -0.00508664,
	}
)
