//
// Created by Jannie De Beer on 2022/09/03.
//
#include <iostream>
#include "gauss-kruger/gausskruger.h"

class GRS80 : public gausskruger::Projection
{
public:
    double flattening() { return 1 / 298.257222101; }
    double equatorialRadius() { return 6378137.0; }
};

class RT90_75_gon_V : public GRS80
{
public:
    double centralMeridian() { return 11 + 18.375 / 60.0; }
    double scale() { return 1.000006; }
    double falseNorthing() { return -667.282; }
    double falseEasting() { return 1500025.141; }
};

int main(int argc, char* argv[])
{
    RT90_75_gon_V projection;
    double lat, lon, northing, easting;
    lat = -33.80479656;
    lon = 18.50534292;
    projection.geodeticToGrid(lat, lon, northing, easting);
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3);
    std::cout << "Northing: " << northing << "\nEasting: " << easting << std::endl;


    return 0;
}
