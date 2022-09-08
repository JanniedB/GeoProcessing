#include <iostream>
#include <iomanip>
#include <spdlog/spdlog.h>
#include "ZAF_ElevationModel.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)
#define radToDeg(angleInRadians) ((angleInRadians) * 180.0 / M_PI)

void Hart94_to_WGS84(int cm, double y, double x, double &lat, double &lng) {
    const double a = 6378137;       // Semi major axis (a)
    const double b = 6356752.314;   // Semi minor axis (b)

//    and may be computed once and for all
//2. Eccentricities
    double e2 = (pow(a, 2) - pow(b, 2)) / pow(a, 2);
    double inv_e2 = (pow(a, 2) - pow(b, 2)) / pow(b, 2);
    double n = (a - b) / (a + b);

//3. Inverse Meridian Arc Length
    double p2 = 3.0 / 2.0 * (n) - 27.0 / 32.0 * pow(n, 3) + 269.0 / 512.0 * pow(n, 5);
    double p4 = 21.0 / 16.0 * pow(n, 2) - 55.0 / 32.0 * pow(n, 4);
    double p6 = 151.0 / 96.0 * pow(n, 3) - 417.0 / 128.0 * pow(n, 5);
    double p8 = 1097.0 / 512.0 * pow(n, 4);
    double p10 = 8011.0 / 2560.0 * pow(n, 5);

//4. Compute Footprint Lattitude
    double a0 = 1.0 / (n + 1.0) * (1.0 + 1.0 / 4.0 * pow(n, 2) + 1.0 / 64.0 * pow(n, 4));
    double foot_bar = x / (a * a0);
    double foot = foot_bar + p2 * sin(2.0 * foot_bar) + p4 * sin(4.0 * foot_bar) + p6 * sin(6.0 * foot_bar) +
                  p8 * sin(8.0 * foot_bar) + p10 * sin(10.0 * foot_bar);

//5. Prime Vertical Raduis of Curvature
    double Nf = a / sqrt(1.0 - e2 * pow(sin(foot), 2));

//5. Use this Value N to evaluate the following     coefficients (like N, all functions of Lat)
    double b1 = 1.0 / (Nf * cos(foot));
    double b2 = tan(foot) / (2 * pow(Nf, 2) * cos(foot));
    double b3 = (1.0 + 2.0 * pow(tan(foot), 2) + inv_e2 * pow(cos(foot), 2)) / (6.0 * pow(Nf, 3) * cos(foot));
    double b4 = (tan(foot) * (5.0 + 6.0 * pow(tan(foot), 2) + inv_e2 * pow(cos(foot), 2))) /
                (24.0 * pow(Nf, 4) * cos(foot));
    double b5 = (5.0 + 28.0 * pow(tan(foot), 2) + 24.0 * pow(tan(foot), 4)) / (120.0 * pow(Nf, 5) * cos(foot));

    double d1 = cos(foot) * (1.0 + inv_e2 * pow(cos(foot), 2));
    double d2 = -1.0 / 2.0 * pow(cos(foot), 2) * tan(foot) * (1.0 + 4.0 * inv_e2 * pow(cos(foot), 2));

//6. Calculate Latitude and longitude in radians
    double rad_lat = -(foot - b2 * d1 * pow(y, 2) + (b4 * d1 + pow(b2, 2) * d2) * pow(y, 4));
    double rad_lng = degToRad(cm) - (b1 * y - b3 * pow(y, 3) + b5 * pow(y, 5));

//7. Calculate Latitude and longitude in degrees (dec.)
    lat = radToDeg(rad_lat);
    lng = radToDeg(rad_lng);

}

// Initialise Static Variables
DEM::ZAF_ElevationModel::common_gaussian DEM::ZAF_ElevationModel::g_t_vars;
std::map<char8_t, char8_t> DEM::ORT::latitude_mapping;
sqlite3 *DEM::ORT::sql_db = nullptr;
int DEM::ORT::ort_instance_count = 0;

int main() {
//    spdlog::info("Welcome to spdlog!");
//    spdlog::error("Some error message with arg:");

    DEM::ZAF_ElevationModel dem_25;
    DEM::ORT ort;


    double lat, lng;

    std::cout << "Hello, World!" << std::endl;

    DEM::ZAF_ElevationModel::Hart94_to_WGS84(23, 64592.212, 2475220.602, lat, lng);

    std::cout << std::setprecision(10) << "latitude[-22.3728395]: " << lat << ", longitude[22.3728395]: " << lng
              << std::endl;

    DEM::ZAF_ElevationModel::Hart94_to_WGS84(19, 45803.274, 3742119.361, lat, lng);

    std::cout << std::setprecision(10) << "latitude[-33.80479657]: " << lat << ", longitude[18.50534291]: " << lng
              << std::endl;

    auto n = ort.read_file("/Users/janniedebeer/CLionProjects/GeoProcessing/data/25m_DEM/2230dd/2230DD25.ORT");


    return 0;
}
