//
// Created by Jannie De Beer on 2022/09/04.
//

#ifndef GEOPROCESSING_ZAF_ELEVATIONMODEL_H
#define GEOPROCESSING_ZAF_ELEVATIONMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <locale>
#include <filesystem>
#include <vector>
#include <map>
#include "sqlite3.h"
#include <spdlog/spdlog.h>
#include <boost/regex.hpp>
#include <iostream>
#include <pqxx/pqxx>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)
#define radToDeg(angleInRadians) ((angleInRadians) * 180.0 / M_PI)

namespace DEM {
    // trim from start (in place)
    static inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    static inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    // trim from both ends (in place)
    static inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }

    // trim from start (copying)
    static inline std::string ltrim_copy(std::string s) {
        ltrim(s);
        return s;
    }

    // trim from end (copying)
    static inline std::string rtrim_copy(std::string s) {
        rtrim(s);
        return s;
    }

    // trim from both ends (copying)
    static inline std::string trim_copy(std::string s) {
        trim(s);
        return s;
    }

    class ORT {
        struct ort_rec_type {
            long x;
            long y;
            float height;
        };
        struct ort_matrix_info {
            unsigned long min_x = INT32_MAX,
                    max_x = 0,
                    min_y = INT32_MAX,
                    max_y = 0;
            unsigned long num_records = 0;
            unsigned long num_cols = 0,
                    num_rows = 0;
        };

        std::string fileName;
        std::string fileExtention;
        std::string filePath;
        unsigned char latitude = 0, longitude = 0;
        static sqlite3 *sql_db;
        static int ort_instance_count;

        static std::map<char8_t, char8_t> latitude_mapping;



    public:
        ORT();

        unsigned long read_file(std::string filename);
        int serialise_file();
        int de_serialise_file();

        virtual ~ORT();


    };

    class ZAF_ElevationModel {
    private:
        struct common_gaussian {
            const double a = 6378137;       // Semi major axis (a)
            const double b = 6356752.314;   // Semi minor axis (b)
            double e2, inv_e2, n;
            double p2, p4, p6, p8, p10;
        };

        struct ort_rec_type {
            long x;
            long y;
            float height;
        };

    public:

        ZAF_ElevationModel();

        static void Hart94_to_WGS84(int cm, double y, double x, double &lat, double &lng);

        static common_gaussian g_t_vars;

        long read_ORT_file(std::string filename);
    };

} // DEM

#endif //GEOPROCESSING_ZAF_ELEVATIONMODEL_H
