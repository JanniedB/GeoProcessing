//
// Created by Jannie De Beer on 2022/09/04.
//


#include "ZAF_ElevationModel.h"

namespace fs = std::filesystem;

namespace DEM {
    ZAF_ElevationModel::ZAF_ElevationModel() {

        // Initialise common variables for Gauss (Hartbeesthoek 94) to WGS84 transformations
        g_t_vars.e2 = (pow(g_t_vars.a, 2) - pow(g_t_vars.b, 2)) / pow(g_t_vars.a, 2);
        g_t_vars.inv_e2 = (pow(g_t_vars.a, 2) - pow(g_t_vars.b, 2)) / pow(g_t_vars.b, 2);
        g_t_vars.n = (g_t_vars.a - g_t_vars.b) / (g_t_vars.a + g_t_vars.b);
        g_t_vars.p2 = 3.0 / 2.0 * (g_t_vars.n) - 27.0 / 32.0 * pow(g_t_vars.n, 3) + 269.0 / 512.0 * pow(g_t_vars.n, 5);
        g_t_vars.p4 = 21.0 / 16.0 * pow(g_t_vars.n, 2) - 55.0 / 32.0 * pow(g_t_vars.n, 4);
        g_t_vars.p6 = 151.0 / 96.0 * pow(g_t_vars.n, 3) - 417.0 / 128.0 * pow(g_t_vars.n, 5);
        g_t_vars.p8 = 1097.0 / 512.0 * pow(g_t_vars.n, 4);
        g_t_vars.p10 = 8011.0 / 2560.0 * pow(g_t_vars.n, 5);
    }

    void ZAF_ElevationModel::Hart94_to_WGS84(int cm, double y, double x, double &lat, double &lng) {

        // Compute Footprint Lattitude
        double a0 = 1.0 / (g_t_vars.n + 1.0) *
                    (1.0 + 1.0 / 4.0 * pow(g_t_vars.n, 2) +
                     1.0 / 64.0 * pow(g_t_vars.n, 4));
        double foot_bar = x / (g_t_vars.a * a0);
        double foot = foot_bar + g_t_vars.p2 * sin(2.0 * foot_bar) +
                      g_t_vars.p4 * sin(4.0 * foot_bar) +
                      g_t_vars.p6 * sin(6.0 * foot_bar) +
                      g_t_vars.p8 * sin(8.0 * foot_bar) +
                      g_t_vars.p10 * sin(10.0 * foot_bar);

        // Prime Vertical Raduis of Curvature
        double Nf = g_t_vars.a / sqrt(1.0 - g_t_vars.e2 * pow(sin(foot), 2));

        // Use this Value N to evaluate the following     coefficients (like N, all functions of Lat)
        double b1 = 1.0 / (Nf * cos(foot));
        double b2 = tan(foot) / (2 * pow(Nf, 2) * cos(foot));
        double b3 = (1.0 + 2.0 * pow(tan(foot), 2) + g_t_vars.inv_e2 * pow(cos(foot), 2)) /
                    (6.0 * pow(Nf, 3) * cos(foot));
        double b4 = (tan(foot) *
                     (5.0 + 6.0 * pow(tan(foot), 2) + g_t_vars.inv_e2 * pow(cos(foot), 2))) /
                    (24.0 * pow(Nf, 4) * cos(foot));
        double b5 = (5.0 + 28.0 * pow(tan(foot), 2) + 24.0 * pow(tan(foot), 4)) / (120.0 * pow(Nf, 5) * cos(foot));

        double d1 = cos(foot) * (1.0 + g_t_vars.inv_e2 * pow(cos(foot), 2));
        double d2 = -1.0 / 2.0 * pow(cos(foot), 2) * tan(foot) *
                    (1.0 + 4.0 * g_t_vars.inv_e2 * pow(cos(foot), 2));

        // Calculate Latitude and longitude in radians
        double rad_lat = -(foot - b2 * d1 * pow(y, 2) + (b4 * d1 + pow(b2, 2) * d2) * pow(y, 4));
        double rad_lng = degToRad(cm) - (b1 * y - b3 * pow(y, 3) + b5 * pow(y, 5));

        // Calculate Latitude and longitude in degrees (dec.)
        lat = radToDeg(rad_lat);
        lng = radToDeg(rad_lng);
    }


    unsigned long ORT::read_file(std::string ort_filename) {
        std::string line_in;
        std::ifstream ort_file(ort_filename);
        std::vector<ort_rec_type> ort_records;
        long t_x = 0, t_y = 0;
        ort_matrix_info ort_info;
        long last_x = INT32_MAX;
        float t_height;
        double lat, lng;

        const boost::regex ort_line_expr{R"(^ *([0-9]*)\.?[0-9]* *([0-9]*)\.?[0-9]* *([0-9]*\.?[0-9]*))"};
        boost::smatch _match_results;

        try {
            fileName = fs::path(ort_filename).filename();
            filePath = fs::path(fileName).parent_path();
            fileExtention = fs::path(ort_filename).extension();

            if (ort_file.is_open()) {
                while (getline(ort_file, line_in)) {
                    // extract the record components
                    if (boost::regex_search(line_in, _match_results, ort_line_expr)) {
                        if (_match_results.size() == 4) {

                            t_y = std::stoi(_match_results[1]);
                            t_x = std::stoi(_match_results[2]);
                            t_height = std::stof(_match_results[3]);


                            DEM::ZAF_ElevationModel::Hart94_to_WGS84(23, t_y, t_x, lat, lng);

                            std::cout << std::setprecision(9) << "south:" << t_y << " east:" << t_x << " latitude: " << lat << ", longitude: " << lng
                                      << std::endl;

                            ort_info.min_x = (ort_info.min_x < t_x) ? ort_info.min_x : t_x;
                            ort_info.min_y = (ort_info.min_y < t_y) ? ort_info.min_y : t_y;
                            ort_info.max_x = (ort_info.max_x > t_x) ? ort_info.max_x : t_x;
                            ort_info.max_y = (ort_info.max_y > t_y) ? ort_info.max_y : t_y;


                            if (last_x == INT32_MAX)
                                last_x = t_x;
                            else if (last_x != t_x) {
                                last_x = t_x;
                                ort_info.num_rows++;
                                ort_info.num_cols = 0;
                            } else
                                ort_info.num_cols++;

                            ort_records.push_back({t_x, t_y, t_height});
                            ort_info.num_records++;
                        }

                    }
                }
                ort_file.close();


                std::cout << "min_x:" << ort_info.min_x << " max_x:" << ort_info.max_x << " min_y:" << ort_info.min_y
                          << " max_y:" << ort_info.max_y
                          << " Num Records:" << ort_info.num_records << " size: " << ort_info.num_rows + 1 << "x"
                          << ort_info.num_cols + 1
                          << std::endl;

                if ((ort_info.num_cols > 0)) {
                    // Create memory block for all the data records
                    float data_out[ort_info.num_rows + 1][ort_info.num_cols + 1];

                    std::cout << "Size of matrix" << sizeof(data_out) << std::endl;

                    memset(data_out, 0, sizeof(data_out));

                    for (auto &element: ort_records) {
//                        std::cout << "(" << (element.x - ort_info.min_x) / 25 << "/"
//                                  << (element.y - ort_info.min_y) / 25 << ") " << element.x << ":" << element.y << "="
//                                  << element.height << std::endl;
                        data_out[(element.x - ort_info.min_x) / 25][(element.y - ort_info.min_y) / 25] = element.height;

                    }
                    std::cout << ".. done .." << std::endl;

                }

            } else {
                spdlog::error("Cannot open file {}{}", fileName, fileExtention);
                return 0;
            }
        }

        catch (const std::overflow_error &e) { // this executes if f() throws std::overflow_error (same type rule)
            spdlog::error(e.what());
        }
        catch (const std::runtime_error &e) { // this executes if f() throws std::underflow_error (base class rule)
            spdlog::error(e.what());
        }
        catch (const std::exception &e) { // this executes if f() throws std::logic_error (base class rule)
            spdlog::error(e.what());
        }
        catch (...) { // this executes if f() throws std::string or int or any other unrelated type
            spdlog::error("Unknown Error");
        }


        return ort_records.size();
    }

    int ORT::serialise_file() {
        return 0;
    }

    int ORT::de_serialise_file() {
        return 0;
    }

    ORT::ORT() {

        ort_instance_count++;

        if (latitude_mapping.empty())
            for (int i = 16; i <= 33; i += 2) {
                latitude_mapping[i] = latitude_mapping[i + 1] = i + 1;
            }

        // Open SQLite database if it was not opened already
        if (sql_db == nullptr)
            if (sqlite3_open("identifier.sqlite", &sql_db)) {
                spdlog::error("Can't open database: {}", sqlite3_errmsg(sql_db));
            } else {
                spdlog::debug("SQLlite database opened successfully");
            }


    }

    ORT::~ORT() {
        ort_instance_count--;
        if (ort_instance_count == 0 && sql_db != nullptr)
            sqlite3_close(sql_db);

    }


} // DEM