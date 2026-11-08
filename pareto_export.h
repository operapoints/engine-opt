#ifndef PARETO_EXPORT_H
#define PARETO_EXPORT_H
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

// Feasibility is established against the original problem before constructing
// exportable rows. No penalty-valued or infeasible fallback is exported.
struct ParetoPoint {
    std::string config;
    double min_isp_s = 0., isp_s = 0., mass_kg = 0., thrust_n = 0., fuel_kg_s = 0.;
    std::vector<double> physical;
};
inline bool valid_point(const ParetoPoint& p) {
    const std::size_t size = p.config == "aftfan" ? 20 : p.config == "turbojet" ? 12 : 0;
    if (!size || p.physical.size() != size) return false;
    for (double value : {p.min_isp_s, p.isp_s, p.mass_kg, p.thrust_n, p.fuel_kg_s})
        if (!std::isfinite(value) || value <= 0.) return false;
    return std::all_of(p.physical.begin(), p.physical.end(),
                       [](double v) { return std::isfinite(v); });
}
inline bool dominates(const ParetoPoint& a, const ParetoPoint& b) {
    return a.config == b.config && a.mass_kg <= b.mass_kg && a.isp_s >= b.isp_s
        && (a.mass_kg < b.mass_kg || a.isp_s > b.isp_s);
}
inline std::vector<ParetoPoint> nondominated(std::vector<ParetoPoint> points) {
    points.erase(std::remove_if(points.begin(), points.end(),
                                [](const ParetoPoint& p) { return !valid_point(p); }), points.end());
    std::vector<ParetoPoint> front;
    for (std::size_t i = 0; i < points.size(); ++i) {
        bool discard = false;
        for (std::size_t j = 0; j < points.size(); ++j) {
            if (dominates(points[j], points[i]) || (j < i
                && points[j].config == points[i].config
                && points[j].mass_kg == points[i].mass_kg
                && points[j].isp_s == points[i].isp_s)) { discard = true; break; }
        }
        if (!discard) front.push_back(points[i]);
    }
    std::sort(front.begin(), front.end(), [](const ParetoPoint& a, const ParetoPoint& b) {
        if (a.config != b.config) return a.config < b.config;
        if (a.isp_s != b.isp_s) return a.isp_s < b.isp_s;
        return a.mass_kg < b.mass_kg;
    });
    return front;
}
inline void write_pareto_csv(const std::string& filename,
                             const std::vector<ParetoPoint>& points, bool front_only) {
    std::ofstream out(filename, std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot write CSV: " + filename);
    out << std::setprecision(17)
        << "config,min_isp_s,isp_s,mass_kg,thrust_n,fuel_kg_s";
    for (int i = 1; i <= 20; ++i) out << ",x" << i;
    out << '\n';
    const auto rows = front_only ? nondominated(points) : points;
    for (const auto& p : rows) {
        if (!valid_point(p)) continue;
        out << p.config << ',' << p.min_isp_s << ',' << p.isp_s << ',' << p.mass_kg
            << ',' << p.thrust_n << ',' << p.fuel_kg_s;
        for (std::size_t i = 0; i < 20; ++i) {
            out << ',';
            if (i < p.physical.size()) out << p.physical[i];
        }
        out << '\n';
    }
    out.flush();
    if (!out) throw std::runtime_error("Failed writing CSV: " + filename);
}
#endif