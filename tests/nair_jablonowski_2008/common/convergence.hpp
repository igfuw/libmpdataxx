#pragma once

#include <fstream>
#include "../../mp3_paper_2018_JCP/moving_vort/transforms.hpp"

template <typename Y = T>
struct stat_t
{
  Y min, max, L1, L2, Li, ord_L1, ord_L2, ord_Li;
  stat_t() : ord_L1{0}, ord_L2{0}, ord_Li{0} {}
};

template <typename Y>
std::ostream& operator<<(std::ostream &os, const stat_t<Y> &s)
{
  auto precise = std::setprecision(6);
  auto rough = std::setprecision(2);
  auto swp = std::setw(10);
  auto swr = std::setw(6);
  return os << std::fixed << precise
            << swp << s.L1 << swp << s.L2 << swp << s.Li
            << rough
            << swr << s.ord_L1 << swr << s.ord_L2 << swr << s.ord_Li
            << precise
            << swp << s.min << swp << s.max;
}

template <typename test_t>
void convergence(test_t test, const std::string &base_name, const T max_cfl)
{
    std::vector<int> nys = { 24, 48, 96};
    std::vector<stat_t<>> stats(nys.size());
    auto test_ny = [test, max_cfl, base_name](const int ny) { return test(base_name, ny, max_cfl);};
    std::transform(nys.begin(), nys.end(), stats.begin(), test_ny);
    
    // calculate pairwise approximations to convergence order in each norm
    for (int i = 1; i < stats.size(); ++i)
    {
      stats[i].ord_L1 = log(stats[i-1].L1 / stats[i].L1) / log(2.0);
      stats[i].ord_L2 = log(stats[i-1].L2 / stats[i].L2) / log(2.0);
      stats[i].ord_Li = log(stats[i-1].Li / stats[i].Li) / log(2.0);
    }

    // stats output
    stat_t<std::string> header;
    header.L1 = "L1";
    header.L2 = "L2";
    header.Li = "Li";
    header.ord_L1 = "oL1";
    header.ord_L2 = "oL2";
    header.ord_Li = "oLi";
    header.min = "min";
    header.max = "max";

    std::ofstream stats_ofs("stats_" + base_name + ".txt");
    stats_ofs << std::setw(4) << "ny" << " " << header << std::endl;
    for (int i = 0; i < stats.size(); ++i)
    {
      stats_ofs << std::setw(4) << nys[i] << " " << stats[i] << std::endl;
    }
}
