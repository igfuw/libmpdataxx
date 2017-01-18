#pragma once

#include <fstream>

template <typename Y = T>
struct stat_t
{
  Y min, max, L1, L2, Li, ord_L1, ord_L2, ord_Li;
};

template <typename Y>
std::ostream& operator<<(std::ostream &os, const stat_t<Y> &s)
{
  auto sw = std::setw(10);
  auto precise = std::setprecision(6);
  auto rough = std::setprecision(2);
  return os << std::fixed << precise
            << sw << s.L1 << sw << s.L2 << sw << s.Li
            << rough
            << sw << s.ord_L1 << sw << s.ord_L2 << sw << s.ord_Li
            << precise
            << sw << s.min << sw << s.max;
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
    header.ord_L1 = "ord_L1";
    header.ord_L2 = "ord_L2";
    header.ord_Li = "ord_Li";
    header.min = "min";
    header.max = "max";

    std::ofstream stats_ofs("stats_" + base_name + ".txt");
    stats_ofs << "ny" << " " << header << std::endl;
    for (int i = 0; i < stats.size(); ++i)
    {
      stats_ofs << nys[i] << " " << stats[i] << std::endl;
    }
}
