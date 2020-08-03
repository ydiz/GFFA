#include <dirent.h>

// FIXME: there might be an "index out of range" error in get_a0() when x is extremly close to 1
// and x_idx_lower == num_high - 1
// But uniform_real_distribution generates random numbers on a half-open range; so it is probably fine.


class Integral_table {
public:
  int num_high;
  int num_low;
  double dividing_line;
  bool initialized = false;
  double high_interval;
  double low_interval;

  double min_k;
  double k_interval;

  std::vector<double> ks;
  std::vector<std::vector<double>> table_high;
  std::vector<std::vector<double>> table_low; // set to 1025 so that interval is 1./1024.
  std::vector<double> x_high;
  std::vector<double> x_low;

  std::string file_prefix;

  // coeff_of_a0 is used to check if it is consistent with coefficient in file
  // example of parameters num_high = 1024, num_low = 1025, dividing_line = 1 / 1024.
  Integral_table(double coeff_of_trUV, const std::string &_file_prefix) : file_prefix(_file_prefix) {

    std::map<double, std::string> files;
    get_ks_files(ks, files);
    // for(auto &x: ks) std::cout << x << "\t" << files[x] << std::endl;

    min_k = ks[0];
    k_interval = ks[1] - ks[0];

    std::ifstream fin(file_prefix + "/"+ files[ks[0]]);
    read_header(fin);
    fin.close();

    table_high.resize(ks.size());
    table_low.resize(ks.size());

    for(int i=0; i<ks.size(); ++i) {
      table_high[i].resize(num_high);
      table_low[i].resize(num_low);
    }

    for(int i=0; i<ks.size(); ++i) {
      double k = ks[i];
      double coeff_of_a0 = 2.0 * coeff_of_trUV * k; // for chekcing if we are reading the correct file

      std::string filename_table_high =  files[k];
      std::string filename_table_low =  files[k];
      filename_table_low.replace(filename_table_low.begin()+6, filename_table_low.begin()+10, "low");

      read_table(table_high[i], filename_table_high, coeff_of_a0);
      read_table(table_low[i], filename_table_low, coeff_of_a0);
    }

    x_high.resize(num_high);
    high_interval = double(1 - dividing_line) / double(num_high - 1); // 1. / 1024.;
    for(int i=0; i<x_high.size(); ++i) x_high[i] = dividing_line + high_interval * i;

    x_low.resize(num_low);
    low_interval = double(dividing_line) / double(num_low - 1); // 1. / 1024. / 1024.;
    for(int i=0; i<x_low.size(); ++i) x_low[i] = low_interval * i;

    std::cout << "===================Integral_table: =========================" << std::endl;
    std::cout << "num_high: " << num_high << std::endl;
    std::cout << "num_low: " << num_low << std::endl;
    std::cout << "dividing_line: " << dividing_line << std::endl;
    for(double x: ks) std::cout << x << "\t"; std::cout << std::endl;
    std::cout << "============================================================" << std::endl;

  }

  void read_table(std::vector<double> &table, const std::string &filename, double coeff_of_a0) {
    std::ifstream fin(file_prefix + "/" + filename);
    double coeff_in_file = read_header(fin);
    assert(std::abs(coeff_in_file - coeff_of_a0) < 1e-5);
    int i = 0;
    while(fin >> table[i]) ++i;
    assert(i==num_low);
    fin.close();
  }

  double read_header(std::ifstream &fin) {
    std::string line;

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    assert(line == "#BEGIN_HEADER");

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    int eq = line.find("=");
    assert(line.substr(0, eq) == "#coefficient");
    double coeff = std::stod(line.substr(eq+1));

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    eq = line.find("=");
    assert(line.substr(0, eq) == "#num_high");
    if(initialized) assert( num_high == std::stoi(line.substr(eq+1)) );
    else num_high = std::stoi(line.substr(eq+1));

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    eq = line.find("=");
    assert(line.substr(0, eq) == "#num_low");
    if(initialized) assert( num_low == std::stoi(line.substr(eq+1)) );
    else num_low = std::stoi(line.substr(eq+1));

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    eq = line.find("=");
    assert(line.substr(0, eq) == "#dividing_line");
    if(initialized) assert( dividing_line == std::stod(line.substr(eq+1)) );
    else dividing_line = std::stod(line.substr(eq+1));

    std::getline(fin, line);
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    assert(line == "#END_HEADER");

    if(!initialized) initialized = true;

    return coeff;
  }

  void get_ks_files(std::vector<double> &ks, std::map<double, std::string> &files) {
    DIR *dir;
    dir = opendir(file_prefix.c_str());
    assert(dir!=NULL); // make sure directory exists
    struct dirent *entry;

    std::string filename;
    while ((entry = readdir (dir)) != NULL) {
      // printf ("%s\n", entry->d_name);
      filename = std::string(entry->d_name);
      if(filename.substr(0, 10) == "table_high") {
        double k = std::stod(filename.substr(13,  filename.find("_beta_") - 13)); // "table_high_k_1.5789_beta_6.0.txt";
        ks.push_back(k);
        files.insert(std::pair<double, std::string>(k, filename));
      }
    }
    closedir (dir);

    std::sort(ks.begin(), ks.end());
  }

  double get_a0(double k, double x) {
    double a0, a0_k_low, a0_k_high;
    int x_idx_lower;
    int k_idx_lower = int( (k - min_k) / k_interval );

    if(x >= dividing_line) {
      x_idx_lower = int( (x - dividing_line) / high_interval );
      // assert(x_idx_lower < num_high); // FIXME: if x==1 -> x_lower==num_high -> index out of range
      a0_k_low = linear_interpolation(x, x_high[x_idx_lower], table_high[k_idx_lower][x_idx_lower], x_high[x_idx_lower+1], table_high[k_idx_lower][x_idx_lower+1]);
      a0_k_high = linear_interpolation(x, x_high[x_idx_lower], table_high[k_idx_lower+1][x_idx_lower], x_high[x_idx_lower+1], table_high[k_idx_lower+1][x_idx_lower+1]);
      a0 = linear_interpolation(k, ks[k_idx_lower], a0_k_low, ks[k_idx_lower+1], a0_k_high);

    }
    else {
      x_idx_lower = int( x / low_interval );
      // assert(x_idx_lower < num_low); // FIXME: if x==dividing_line -> x_lower==num_low -> index out of range
      a0_k_low = linear_interpolation(x, x_low[x_idx_lower], table_low[k_idx_lower][x_idx_lower], x_low[x_idx_lower+1], table_low[k_idx_lower][x_idx_lower+1]);
      a0_k_high = linear_interpolation(x, x_low[x_idx_lower], table_low[k_idx_lower+1][x_idx_lower], x_low[x_idx_lower+1], table_low[k_idx_lower+1][x_idx_lower+1]);
      a0 = linear_interpolation(k, ks[k_idx_lower], a0_k_low, ks[k_idx_lower+1], a0_k_high);

    }
    return a0;
  }

  double linear_interpolation(double x, double x_lower, double y_lower, double x_upper, double y_upper) {
    return ( (x - x_lower) * y_upper +  (x_upper - x) * y_lower ) / (x_upper - x_lower);
  }

};
