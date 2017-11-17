#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

enum LogEntry {
  LOGENTRY_DENSITY = 0,
  LOGENTRY_VELOCITY,
  LOGENTRY_PRESSURE,
  LOGENTRY_NFRAC,
  NUMBER_OF_LOGENTRIES
};

struct Entry {
  unsigned long _previous_entry;
  unsigned short _index;
  LogEntry _entry;
  double _time;
  double _value;
};

Entry read_entry(unsigned long position, std::ifstream &file) {
  file.seekg(position);
  Entry entry;
  file.read(reinterpret_cast<char *>(&entry._previous_entry),
            sizeof(unsigned long));
  file.read(reinterpret_cast<char *>(&entry._index), sizeof(unsigned short));
  file.read(reinterpret_cast<char *>(&entry._entry), sizeof(int));
  file.read(reinterpret_cast<char *>(&entry._time), sizeof(double));
  file.read(reinterpret_cast<char *>(&entry._value), sizeof(double));
  return entry;
}

void print_entry(const Entry &entry, std::ostream &out) {
  out << "Entry:\n";
  out << entry._previous_entry << "\n";
  out << entry._index << "\n";
  out << entry._entry << "\n";
  out << entry._time << "\n";
  out << entry._value << "\n";
}

std::map<double, double> get_quantity(LogEntry type, unsigned short index,
                                      std::ifstream &ifile) {
  const unsigned int lastpos = ifile.tellg();
  // find first entry for particle
  unsigned long pos = ifile.tellg();
  pos -= 30;
  Entry entry = read_entry(pos, ifile);
  while (entry._index != index) {
    pos -= 30;
    entry = read_entry(pos, ifile);
  }
  // bingo: found last entry
  // now try to find the last entry containing the quantity of interest
  // (this is guaranteed to exist)
  long next_pos = pos - entry._previous_entry;
  while (entry._entry != type) {
    pos = next_pos;
    entry = read_entry(pos, ifile);
    next_pos = pos - entry._previous_entry;
  }
  // bingo: found last entry of interest
  std::map<double, double> table;
  while (next_pos >= 0 && pos != next_pos && entry._index == index) {
    if (entry._entry == type) {
      table[entry._time] = entry._value;
      //      print_entry(entry, std::cout);
    }
    pos = next_pos;
    entry = read_entry(pos, ifile);
    next_pos = pos - entry._previous_entry;
  }
  if (entry._index == index && entry._entry == type) {
    table[entry._time] = entry._value;
    //    print_entry(entry, std::cout);
  }

  // reset file
  ifile.seekg(lastpos);
  return table;
}

int main(int argc, char **argv) {

  std::ifstream ifile("logfile.dat", std::ios::binary | std::ios::ate);
  const unsigned long lastpos = ifile.tellg();
  Entry entry = read_entry(lastpos - 30, ifile);
  const double endtime = entry._time;

  std::vector<std::vector<double>> snaps(10, std::vector<double>(1000, 0.));

  const double snaptime = 0.1 * endtime;
  for (unsigned int i = 0; i < 1000; ++i) {
    std::map<double, double> table = get_quantity(LOGENTRY_DENSITY, i, ifile);
    for (unsigned int j = 0; j < 10; ++j) {
      // find the time keys surrounding the snapshot time
      const double t = j * snaptime;
      auto upper_bound = table.upper_bound(t);
      auto lower_bound = upper_bound--;
      const double t0 = lower_bound->first;
      const double t1 = upper_bound->first;
      const double v0 = lower_bound->second;
      const double v1 = upper_bound->second;
      const double tfac = (t - t0) / (t1 - t0);
      snaps[j][i] = (1. - tfac) * v0 + tfac * v1;
    }
  }

  for (unsigned int j = 0; j < 10; ++j) {
    std::stringstream filename;
    filename << "logsnap_";
    filename.fill('0');
    filename.width(4);
    filename << j;
    filename << ".txt";
    std::ofstream ofile(filename.str());
    for (unsigned int i = 0; i < 1000; ++i) {
      ofile << i * 0.001 << "\t" << snaps[j][i] << "\n";
    }
  }

  return 0;
}
