#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "misc.hpp"

bool IsFileExist(std::string fname)
{
  std::ifstream ifile(fname.c_str());
  return ifile.is_open();
}

std::string DecommentFile(std::string fname)
{
  std::stringstream msg;
  if (!IsFileExist(fname)) {
    msg << "### FATAL ERROR in DecommentFile. File doesn't exist." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  std::ifstream file(fname.c_str(), std::ios::in);
  std::string ss;
  char c;
  while (file) {
    file.get(c);
    if (c == '#') {
      while (c != '\n' && file) file.get(c);
      continue;
    }
    ss += c;
  }
  return ss;
}

void SplitString(std::string str, std::vector<std::string>& result)
{
  std::istringstream ss(str);
  result.clear();

  do {
    std::string sub;
    ss >> sub;
    result.push_back(sub);
  } while (ss);
}
