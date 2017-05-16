#ifndef MISC_HPP
#define MISC_HPP
#include <string>
#include <vector>

bool IsFileExist(std::string fname);
std::string DecommentFile(std::string fname);
void SplitString(std::string str, std::vector<std::string>& result);

#endif
