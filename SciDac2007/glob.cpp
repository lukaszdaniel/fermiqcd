#ifndef _WIN64
#include "glob.h"
#endif
#include <string>
#include <iostream>
#include <vector>


std::vector<std::string> glob(std::string pattern)
{
  std::vector<std::string> v;
  glob_t pglob;
  pglob.gl_offs = 2;
  if (glob(pattern.c_str(), 0, 0, &pglob) != 0)
    v.push_back("?");
  else
    for (unsigned int i = 0; i < pglob.gl_pathc; i++)
      v.push_back(std::string(pglob.gl_pathv[i]));
  globfree(&pglob);
  return v;
}

int main()
{
  std::vector<std::string> v = glob("*.vtk");
  for (int i = 0; i < v.size(); i++)
    std::cout << v[i] << std::endl;
  return 0;
}
