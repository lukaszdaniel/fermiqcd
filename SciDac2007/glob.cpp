#include <functional>
#include <filesystem>
#include <string>
#include <iostream>
#include <vector>
#ifndef _WIN32
#include <glob.h>
#endif

namespace fs = std::filesystem;

std::vector<std::string> glob(const std::string &pattern)
{
  std::vector<std::string> v;

#ifndef _WIN32
  glob_t pglob;
  pglob.gl_offs = 2;

  if (glob(pattern.c_str(), 0, nullptr, &pglob) != 0)
  {
    v.push_back("?");
  }
  else
  {
    for (size_t i = 0; i < pglob.gl_pathc; ++i)
      v.emplace_back(pglob.gl_pathv[i]);
  }

  globfree(&pglob);

#else
  fs::path p(pattern);
  fs::path dir = p.parent_path();
  std::string mask = p.filename().string();

  if (dir.empty())
    dir = ".";

  if (!fs::exists(dir) || !fs::is_directory(dir))
  {
    v.push_back("?");
    return v;
  }

  std::function<bool(const char *, const char *)> match = [&](const char *str, const char *pat) -> bool
  {
    while (*pat)
    {
      if (*pat == '*')
      {
        ++pat;
        if (!*pat)
          return true;
        while (*str)
          if (match(str++, pat))
            return true;
        return false;
      }
      else if (*pat == '?')
      {
        if (!*str)
          return false;
        ++str;
        ++pat;
      }
      else
      {
        if (*str != *pat)
          return false;
        ++str;
        ++pat;
      }
    }
    return *str == '\0';
  };

  bool found = false;
  for (const auto &entry : fs::directory_iterator(dir))
  {
    const std::string name = entry.path().filename().string();
    if (match(name.c_str(), mask.c_str()))
    {
      v.push_back(entry.path().string());
      found = true;
    }
  }

  if (!found)
    v.push_back("?");

#endif

  return v;
}

int main()
{
  std::vector<std::string> v = glob("*.vtk");
  for (size_t i = 0; i < v.size(); i++)
    std::cout << v[i] << std::endl;
  return 0;
}
