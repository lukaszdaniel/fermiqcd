#ifndef MY_PARAMETERS_
#define MY_PARAMETERS_

#include <string>
#include <map>
#include <iostream>
#include <cstdio>
#include <cstring>

namespace MDP
{
  class parameter
  {
  private:
    std::map<std::string, double> _d_rep;
    std::map<std::string, int> _i_rep;
    std::map<std::string, std::string> _s_rep;

    int fscanf(FILE *fin)
    {
      constexpr int LENGTH = 64;
      constexpr int LINE_LENGTH = 256;
      constexpr char white_char[] = " \t\n";
      char line[LINE_LENGTH];

      char type[LENGTH];
      char name[LENGTH];
      char value[LENGTH];

      int line_number = 0;

      while (fgets(line, LINE_LENGTH, fin))
      {
        if (strspn(line, white_char) != strlen(line) && line[0] != '#')
        {
          ++line_number;

          if (sscanf(line, "%s %s %s\n", type, name, value) == 3)
          {
            if (strncmp(type, "int", LENGTH) == 0)
              _i_rep[std::string(name)] = atoi(value);
            else if (strncmp(type, "double", LENGTH) == 0)
              _d_rep[std::string(name)] = atof(value);
            else if (strncmp(type, "string", LENGTH) == 0)
              _s_rep[std::string(name)] = std::string(value);
            else
            {
              std::cerr << "unknown type " << type << std::endl;
              exit(1);
            }
          }
          else
          {
            std::cerr << "error reading line " << line_number << std::endl;
            exit(1);
          }
        }
      }
      return 0;
    }

  public:
    parameter(FILE *fin)
    {
      fscanf(fin);
    }

    int fprint(FILE *fout)
    {
      for (const auto &p : _i_rep)
      {
        fprintf(fout, "int  %s =  %d\n", (p.first).c_str(), p.second);
      }

      for (const auto &p : _d_rep)
      {
        fprintf(fout, "double  %s =  %.12g\n", (p.first).c_str(), p.second);
      }

      for (const auto &p : _s_rep)
      {
        fprintf(fout, "string  %s  =  %s\n", (p.first).c_str(), (p.second).c_str());
      }

      fflush(fout);
      return 0;
    }

    parameter(const char *file)
    {
      FILE *fin;

      if ((fin = fopen(file, "r")) == NULL)
      {
        std::cerr << "cannot open file " << file << std::endl;
        exit(1);
      }
      fscanf(fin);
      fclose(fin);
    }

    double d(std::string name)
    {
      const auto &p = _d_rep.find(name);

      if (p == _d_rep.end())
      {
        std::cerr << "double parameter " << name << " was not defined" << std::endl;
        exit(1);
      }
      return p->second;
    }

    int i(std::string name)
    {
      const auto &p = _i_rep.find(name);

      if (p == _i_rep.end())
      {
        std::cerr << "int parameter " << name << " was not defined" << std::endl;
        exit(1);
      }
      return p->second;
    }

    std::string s(std::string name)
    {
      const auto &p = _s_rep.find(name);

      if (p == _s_rep.end())
      {
        std::cerr << "string parameter " << name << " was not defined" << std::endl;
        exit(1);
      }
      return p->second;
    }

    bool defined_d(std::string name)
    {
      const auto &p = _d_rep.find(name);
      return (p != _d_rep.end());
    }

    bool defined_i(std::string name)
    {
      const auto &p = _i_rep.find(name);
      return (p != _i_rep.end());
    }

    bool defined_s(std::string name)
    {
      const auto &p = _s_rep.find(name);
      return (p != _s_rep.end());
    }
  };
} // namespace MDP

#endif /* MY_PARAMETERS_ */
