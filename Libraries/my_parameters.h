#ifndef MY_PARAMETERS_
#define MY_PARAMETERS_

#include <string>
#include <map>
#include <iostream>
#include <cstdio>
#include <cstdlib>
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

          if (sscanf(line, "%s %s %s", type, name, value) == 3)
          {
            if (strncmp(type, "int", LENGTH) == 0)
            {
              _i_rep[std::string(name)] = std::atoi(value);
            }
            else if (strncmp(type, "double", LENGTH) == 0)
            {
              _d_rep[std::string(name)] = std::atof(value);
            }
            else if (strncmp(type, "string", LENGTH) == 0)
            {
              _s_rep[std::string(name)] = std::string(value);
            }
            else
            {
              std::cerr << "Unknown type " << type << std::endl;
              throw std::runtime_error("Unknown parameter type");
            }
          }
          else
          {
            std::cerr << "Error reading line " << line_number << std::endl;
            throw std::runtime_error("Error reading configuration file");
          }
        }
      }
      return 0;
    }

  public:
    explicit parameter(FILE *fin)
    {
      if (!fin)
      {
        std::cerr << "File pointer is null." << std::endl;
        throw std::invalid_argument("Invalid file pointer");
      }
      fscanf(fin);
    }

    explicit parameter(const char *file)
    {
      FILE *fin = fopen(file, "r");
      if (!fin)
      {
        std::cerr << "Cannot open file " << file << std::endl;
        throw std::ios_base::failure("File not found");
      }
      fscanf(fin);
      fclose(fin);
    }

    ~parameter() = default;

    int fprint(FILE *fout) const
    {
      if (!fout)
      {
        std::cerr << "File pointer is null." << std::endl;
        throw std::invalid_argument("Invalid file pointer");
      }

      for (const auto &p : _i_rep)
      {
        fprintf(fout, "int    %s =  %d\n", p.first.c_str(), p.second);
      }

      for (const auto &p : _d_rep)
      {
        fprintf(fout, "double %s =  %.12g\n", p.first.c_str(), p.second);
      }

      for (const auto &p : _s_rep)
      {
        fprintf(fout, "string %s =  %s\n", p.first.c_str(), p.second.c_str());
      }

      fflush(fout);
      return 0;
    }

    double d(const std::string &name) const
    {
      return get_param(_d_rep, name, "double");
    }

    int i(const std::string &name) const
    {
      return get_param(_i_rep, name, "int");
    }

    std::string s(const std::string &name) const
    {
      return get_param(_s_rep, name, "string");
    }

    bool defined_d(const std::string &name) const
    {
      return _d_rep.find(name) != _d_rep.end();
    }

    bool defined_i(const std::string &name) const
    {
      return _i_rep.find(name) != _i_rep.end();
    }

    bool defined_s(const std::string &name) const
    {
      return _s_rep.find(name) != _s_rep.end();
    }

  private:
    template <typename T>
    T get_param(const std::map<std::string, T> &rep, const std::string &name, const char *type) const
    {
      const auto &p = rep.find(name);
      if (p == rep.end())
      {
        std::cerr << type << " parameter " << name << " was not defined" << std::endl;
        throw std::runtime_error("Parameter not defined");
      }
      return p->second;
    }
  };
} // namespace MDP

#endif /* MY_PARAMETERS_ */
