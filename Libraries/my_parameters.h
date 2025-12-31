#ifndef MY_PARAMETERS_
#define MY_PARAMETERS_

#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace MDP
{
  class Parameter
  {
  public:
    explicit Parameter(const std::string &file)
    {
      std::ifstream fin(file);
      if (!fin)
        throw std::ios_base::failure("Cannot open file: " + file);

      read(fin);
    }

    explicit Parameter(std::istream &in)
    {
      read(in);
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
#if __cplusplus >= 202002L
      return _d_rep.contains(name);
#else
      return _d_rep.find(name) != _d_rep.end();
#endif
    }

    bool defined_i(const std::string &name) const
    {
#if __cplusplus >= 202002L
      return _i_rep.contains(name);
#else
      return _i_rep.find(name) != _i_rep.end();
#endif
    }

    bool defined_s(const std::string &name) const
    {
#if __cplusplus >= 202002L
      return _s_rep.contains(name);
#else
      return _s_rep.find(name) != _s_rep.end();
#endif
    }

    friend std::ostream &operator<<(std::ostream &os, const Parameter &p)
    {
      for (const auto &[name, value] : p._i_rep)
        os << "int    " << name << " =  " << value << '\n';

      for (const auto &[name, value] : p._d_rep)
        os << "double " << name << " =  " << std::setprecision(12) << value << '\n';

      for (const auto &[name, value] : p._s_rep)
        os << "string " << name << " =  " << value << '\n';

      return os;
    }

  private:
    std::map<std::string, double> _d_rep;
    std::map<std::string, int> _i_rep;
    std::map<std::string, std::string> _s_rep;

    void read(std::istream &in)
    {
      std::string line;
      std::size_t line_number = 0;

      while (std::getline(in, line))
      {
        ++line_number;

        auto pos = line.find_first_not_of(" \t");
        if (pos == std::string::npos || line[pos] == '#')
          continue;

        std::istringstream iss(line);
        std::string type, name, value;

        if (!(iss >> type >> name >> value))
          throw std::runtime_error("Syntax error in line " + std::to_string(line_number));

        if (type == "int")
          _i_rep[name] = std::stoi(value);
        else if (type == "double")
          _d_rep[name] = std::stod(value);
        else if (type == "string")
          _s_rep[name] = value;
        else
          throw std::runtime_error("Unknown parameter type: " + type);
      }
    }

    template <typename T>
    static T get_param(const std::map<std::string, T> &rep,
                       const std::string &name,
                       const char *type)
    {
      auto it = rep.find(name);
      if (it == rep.end())
        throw std::runtime_error(std::string(type) + " parameter '" + name + "' not defined");
      return it->second;
    }
  };
} // namespace MDP

#endif /* MY_PARAMETERS_ */
