/////////////////////////////////////////////////////////////////
/// @file mdp_parameter.h
/// @version 2009-12-21
/// @author 2006 Piotr Bialas, 2009 onwards Lukasz Daniel
///
/// Contains class mdp_parameter
///
/// Distributed under GPL3 License
/// Read attached license in file mdp_license_v3.txt
/// This file cannot be distributed without file mdp_license_v3.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PARAMETER_
#define MDP_PARAMETER_

#include <map>
#include <string>
#include <variant>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace MDP
{
  class mdp_parameter
  {
  public:
    explicit mdp_parameter(const std::string &file)
    {
      std::ifstream fin(file);
      if (!fin)
        throw std::ios_base::failure("Cannot open file: " + file);

      read(fin);
    }

    explicit mdp_parameter(std::istream &in)
    {
      read(in);
    }

    /** @brief Check if parameter is defined
     *
     * @param name Parameter name
     * @return True if parameter is defined regardless of its type, false otherwise
     */
    bool defined(const std::string &name) const
    {
#if __cplusplus >= 202002L
      return m_map.contains(name);
#else
      return m_map.find(name) != m_map.end();
#endif
    }

    /** @brief Check if parameter is defined as double
     *
     * This is a convenience wrapper around defined_as<double>.
     *
     * @param name Parameter name
     * @return True if parameter is defined as double, false otherwise
     */
    inline bool defined_d(const std::string &name) const
    {
      return defined_as<double>(name);
    }

    /** @brief Check if parameter is defined as int
     *
     * This is a convenience wrapper around defined_as<int>.
     *
     * @param name Parameter name
     * @return True if parameter is defined as int, false otherwise
     */
    inline bool defined_i(const std::string &name) const
    {
      return defined_as<int>(name);
    }

    /** @brief Check if parameter is defined as string
     *
     * This is a convenience wrapper around defined_as<string>.
     *
     * @param name Parameter name
     * @return True if parameter is defined as string, false otherwise
     */
    inline bool defined_s(const std::string &name) const
    {
      return defined_as<std::string>(name);
    }

    template <typename T>
    bool defined_as(const std::string &name) const
    {
      auto it = m_map.find(name);
      if (it == m_map.end())
        return false;

      return std::holds_alternative<T>(it->second);
    }

    /** @brief Get parameter as double
     *
     * This is a convenience wrapper around get_as<double>.
     *
     * @param name Parameter name
     * @return Parameter value as double
     */
    inline double d(const std::string &name) const
    {
      return get_as<double>(name);
    }

    /** @brief Get parameter as int
     *
     * This is a convenience wrapper around get_as<int>.
     *
     * @param name Parameter name
     * @return Parameter value as int
     */
    inline int i(const std::string &name) const
    {
      return get_as<int>(name);
    }

    /** @brief Get parameter as string
     *
     * This is a convenience wrapper around get_as<string>.
     *
     * @param name Parameter name
     * @return Parameter value as string
     */
    inline std::string s(const std::string &name) const
    {
      return get_as<std::string>(name);
    }

    /** @brief Get parameter
     *
     * @param name Parameter name
     * @return Parameter value as variant<int, double, string>
     */
    const auto &get(const std::string &name) const
    {
      auto it = m_map.find(name);
      if (it == m_map.end())
        throw std::runtime_error("Parameter '" + name + "' not defined");

      return it->second;
    }

    /** @brief Get parameter
     *
     * @param name Parameter name
     * @return Parameter value as type T
     */
    template <typename T>
    T get_as(const std::string &name) const
    {
      auto it = m_map.find(name);
      if (it == m_map.end())
        throw std::runtime_error("Parameter '" + name + "' not defined");

      if (!std::holds_alternative<T>(it->second))
        throw std::runtime_error(
            "Parameter '" + name + "' has different type");

      return std::get<T>(it->second);
    }

    friend std::ostream &operator<<(std::ostream &os, const mdp_parameter &p)
    {
      for (const auto &[name, value] : p.m_map)
      {
        std::visit([&](const auto &v)
                   {
        using T = std::decay_t<decltype(v)>;

        if constexpr (std::is_same_v<T, int>)
          os << "int    ";
        else if constexpr (std::is_same_v<T, double>)
          os << "double ";
        else
          os << "string ";

        os << name << " = ";

        if constexpr (std::is_same_v<T, double>)
          os << std::setprecision(12) << v;
        else
          os << v;

        os << '\n'; }, value);
      }
      return os;
    }

  private:
    using ParamValue = std::variant<int, double, std::string>;
    std::map<std::string, ParamValue> m_map;

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
          m_map[name] = std::stoi(value);
        else if (type == "double")
          m_map[name] = std::stod(value);
        else if (type == "string")
          m_map[name] = value;
        else
          throw std::runtime_error("Unknown parameter type: " + type + " in line " + std::to_string(line_number));
      }
    }
  };
} // namespace MDP

#endif /* MDP_PARAMETER_ */
