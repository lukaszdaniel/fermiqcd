/////////////////////////////////////////////////////////////////
/// @file mdp_postscript.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Yes...MDP can print and draw in postscript
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_POSTSCRIPT_
#define MDP_POSTSCRIPT_

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

namespace MDP
{
  /// @brief to output and draw in postscript
  ///
  /// Example:
  /// @verbatim
  ///    mdp_postscript ps("test.eps");
  ///    ps.color(0.2,0.2,0.7);
  ///    ps.line(0,0,  5,5);
  ///    ps.font("Times-Roman", 12);
  ///    ps.print(5,5,"a line from (0,0) to here");
  /// @endverbatim
  class mdp_postscript
  {
  private:
    std::ofstream m_fp;
    float m_red;
    float m_green;
    float m_blue;
    float m_scale;

    void open(const std::string &filename)
    {
      std::cout << "Making frame: " << filename << std::endl;
      m_fp.open(filename, std::ios::out);
      if (m_fp.is_open())
        m_fp << "%!PS-Adobe-3.0 EPSF-3.0\n";
      m_fp.flush();
    }

    void close()
    {
      m_fp << "showpage\n";
      m_fp << "%%%%Trailer\n";
      m_fp << "%%%%EOF\n";
      m_fp.flush();
      m_fp.close();
    }

  public:
    mdp_postscript() : m_red(0), m_green(0), m_blue(0), m_scale(1)
    {
    }

    explicit mdp_postscript(const std::string &filename) : m_red(0), m_green(0), m_blue(0), m_scale(1)
    {
      open(filename);
    }

    ~mdp_postscript()
    {
      if (m_fp.is_open())
        close();
    }

    /** @brief Set paper size
     */
    void size(float x0, float y0, float x1, float y1)
    {
      m_fp << "%%%%BoundingBox: "
            << static_cast<int>(m_scale * x0) << " "
            << static_cast<int>(m_scale * y0) << " "
            << static_cast<int>(m_scale * x1) << " "
            << static_cast<int>(m_scale * y1) << "\n";
      m_fp << "%%%%EndComments\n";
      m_fp.flush();
    }

    /** @brief Draw a line segment
     */
    void line(float x0, float y0, float x1, float y1)
    {
      m_fp << std::fixed << std::setprecision(2)
            << m_red << " " << m_green << " " << m_blue << " setrgbcolor\n"
            << m_scale * x0 << " " << m_scale * y0 << " moveto\n"
            << m_scale * x1 << " " << m_scale * y1 << " lineto\n"
            << "stroke\n";
      m_fp.flush();
    }

    /** @brief Draw a box
     */
    void box(float x0, float y0, float x1, float y1, bool fill = false)
    {
      if (fill)
        m_fp << "stroke\n";

      m_fp << std::fixed << std::setprecision(2)
            << m_scale * x0 << " " << m_scale * y0 << " moveto\n"
            << m_scale * x0 << " " << m_scale * y1 << " lineto\n"
            << m_scale * x1 << " " << m_scale * y1 << " lineto\n"
            << m_scale * x1 << " " << m_scale * y0 << " lineto\n"
            << m_scale * x0 << " " << m_scale * y0 << " lineto\n";

      if (fill)
        m_fp << "fill\n";
      m_fp << "stroke\n";
      m_fp.flush();
    }

    /** @brief Draw an arc
     */
    void arc(float x0, float y0, float r, float alpha, float beta)
    {
      m_fp << std::fixed << std::setprecision(2)
            << m_scale * x0 << " " << m_scale * y0 << " "
            << m_scale * r << " " << alpha << " " << beta << " arc\n"
            << "stroke\n";
      m_fp.flush();
    }

    /** @brief Draw a circle
     */
    void circle(float x0, float y0, float r, int fill = 0)
    {
      constexpr int BOLD = 10;
      m_fp << std::fixed << std::setprecision(2);

      if (fill != BOLD)
        m_fp << m_scale * x0 << " " << m_scale * y0 << " " << r << " 0 360 arc\n";

      if (fill == 1)
        m_fp << "fill\n";

      if (fill == BOLD)
      {
        m_fp << "1 -0.1 0.1 {\n"
              << "/r exch def\n"
              << m_red << " 1 r 0.5 mul sub mul "
              << m_green << " 1 r 0.5 mul sub mul "
              << m_blue << " 1 r 0.5 mul sub mul setrgbcolor\n"
              << "stroke\n"
              << m_scale * x0 << " " << m_scale * y0 << " " << m_scale * r << " 0 360 arc fill\n"
              << "} for\n";
      }
      m_fp.flush();
    }

    /** @brief Set line width
     */
    void pen(float size)
    {
      m_fp << std::fixed << std::setprecision(2)
            << size << " setlinewidth\n";
      m_fp.flush();
    }

    /** @brief Set default colours
     *
     * @note colours are numbers in [0,1] black=(0,0,0) white=(1,1,1)
     */
    void color(float r, float g, float b)
    {
      m_red = r;
      m_green = g;
      m_blue = b;
      m_fp << std::fixed << std::setprecision(2)
            << r << " " << g << " " << b << " setrgbcolor\n";
      m_fp.flush();
    }

    /** @brief Set default font
     */
    void font(const std::string &text, int size)
    {
      m_fp << "/" << text << " findfont\n"
            << size << " scalefont\nsetfont\n";
      m_fp.flush();
    }

    /** @brief Print text
     */
    void print(float x0, float y0, const std::string &text)
    {
      m_fp << std::fixed << std::setprecision(2)
            << m_scale * x0 << " " << m_scale * y0 << " moveto\n"
            << "(" << text << ") show\n"
            << "stroke\n";
      m_fp.flush();
    }
  };
} // namespace MDP

#endif /* MDP_POSTSCRIPT_ */
