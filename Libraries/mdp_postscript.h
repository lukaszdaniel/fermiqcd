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

#include <cstdio>

namespace MDP
{
  /// @brief to output and draw in postscript
  ///
  /// Example:
  /// @verbatim
  ///    mdp_postscript ps("test.ps");
  ///    ps.color(0.2,0.2,0.7);
  ///    ps.line(0,0,  5,5);
  ///    ps.font("Times-Roman", 12);
  ///    ps.print(5,5,"a line from (0,0) to here");
  /// @endverbatim
  class mdp_postscript
  {
  private:
    FILE *m_fp;
    float m_red;
    float m_green;
    float m_blue;
    float m_scale;

    FILE *open(const std::string &filename)
    {
      std::cout << "Making frame: " << filename << std::endl;
      m_fp = fopen(filename.c_str(), "w");
      if (m_fp)
        fprintf(m_fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
      fflush(m_fp);
      return m_fp;
    }

    void close()
    {
      fprintf(m_fp, "showpage\n");
      fprintf(m_fp, "%%%%Trailer\n");
      fprintf(m_fp, "%%%%EOF\n");
      fflush(m_fp);
      fclose(m_fp);
      m_fp = nullptr;
    }

  public:
    mdp_postscript() : m_fp(nullptr), m_red(0), m_green(0), m_blue(0), m_scale(1)
    {
    }

    mdp_postscript(const std::string &filename) : m_fp(nullptr), m_red(0), m_green(0), m_blue(0), m_scale(1)
    {
      open(filename);
    }

    virtual ~mdp_postscript()
    {
      if (m_fp)
        close();
    }

    void size(float x0, float y0, float x1, float y1)
    {
      fprintf(m_fp, "%%%%BoundingBox: %.0f %.0f %.0f %.0f\n", m_scale * x0, m_scale * y0, m_scale * x1, m_scale * y1);
      fprintf(m_fp, "%%%%EndComments\n");
      fflush(m_fp);
    }

    void line(float x0, float y0, float x1, float y1)
    {
      fprintf(m_fp, "%.2f %.2f %.2f setrgbcolor\n", m_red, m_green, m_blue);
      fprintf(m_fp, "%.2f %.2f moveto\n", m_scale * x0, m_scale * y0);
      fprintf(m_fp, "%.2f %.2f lineto\n", m_scale * x1, m_scale * y1);
      fprintf(m_fp, "stroke\n");
      fflush(m_fp);
    }

    void box(float x0, float y0, float x1, float y1, int fill = 0)
    {
      if (fill == 1)
        fprintf(m_fp, "stroke\n");
      fprintf(m_fp, "%.2f %.2f moveto\n", m_scale * x0, m_scale * y0);
      fprintf(m_fp, "%.2f %.2f lineto\n", m_scale * x0, m_scale * y1);
      fprintf(m_fp, "%.2f %.2f lineto\n", m_scale * x1, m_scale * y1);
      fprintf(m_fp, "%.2f %.2f lineto\n", m_scale * x1, m_scale * y0);
      fprintf(m_fp, "%.2f %.2f lineto\n", m_scale * x0, m_scale * y0);
      if (fill == 1)
        fprintf(m_fp, "fill\n");
      fprintf(m_fp, "stroke\n");
      fflush(m_fp);
    }

    void arc(float x0, float y0, float r, float alpha, float beta)
    {
      fprintf(m_fp, "%.2f %.2f %.2f %.2f %.2f arc\n", m_scale * x0, m_scale * y0, m_scale * r, alpha, beta);
      fprintf(m_fp, "stroke\n");
      fflush(m_fp);
    }

    void circle(float x0, float y0, float r, int fill = 0)
    {
      constexpr int BOLD = 10;
      if (fill != BOLD)
        fprintf(m_fp, "%.2f %.2f %.2f 0 360 arc\n", m_scale * x0, m_scale * y0, r);
      if (fill == 1)
        fprintf(m_fp, "fill\n");
      if (fill == BOLD)
      {
        fprintf(m_fp, "1 -0.1 0.1 {\n");
        fprintf(m_fp, "/r exch def\n");
        fprintf(m_fp, "%.2f 1 r 0.5 mul sub mul ", m_red);
        fprintf(m_fp, "%.2f 1 r 0.5 mul sub mul ", m_green);
        fprintf(m_fp, "%.2f 1 r 0.5 mul sub mul ", m_blue);
        fprintf(m_fp, "setrgbcolor\n");
        fprintf(m_fp, "stroke\n");
        fprintf(m_fp, "%.2f %.2f %.2f r mul 0 360 arc fill\n", m_scale * x0, m_scale * y0, m_scale * r);
        fprintf(m_fp, "} for\n");
      };
      fflush(m_fp);
    }

    void pen(float size)
    {
      fprintf(m_fp, "%.2f setlinewidth\n", size);
    }

    // colors are numbers in [0,1] black=(0,0,0) white=(1,1,1)
    void color(float r, float g, float b)
    {
      fprintf(m_fp, "%.2f %.2f %.2f setrgbcolor\n", r, g, b);
      m_red = r;
      m_green = g;
      m_blue = b;
      fflush(m_fp);
    }

    void font(const char *text, int size)
    {
      fprintf(m_fp, "/%s findfont\n%i scalefont\nsetfont\n", text, size);
      fflush(m_fp);
    }

    void print(float x0, float y0, const char text[])
    {
      fprintf(m_fp, "%.2f %.2f moveto\n", m_scale * x0, m_scale * y0);
      fprintf(m_fp, "(%s) show\n", text);
      fprintf(m_fp, "stroke\n");
      fflush(m_fp);
    }
  };
} // namespace MDP

#endif /* MDP_POSTSCRIPT_ */
