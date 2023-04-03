#include "mdp_all.h"
#include <string>
#include <map>
#include <vector>


class Step
{
public:
  std::string algorithm;
  std::map<std::string, float> parameters;
};

void parse(int argc, char **argv)
{
  std::vector<Step> steps;
  Step step;
  std::string s;
  int i, j, k;

  for (i = 1; i < argc; i++)
  {
    std::string s = std::string(argv[i]);
    if (argv[i][0] == '-')
    {
      j = s.find(":");
      if (j < 0)
      {
      }
    }

    std::cout << argv[i] << std::endl;
  }
}

int main(int argc, char **argv)
{
  parse(argc, argv);
  return 0;
}
