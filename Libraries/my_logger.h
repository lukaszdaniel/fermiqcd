#ifndef MY_LOGGER_
#define MY_LOGGER_

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cassert>

namespace MDP
{
  class Logger
  {
  public:
    explicit Logger(const std::string &target,
                    int argc,
                    char *argv[],
                    const std::string &id)
    {
      open(target);
      header(argc, argv, id);
    }

    Logger(const Logger &) = delete;
    Logger &operator=(const Logger &) = delete;
    Logger(Logger &&) = default;
    Logger &operator=(Logger &&) = default;

    ~Logger() noexcept
    {
      try
      {
        footer();
        close();
      }
      catch (...)
      {
        // destructors must not throw
      }
    }

    template <typename T>
    Logger &operator<<(const T &value)
    {
      assert(out_);
      (*out_) << value;
      return *this;
    }

    Logger &flush()
    {
      assert(out_);
      out_->flush();
      return *this;
    }

  private:
    std::ostream *out_ = nullptr;
    std::ofstream file_;

    std::chrono::system_clock::time_point time_start_;
    std::clock_t cpu_start_{};

    bool cpu_ok_ = false;
    bool footer_done_ = false;

    static std::string format_time(const std::chrono::system_clock::time_point &tp)
    {
      std::time_t t = std::chrono::system_clock::to_time_t(tp);
      std::tm tm{};

      if (auto *ptm = std::localtime(&t))
        tm = *ptm;
      else
        return "time not available";

      std::ostringstream oss;
      oss << std::put_time(&tm, "%d.%m.%y %H:%M:%S");
      return oss.str();
    }

    void open(const std::string &target)
    {
      if (target == "stderr")
      {
        out_ = &std::cerr;
      }
      else if (target == "stdout")
      {
        out_ = &std::cout;
      }
      else
      {
        file_.open(target);
        if (!file_)
        {
          std::cerr << "cannot open log file " << target << "\n";
          std::abort();
        }
        out_ = &file_;
      }
    }

    void close()
    {
      if (file_.is_open())
        file_.close();
    }

    void header(int argc, char *argv[], const std::string &id)
    {
      assert(out_);

      time_start_ = std::chrono::system_clock::now();
      cpu_start_ = std::clock();
      cpu_ok_ = (cpu_start_ != static_cast<std::clock_t>(-1));

      (*out_) << "# program    " << argv[0] << "   " << id << "\n";
      (*out_) << "# time start " << format_time(time_start_) << "\n";

      for (int i = 0; i < argc; ++i)
        (*out_) << "# argv " << std::setw(2) << i << " " << argv[i] << "\n";

      (*out_) << "# header end\n";
      out_->flush();
    }

    void footer()
    {
      if (!out_ || footer_done_)
        return;

      footer_done_ = true;

      const auto time_finish = std::chrono::system_clock::now();
      const auto cpu_finish = std::clock();

      if (cpu_ok_)
      {
        const double cpu_time = static_cast<double>(cpu_finish - cpu_start_) / CLOCKS_PER_SEC;
        (*out_) << "# CPU time  " << std::setw(12) << cpu_time << " sec\n";
      }

      const double wall_time = std::chrono::duration<double>(time_finish - time_start_).count();
      (*out_) << "# user time " << std::setw(12) << wall_time << " sec\n";
      (*out_) << "# time stop " << format_time(time_finish) << "\n";
      out_->flush();
    }
  };
} // namespace MDP

#endif /* MY_LOGGER_ */
