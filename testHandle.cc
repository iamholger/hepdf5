#include <iostream>

#include "Handle.h"

using namespace hepdf5;

int main(int argc, char ** argv)
{
    auto handle = Handle(std::string(argv[1]));
    std::cerr << "num events: " << handle.nevents() << "\n";

    //for (auto s : handle.start_p(0,100)) std::cerr << s << "\n";
    //for (auto s : handle.getZ(0,100, "start_p")) std::cerr << s << "\n";
    //for (auto s : handle.get<int>(0, 100, std::string("start_p"))) std::cerr << s << "\n";
    //for (auto s : handle.get<double>(0, 100, std::string("vertex/x"))) std::cerr << s << "\n";

    readEvents(handle, 0, 100);
  return 0;
}
