#include <iostream>

#include "HepMC3/WriterAscii.h"
#include "Handle.h"

using namespace hepdf5;

int main(int argc, char ** argv)
{
    auto handle = Handle(std::string(argv[1]));
    //std::cerr << "num events: " << handle.nevents() << "\n";

    size_t numread=0;
    size_t nread=atoi(argv[2]);

  std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();;
    //WriterAscii writer("blabla.hepmc3", run);
  ////for (auto g : genevents) Print::listing(g);
  //for (auto g : genevents) writer.write_event(g);
    GenEvent _evt(Units::GEV,Units::MM);//, run);
    while (numread<handle.nevents())
    {
      if (handle.nevents() - numread < nread) nread = handle.nevents() - numread;
      handle.fillBuffers(numread, nread);

      for (size_t i=0; i<nread;++i)
      {
        handle.fillEvent(i, _evt);
        //writer.write_event(_evt);
      }
      //if (numread==0) 
        //Print::listing(_evt);
      numread+=nread;
    }

    return 0;
}
