#pragma once

namespace hepdf5
{
  class Communicator
  {
    public:
      int size() {return 1;}
      int rank() {return 0;}
  };
}
