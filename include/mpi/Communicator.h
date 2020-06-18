#pragma once

#include "mpi.h"

namespace hepdf5
{


  inline int getSize()
  {
    // TODO do I want to MPI_Init here???
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }

  inline int getRank()
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  class Communicator
  {
    public:
      Communicator() :  size_(getSize()), rank_(getRank()) {}
      
      int size() const {return size_;}
      int rank() const {return rank_;}
    
    private:
      int const size_;
      int const rank_;
  };
}
