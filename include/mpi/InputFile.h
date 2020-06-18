#pragma once

#include "mpi.h"
#include <highfive/H5File.hpp>
using namespace HighFive;

namespace hepdf5
{

  class InputFile
  {
    public:
      InputFile(const string & fname) : file_(File(fname, File::ReadOnly, HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL))) {}
      
      DataSet getDataSet(const string & dsname) const {return file_.getDataSet(dsname);}
    
    private:
      HighFive::File file_;
  };
}
