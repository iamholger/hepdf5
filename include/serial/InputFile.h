#pragma once

#include <highfive/H5File.hpp>
using namespace HighFive;

namespace hepdf5
{

  class InputFile
  {
    public:
      InputFile(const string & fname) : file_(File(fname, File::ReadOnly)) {}
      
      DataSet getDataSet(const string & dsname) const {return file_.getDataSet(dsname);}
    
    private:
      HighFive::File file_;
  };
}

