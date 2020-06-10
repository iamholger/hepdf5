#pragma once

#include <unordered_map>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
using namespace HighFive;

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenParticle.h"
using namespace HepMC3;

namespace hepdf5
{
    template<typename T>
    std::vector<T> read_vector(size_t offset, size_t readsize, DataSet const & ds)
    {
        std::vector<T> vec;
        vec.reserve(readsize);
        ds.select({offset}, {readsize}).read(vec);
        return vec;
    }
    

    class Handle
    {
        public:
            explicit Handle(
                    std::string const & fname, 
                    std::vector<size_t> const & cols = {}
                    ) : cols_(cols), file_(File(fname, File::ReadOnly)),
                P_START      ( file_.getDataSet("event/start_p")   ),
                V_START      ( file_.getDataSet("event/start_v")   ),
                P_N          ( file_.getDataSet("event/npart")     ),
                V_N          ( file_.getDataSet("event/nvert")     ),
                WGT          ( file_.getDataSet("event/weight")    ),
                V_ID         ( file_.getDataSet("vertex/id")       ),
                V_STATUS     ( file_.getDataSet("vertex/status")   ),
                V_X          ( file_.getDataSet("vertex/x")        ),
                V_Y          ( file_.getDataSet("vertex/y")        ),
                V_Z          ( file_.getDataSet("vertex/z")        ),
                V_T          ( file_.getDataSet("vertex/t")        ),
                P_ID         ( file_.getDataSet("particle/pid")    ),
                P_STATUS     ( file_.getDataSet("particle/status") ),
                P_PX         ( file_.getDataSet("particle/px")     ),
                P_PY         ( file_.getDataSet("particle/py")     ),
                P_PZ         ( file_.getDataSet("particle/pz")     ),
                P_E          ( file_.getDataSet("particle/E")      ),
                P_M          ( file_.getDataSet("particle/m")      ),
                P_VTX_END    ( file_.getDataSet("particle/vtx_end")),
                P_VTX_PROD   ( file_.getDataSet("particle/vtx_prod")) {};

            long int nevents() {return H5Sget_simple_extent_npoints(H5Dget_space(P_N.getId()));};
            size_t nweights() {return H5Sget_simple_extent_npoints(H5Dget_space(file_.getDataSet("weight_names").getId()));};
            

            template<typename T>
            std::vector<T> get(size_t offset, size_t readsize, std::string const & dsname)
            {
                if      (dsname == "event/npart"      ) return read_vector<T>(offset, readsize, this->P_N);
                else if (dsname == "event/nvert"      ) return read_vector<T>(offset, readsize, this->V_N);
                else if (dsname == "event/start_p"    ) return read_vector<T>(offset, readsize, this->P_START);
                else if (dsname == "event/start_v"    ) return read_vector<T>(offset, readsize, this->V_START);
                else if (dsname == "vertex/x"         ) return read_vector<T>(offset, readsize, this->V_X);
                else if (dsname == "vertex/y"         ) return read_vector<T>(offset, readsize, this->V_Y);
                else if (dsname == "vertex/z"         ) return read_vector<T>(offset, readsize, this->V_Z);
                else if (dsname == "vertex/t"         ) return read_vector<T>(offset, readsize, this->V_T);
                else if (dsname == "vertex/id"        ) return read_vector<T>(offset, readsize, this->V_ID);
                else if (dsname == "vertex/status"    ) return read_vector<T>(offset, readsize, this->V_STATUS);
                else if (dsname == "particle/pid"     ) return read_vector<T>(offset, readsize, this->P_ID      );
                else if (dsname == "particle/status"  ) return read_vector<T>(offset, readsize, this->P_STATUS  );
                else if (dsname == "particle/px"      ) return read_vector<T>(offset, readsize, this->P_PX      );
                else if (dsname == "particle/py"      ) return read_vector<T>(offset, readsize, this->P_PY      );
                else if (dsname == "particle/pz"      ) return read_vector<T>(offset, readsize, this->P_PZ      );
                else if (dsname == "particle/E"       ) return read_vector<T>(offset, readsize, this->P_E       );
                else if (dsname == "particle/m"       ) return read_vector<T>(offset, readsize, this->P_M       );
                else if (dsname == "particle/vtx_end" ) return read_vector<T>(offset, readsize, this->P_VTX_END );
                else if (dsname == "particle/vtx_prod") return read_vector<T>(offset, readsize, this->P_VTX_PROD);
                else 
                {
                    std::cerr<< dsname << "\n";
                    throw 20;
                }
            }
    
            std::vector<std::vector<double> > weights(size_t offset, size_t readsize)//; // TODO augment with columns
            {
                vector<vector<double> > e_weights;
                e_weights.reserve(readsize);
                this->WGT.select({offset, 0}, {readsize, this->nweights()}).read(e_weights);
                return e_weights;
            }

        private:
            std::vector<size_t> const cols_;
            HighFive::File const file_;
            DataSet const P_START;
            DataSet const V_START;
            DataSet const P_N;
            DataSet const V_N;
            DataSet const WGT;
            DataSet const V_ID;
            DataSet const V_STATUS;
            DataSet const V_X;
            DataSet const V_Y;
            DataSet const V_Z;
            DataSet const V_T;
            DataSet const P_ID;
            DataSet const P_STATUS;
            DataSet const P_PX; 
            DataSet const P_PY;
            DataSet const P_PZ;
            DataSet const P_E;
            DataSet const P_M;
            DataSet const P_VTX_END;
            DataSet const P_VTX_PROD;
    };
    
    std::vector<GenEvent> readEvents(Handle & handle,  size_t first_event, size_t n_events);
}
