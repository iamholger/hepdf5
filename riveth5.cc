//#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
using namespace HighFive;

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenParticle.h"

#include "YODA/ReaderYODA.h"
#include "YODA/WriterYODA.h"
#include "YODA/AnalysisObject.h"

#include "Rivet/AnalysisHandler.hh"
//#include "Rivet/Log.hh"
using namespace HepMC3;


namespace hepdf5 {
   

    // Read function, returns an Events struct --- this is for the new structure
    void readEvents(Group& g_event, Group& g_vertex, Group& g_particle,  size_t first_event, size_t n_events)
    {
        DataSet P_START      = g_event.getDataSet("start_p");
        DataSet V_START      = g_event.getDataSet("start_v");
        DataSet P_N          = g_event.getDataSet("npart");
        DataSet V_N          = g_event.getDataSet("nvert");
        DataSet WGT          = g_event.getDataSet("weight");
        DataSet V_ID         = g_vertex.getDataSet("id");
        DataSet V_STATUS     = g_vertex.getDataSet("status");
        DataSet P_ID         = g_particle.getDataSet("pid");
        DataSet P_STATUS     = g_particle.getDataSet("status");
        DataSet P_PX         = g_particle.getDataSet("px");
        DataSet P_PY         = g_particle.getDataSet("py");
        DataSet P_PZ         = g_particle.getDataSet("pz");
        DataSet P_E          = g_particle.getDataSet("E");
        DataSet P_M          = g_particle.getDataSet("m");
        DataSet P_VTX_END    = g_particle.getDataSet("vtx_end");
        DataSet P_VTX_PROD   = g_particle.getDataSet("vtx_prod");
        
        std::vector<size_t> v_start, p_start;
        std::vector<int>    v_n, p_n, v_id, p_id, v_status, p_status, p_vtx_end, p_vtx_prod;
        std::vector<double> p_px, p_py, p_pz, p_e, p_m, wgt;
       
        std::vector<size_t> offset_e   = {first_event};
        std::vector<size_t> readsize_e = {n_events}   ;

        // Indexing info first
        P_START.select(offset_e, readsize_e).read(p_start);
        V_START.select(offset_e, readsize_e).read(v_start);
        p_n.reserve(n_events);
        v_n.reserve(n_events);
        P_N.select(offset_e, readsize_e).read(p_n);
        V_N.select(offset_e, readsize_e).read(v_n);
        wgt.reserve(n_events);
        WGT.select(offset_e, readsize_e).read(wgt);
        
        // First particle and vertex, some arithmatic to see how much reading is necessary
        std::vector<size_t> p_offset   = {p_start.front()};
        std::vector<size_t> v_offset   = {v_start.front()};
        
        size_t p_read = p_start.back() - p_start.front() + p_n.back();
        size_t v_read = v_start.back() - v_start.front() + v_n.back();
        std::vector<size_t> p_readsize = {p_read};
        std::vector<size_t> v_readsize = {v_read};

        p_id      .reserve(p_read);
        p_status  .reserve(p_read);
        p_px      .reserve(p_read);
        p_py      .reserve(p_read);
        p_pz      .reserve(p_read);
        p_e       .reserve(p_read);
        p_m       .reserve(p_read);
        p_vtx_end .reserve(p_read);
        p_vtx_prod.reserve(p_read);
        v_id      .reserve(v_read);
        v_status  .reserve(v_read);
        
        // Read dataset extends to vectors
        P_ID      .select(p_offset, p_readsize).read(p_id      );
        P_STATUS  .select(p_offset, p_readsize).read(p_status  );
        P_PX      .select(p_offset, p_readsize).read(p_px      );
        P_PY      .select(p_offset, p_readsize).read(p_py      );
        P_PZ      .select(p_offset, p_readsize).read(p_pz      );
        P_E       .select(p_offset, p_readsize).read(p_e       );
        P_M       .select(p_offset, p_readsize).read(p_m       );
        P_VTX_END .select(p_offset, p_readsize).read(p_vtx_end );
        P_VTX_PROD.select(p_offset, p_readsize).read(p_vtx_prod);
        V_ID      .select(v_offset, v_readsize).read(v_id      );
        V_STATUS  .select(v_offset, v_readsize).read(v_status  );



        //Rivet::Log::setLevel("Rivet.AnalysisHandler", Rivet::Log::TRACE);
  //double t0 = MPI_Wtime();
        size_t idxV(0), idxP(0);
        std::vector<GenVertexPtr> V;
        std::unordered_map<int,int> vmap;
        auto ah = new Rivet::AnalysisHandler;
        ah->addAnalysis("MC_GENERIC");
        ah->setIgnoreBeams(true);
        // Event loop
        for (size_t nev=0; nev<n_events; ++nev)
        {
          for (unsigned int iv=idxV; iv < (v_n[nev]+idxV); ++iv)
          {
            GenVertexPtr v = std::make_shared<GenVertex>();
            v->set_id(        v_id[iv]);
            v->set_status(v_status[iv]);
            V.push_back(v);
            vmap.insert({v_id[iv], iv-idxV});
          }
          
          for (unsigned int ip=idxP; ip < (p_n[nev]+idxP); ++ip)
          {
            GenParticlePtr p = std::make_shared<GenParticle>(
                FourVector(p_px[ip], p_py[ip], p_pz[ip], p_e[ip]),
                p_id[ip], p_status[ip]);
            //p->set_generated_mass(p_m[ip]);
            if (p_vtx_end[ip]  !=0) V[vmap[ p_vtx_end[ip]]]->add_particle_in(p);
            if (p_vtx_prod[ip] !=0) V[vmap[p_vtx_prod[ip]]]->add_particle_out(p);

          }

          idxV += v_n[nev];
          idxP += p_n[nev];
          
          GenEvent _evt(Units::GEV,Units::MM);
          _evt.set_event_number(nev);
          for (auto v: V) _evt.add_vertex(v);
          //_evt.weights().push_back(wgt[nev]);
          _evt.weights()= {wgt[nev]};
          V.clear();
          vmap.clear();
          //const GenEvent evt(_evt);

          //if (nev==0) {
            //ah->setWeightNames(_evt);
          //}
      
          try {ah->analyze( {_evt} ) ;} catch (const std::exception& e)
          {
          }



        }
        ah->setCrossSection(1* 1.0E9, 0.1 * 1e9 );
        ah->finalize();
        ah->writeData("bla.yoda");
        }
  //std::cout << "Processing loop took: " << t1-t0 << "\n";
}

int main(int argc, char ** argv)
{
  int rank, size;
  //MPI_Init(&argc, &argv);
  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  File file(argv[1], File::ReadOnly);//, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
  Group g_part = file.getGroup("particle");
  Group g_evt  = file.getGroup("event");
  Group g_vtx  = file.getGroup("vertex");
  hid_t dspace = H5Dget_space(file.getDataSet("event/npart").getId());
 
  unsigned int nread = 0;
  long int nEvents =  H5Sget_simple_extent_npoints(dspace);
  std::cout << "File contains " << nEvents << " events\n";
  //if (argc>=2) nread=atoi(argv[2]);
  //if (nread < nEvents) nEvents=nread;
  std::cout << "Total number of ranks: " << size << "\n";

  //double t0 = MPI_Wtime();
  hepdf5::readEvents(g_evt, g_vtx, g_part, 0, nEvents);
  //double t1 = MPI_Wtime();
  //std::cout << "Processing took: " << t1-t0 << "\n";

  return 0;
}
