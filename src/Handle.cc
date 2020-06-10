#include "Handle.h"

using namespace hepdf5;


//long int Handle::nevents()
//{
  //hid_t dspace = H5Dget_space(file_.P_N.getId());
  //return H5Sget_simple_extent_npoints(H5Dget_space(file_.P_N.getId()));
//}








vector<GenEvent> hepdf5::readEvents(Handle & handle,  size_t first_event, size_t n_events)
//void hepdf5::readEvents(Handle & handle,  size_t first_event, size_t n_events)
{
    // Weights
    //vector<vector<double> > e_weights;
    //e_weights.reserve(n_events);
    // NOTE this could be augmented to read only a central weight
    //WGT.select({first_event, 0}, {n_events, wnames.size()}).read(e_weights);

    // Event-wise quantities
    auto const p_start = handle.get<size_t>(first_event, n_events, "event/start_p");
    auto const v_start = handle.get<size_t>(first_event, n_events, "event/start_v");
    auto const p_n     = handle.get<int>   (first_event, n_events, "event/npart"  );
    auto const v_n     = handle.get<int>   (first_event, n_events, "event/nvert"  );

    // Some arithmetic to see how much reading is necessary for vertices
    size_t v_offset  = v_start.front();
    size_t v_readsize = v_start.back() - v_start.front() + v_n.back();
    //std::vector<size_t> v_readsize = {v_read};

    // Vertex quantities
    auto const v_id      = handle.get<int>    (v_offset, v_readsize, "vertex/id");
    auto const v_status  = handle.get<int>    (v_offset, v_readsize, "vertex/status");
    auto const v_x       = handle.get<double> (v_offset, v_readsize, "vertex/x");
    auto const v_y       = handle.get<double> (v_offset, v_readsize, "vertex/y");
    auto const v_z       = handle.get<double> (v_offset, v_readsize, "vertex/z");
    auto const v_t       = handle.get<double> (v_offset, v_readsize, "vertex/t");

    //// Some arithmetic to see how much reading is necessary for particles
    size_t p_offset   = p_start.front();
    size_t p_readsize = p_start.back() - p_start.front() + p_n.back();

    auto const p_id       = handle.get<int>    (p_offset, p_readsize, "particle/pid");    
    auto const p_status   = handle.get<int>    (p_offset, p_readsize, "particle/status");
    auto const p_px       = handle.get<double> (p_offset, p_readsize, "particle/px");
    auto const p_py       = handle.get<double> (p_offset, p_readsize, "particle/py");
    auto const p_pz       = handle.get<double> (p_offset, p_readsize, "particle/px");
    auto const p_e        = handle.get<double> (p_offset, p_readsize, "particle/E");
    auto const p_m        = handle.get<double> (p_offset, p_readsize, "particle/m");
    auto const p_vtx_end  = handle.get<int>    (p_offset, p_readsize, "particle/vtx_end");
    auto const p_vtx_prod = handle.get<int>    (p_offset, p_readsize, "particle/vtx_prod");

    auto const wgt = handle.weights(first_event, n_events); 


    size_t idxV(0), idxP(0);
    std::vector<GenVertexPtr> V;
    std::unordered_map<int,int> vmap;
        
    std::vector<GenEvent> EVENTS;
    EVENTS.reserve(n_events);
        //
        for (size_t nev=0; nev<n_events; ++nev)
        {
          for (unsigned int iv=idxV; iv < (v_n[nev]+idxV); ++iv)
          {
            GenVertexPtr v = std::make_shared<GenVertex>();
            v->set_position({v_x[iv], v_y[iv], v_z[iv], v_t[iv]});
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

            //Print::line(p);
          }

          idxV += v_n[nev];
          idxP += p_n[nev];
          
          GenEvent _evt(Units::GEV,Units::MM);
          _evt.weights() = wgt[nev];
          //_evt.run_info()->set_weight_names(wnames);
          _evt.set_event_number(nev);
          for (auto v: V) _evt.add_vertex(v);
          _evt.add_attribute("cycles", std::make_shared<IntAttribute>(1));
          //Print::listing(_evt);
          //_evt.weights().push_back(wgt[nev]);
          V.clear();
          vmap.clear();
          EVENTS.push_back(_evt);
        }
        return EVENTS;
    }
  //std::cout << "Processing loop took: " << t1-t0 << "\n";
//}

////Handle::Handle(std::string const & fname, std::vector<size_t> const & cols) : cols_(cols), file_(File(fname, File::ReadOnly))
////{
    //////this->file_ =File(fname, File::ReadOnly);
    //////this->cols_ = cols;
  //////file(argv[1], File::ReadOnly);//, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
    
////}
////


