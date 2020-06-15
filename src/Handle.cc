#include "Handle.h"
using namespace hepdf5;

void Handle::fillBuffers(size_t first_event, size_t n_events)
{
    this->idxV = 0;
    this->idxP = 0;
    read_vector<size_t>(first_event, n_events, this->P_START, this->p_start);
    read_vector<size_t>(first_event, n_events, this->V_START, this->v_start);
    read_vector<int>   (first_event, n_events, this->P_N    , this->p_n    );
    read_vector<int>   (first_event, n_events, this->V_N    , this->v_n    );

    // Some arithmetic to see how much reading is necessary for vertices
    size_t v_offset  = v_start.front();
    size_t v_readsize = v_start.back() - v_start.front() + v_n.back();

    read_vector<int>    (v_offset, v_readsize, this->V_ID    , this->v_id    );
    read_vector<int>    (v_offset, v_readsize, this->V_STATUS, this->v_status);
    read_vector<double> (v_offset, v_readsize, this->V_X     , this->v_x     );
    read_vector<double> (v_offset, v_readsize, this->V_Y     , this->v_y     );
    read_vector<double> (v_offset, v_readsize, this->V_Z     , this->v_z     );
    read_vector<double> (v_offset, v_readsize, this->V_T     , this->v_t     );

    //// Some arithmetic to see how much reading is necessary for particles
    size_t p_offset   = this->p_start.front();
    size_t p_readsize = this->p_start.back() - this->p_start.front() + p_n.back();

    read_vector<int>    (p_offset, p_readsize, this->P_ID      , this->p_id      );
    read_vector<int>    (p_offset, p_readsize, this->P_STATUS  , this->p_status  );
    read_vector<double> (p_offset, p_readsize, this->P_PX      , this->p_px      );
    read_vector<double> (p_offset, p_readsize, this->P_PY      , this->p_py      );
    read_vector<double> (p_offset, p_readsize, this->P_PZ      , this->p_pz      );
    read_vector<double> (p_offset, p_readsize, this->P_E       , this->p_e       );
    read_vector<double> (p_offset, p_readsize, this->P_M       , this->p_m       );
    read_vector<int>    (p_offset, p_readsize, this->P_VTX_END , this->p_vtx_end );
    read_vector<int>    (p_offset, p_readsize, this->P_VTX_PROD, this->p_vtx_prod);

    this->fillWeights(first_event, n_events);
}

void Handle::fillEvent(size_t ievent, GenEvent & evt)
{
    evt.clear();
    evt.reserve(this->p_n[ievent], this->v_n[ievent]);

    std::vector<GenVertexPtr> V;
    for (unsigned int iv=idxV; iv < (this->v_n[ievent]+idxV); ++iv)
        {
        GenVertexPtr v = std::make_shared<GenVertex>();
        v->set_position({this->v_x[iv], this->v_y[iv], this->v_z[iv], this->v_t[iv]});
        v->set_id(    this->v_id[iv]    );
        v->set_status(this->v_status[iv]);
        V.push_back(v);
    }

    for (unsigned int ip=idxP; ip < (this->p_n[ievent]+idxP); ++ip)
    {
        GenParticlePtr p = std::make_shared<GenParticle>(
            FourVector(this->p_px[ip], this->p_py[ip], this->p_pz[ip], this->p_e[ip]),
            this->p_id[ip], 
            this->p_status[ip]);
        //p->set_generated_mass(p_m[ip]);
        if (this->p_vtx_end[ip]  !=0) V[abs(this->p_vtx_end[ip])-1]->add_particle_in(p);
        if (this->p_vtx_prod[ip] !=0) V[abs(this->p_vtx_prod[ip])-1]->add_particle_out(p);
    }
    this->idxV += this->v_n[ievent];
    this->idxP += this->p_n[ievent];

  for (auto v: V) evt.add_vertex(v);
  evt.add_attribute("cycles", std::make_shared<IntAttribute>(1));
  evt.weights() = this->wgt[ievent];
}

//vector<GenEvent> hepdf5::readEvents(Handle & handle,  size_t first_event, size_t n_events)
//{
    ////std::cerr << first_event << "\n";
    //handle.fillBuffers(first_event, n_events);
    ////for (auto p : handle.p_start) std::cerr << " " << p;
    ////std::cerr <<"\n";
    //// Weights
    ////vector<vector<double> > e_weights;
    ////e_weights.reserve(n_events);
    //// NOTE this could be augmented to read only a central weight
    ////WGT.select({first_event, 0}, {n_events, wnames.size()}).read(e_weights);

    //// Event-wise quantities
    ////auto const p_start = read_vector<size_t>(first_event, n_events, handle.P_START);
    ////auto const v_start = read_vector<size_t>(first_event, n_events, handle.V_START);
    ////auto const p_n     = read_vector<int>   (first_event, n_events, handle.P_N);
    ////auto const v_n     = read_vector<int>   (first_event, n_events, handle.V_N);

    ////// Some arithmetic to see how much reading is necessary for vertices
    ////size_t v_offset  = v_start.front();
    ////size_t v_readsize = v_start.back() - v_start.front() + v_n.back();

    ////// Vertex quantities
    ////auto const v_id      = read_vector<int>    (v_offset, v_readsize, handle.V_ID    );
    ////auto const v_status  = read_vector<int>    (v_offset, v_readsize, handle.V_STATUS);
    ////auto const v_x       = read_vector<double> (v_offset, v_readsize, handle.V_X     );
    ////auto const v_y       = read_vector<double> (v_offset, v_readsize, handle.V_Y     );
    ////auto const v_z       = read_vector<double> (v_offset, v_readsize, handle.V_Z     );
    ////auto const v_t       = read_vector<double> (v_offset, v_readsize, handle.V_T     );

    //////// Some arithmetic to see how much reading is necessary for particles
    ////size_t p_offset   = handle.p_start.front();
    ////size_t p_readsize = handle.p_start.back() - handle.p_start.front() + p_n.back();

    ////auto const p_id       = read_vector<int>    (p_offset, p_readsize, handle.P_ID      );    
    ////auto const p_status   = read_vector<int>    (p_offset, p_readsize, handle.P_STATUS  );
    ////auto const p_px       = read_vector<double> (p_offset, p_readsize, handle.P_PX      );
    ////auto const p_py       = read_vector<double> (p_offset, p_readsize, handle.P_PY      );
    ////auto const p_pz       = read_vector<double> (p_offset, p_readsize, handle.P_PZ      );
    ////auto const p_e        = read_vector<double> (p_offset, p_readsize, handle.P_E       );
    ////auto const p_m        = read_vector<double> (p_offset, p_readsize, handle.P_M       );
    ////auto const p_vtx_end  = read_vector<int>    (p_offset, p_readsize, handle.P_VTX_END );
    ////auto const p_vtx_prod = read_vector<int>    (p_offset, p_readsize, handle.P_VTX_PROD);

    //auto const wgt = handle.weights(first_event, n_events); 


    //size_t idxV(0), idxP(0);
    //std::vector<GenVertexPtr> V;
    //V.reserve(std::max_element(handle.v_n.begin(), handle.v_n.end())[0]);

    //std::vector<GenEvent> EVENTS;
    //EVENTS.reserve(n_events);
        ////
    //for (size_t nev=0; nev<n_events; ++nev)
    //{
      //V.resize(handle.v_n[nev]);
      //int cv(0);
      //for (unsigned int iv=idxV; iv < (handle.v_n[nev]+idxV); ++iv)
      //{
        //GenVertexPtr v = std::make_shared<GenVertex>();
        //v->set_position({handle.v_x[iv], handle.v_y[iv], handle.v_z[iv], handle.v_t[iv]});
        //v->set_id(        handle.v_id[iv]);
        //v->set_status(handle.v_status[iv]);
        //V[cv]=v;
        //cv++;
      //}
     
      ////parts.reserve(p_n[nev]); 
      //for (unsigned int ip=idxP; ip < (handle.p_n[nev]+idxP); ++ip)
      //{
        //GenParticlePtr p = std::make_shared<GenParticle>(
            //FourVector(handle.p_px[ip], handle.p_py[ip], handle.p_pz[ip], handle.p_e[ip]),
            //handle.p_id[ip], handle.p_status[ip]);
        ////p->set_generated_mass(p_m[ip]);
        //if (handle.p_vtx_end[ip]  !=0) V[abs(handle.p_vtx_end[ip])-1]->add_particle_in(p);
        //if (handle.p_vtx_prod[ip] !=0) V[abs(handle.p_vtx_prod[ip])-1]->add_particle_out(p);
      //}

      //idxV += handle.v_n[nev];
      //idxP += handle.p_n[nev];
      
      ////std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();;
      ////run->set_weight_names(wnames);
      //GenEvent _evt(Units::GEV,Units::MM);//, run);
      //_evt.weights() = wgt[nev];
      ////_evt.run_info()->set_weight_names(wnames);
      //_evt.set_event_number(nev);
      //_evt.reserve(handle.p_n[nev], handle.v_n[nev]);
      //for (auto v: V) _evt.add_vertex(v);
      //_evt.add_attribute("cycles", std::make_shared<IntAttribute>(1));
      //EVENTS.push_back(_evt);
    //}
    //return std::move(EVENTS);
//}
