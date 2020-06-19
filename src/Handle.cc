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
    evt.clear(); // Absolutely vital to prevent leakage
    evt.reserve(this->p_n[ievent], this->v_n[ievent]);


    // TODO figure out best way have dynamic size container holder shared ptrs.
    std::vector<GenVertexPtr> V;
    //V.reserve(v_n[ievent]);
    //std::cerr << ievent << " : " << v_n[ievent] << "\n";
    for (unsigned int iv=idxV; iv < (this->v_n[ievent]+idxV); ++iv)
        {
        GenVertexPtr v = std::make_shared<GenVertex>();
        v->set_position({this->v_x[iv], this->v_y[iv], this->v_z[iv], this->v_t[iv]});
        v->set_id(    this->v_id[iv]    );
        v->set_status(this->v_status[iv]);
        //V.push_back(v);
        V.push_back(std::move(v));
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
